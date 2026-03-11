#' Diagnostic Plots
#'
#' Generic function for creating diagnostic plots.
#'
#' @param x An object for which diagnostic plots should be produced.
#' @param ... Additional arguments passed to methods.
#' @export
diagplot <- function(x, ...) UseMethod("diagplot")

# --- Internal Helper Functions (Not Exported) ---

generate_mu_grid <- function(yi, vi, M, X = NULL) {
  v_typ <- median(vi)

  if (is.null(X) || ncol(as.matrix(X)) == 0) {
    mu_rob <- median(yi)
    res <- yi - mu_rob
    var_total_rob <- mad(yi)^2
    theta_base <- mu_rob
  } else {
    X_mat <- as.matrix(X)
    fit <- lm(yi ~ X_mat, weights = 1 / vi)
    beta_hat <- coef(fit)[-1]
    beta_hat[is.na(beta_hat)] <- 0

    res_base <- yi - as.numeric(X_mat %*% beta_hat)
    mu_rob <- median(res_base)
    res <- res_base - mu_rob
    var_total_rob <- mad(res_base)^2
    theta_base <- mu_rob
  }

  tau2_rob <- max(0, var_total_rob - v_typ)
  shrinkage <- tau2_rob / (tau2_rob + vi)
  theta_tilde <- theta_base + shrinkage * res

  pad <- sqrt(tau2_rob)
  if (pad == 0) pad <- sqrt(v_typ)

  lower_bound <- min(theta_tilde) - 1.5 * pad
  upper_bound <- max(theta_tilde) + 1.5 * pad

  return(seq(lower_bound, upper_bound, length.out = M))
}

pgm_objective <- function(eta, yi, vi, X, mu_c, tau_c, lambda, DtD) {
  n <- length(yi)
  M <- length(mu_c)
  p <- ncol(X)
  alpha_0_index <- floor((M - 1) / 2) + 1

  alpha <- eta[1:(M - 1)]
  beta <- if (p > 0) eta[M:(M - 1 + p)] else NULL

  alpha_tilde <- numeric(M)
  alpha_tilde[-alpha_0_index] <- alpha
  alpha_tilde[alpha_0_index] <- 0

  max_alpha <- max(alpha_tilde)
  exp_alpha <- exp(alpha_tilde - max_alpha)
  w <- exp_alpha / sum(exp_alpha)

  mu_mat <- matrix(mu_c, nrow = n, ncol = M, byrow = TRUE)
  if (p > 0) mu_mat <- mu_mat + as.numeric(X %*% beta)

  sd2 <- vi + tau_c^2
  sd_mat <- matrix(sqrt(sd2), nrow = n, ncol = M)

  f_mat <- dnorm(yi, mean = mu_mat, sd = sd_mat)
  l_vec <- as.numeric(f_mat %*% w)
  l <- sum(log(l_vec))

  pen <- (lambda / 2) * as.numeric(t(alpha_tilde) %*% DtD %*% alpha_tilde)

  P_mat <- t(t(f_mat) * w) / l_vec
  s_tilde <- colSums(P_mat) - n * w

  grad_pen_tilde <- as.numeric(lambda * (DtD %*% alpha_tilde))
  grad_alpha <- (s_tilde - grad_pen_tilde)[-alpha_0_index]

  if (p > 0) {
    R_mat <- (yi - mu_mat) / sd2
    E_vec <- rowSums(P_mat * R_mat)
    grad_beta <- as.numeric(t(X) %*% E_vec)
    grad_obj <- c(grad_alpha, grad_beta)
  } else {
    grad_obj <- grad_alpha
  }

  J_tilde <- crossprod(P_mat) - (n * tcrossprod(w)) - diag(s_tilde)
  H_pen_tilde <- lambda * DtD
  H_alphaalpha <- -(J_tilde + H_pen_tilde)[-alpha_0_index, -alpha_0_index, drop = FALSE]

  if (p > 0) {
    W_vec <- rowSums(P_mat / sd2) - rowSums(P_mat * R_mat^2) + E_vec^2
    H_betabeta <- -crossprod(X, W_vec * X)

    H_alphabeta <- t(P_mat * (R_mat - E_vec)) %*% X
    H_alphabeta <- H_alphabeta[-alpha_0_index, , drop = FALSE]

    hess_obj <- rbind(cbind(H_alphaalpha, H_alphabeta), cbind(t(H_alphabeta), H_betabeta))
  } else {
    hess_obj <- H_alphaalpha
  }

  return(list(value = l - pen, gradient = grad_obj, hessian = hess_obj))
}

fit_pgm_path <- function(yi, vi, X, mu_c, tau_c, alpha_init, beta_init, lambda_grid, choose_lambda, DtD, ...) {
  n <- length(yi)
  M <- length(mu_c)
  p <- ncol(X)
  alpha_0_index <- floor((M - 1) / 2) + 1

  trust_args <- list(...)
  if (is.null(trust_args$rinit)) trust_args$rinit <- 1
  if (is.null(trust_args$rmax)) trust_args$rmax <- 100

  parinit <- if (p > 0) c(alpha_init, beta_init) else alpha_init

  results_list <- lapply(lambda_grid, function(lambda) {
    call_args <- c(
      list(
        objfun = pgm_objective, parinit = parinit, minimize = FALSE,
        yi = yi, vi = vi, X = X, mu_c = mu_c, tau_c = tau_c, lambda = lambda, DtD = DtD
      ),
      trust_args
    )

    opt <- tryCatch(do.call(trust::trust, call_args), error = function(e) NULL)

    if (is.null(opt) || !opt$converged) {
      return(list(ic = Inf, alpha = rep(NA, length(alpha_init)), beta = rep(NA, p), edf = NA, Jp_inv = NULL, J = NULL))
    }

    alpha <- opt$argument[1:(M - 1)]
    beta <- if (p > 0) opt$argument[M:(M - 1 + p)] else NULL

    Jp <- -opt$hessian
    Jp_inv <- tryCatch(chol2inv(chol(Jp)), error = function(e) NULL)

    if (is.null(Jp_inv)) {
      return(list(ic = Inf, alpha = rep(NA, length(alpha_init)), beta = rep(NA, p), edf = NA, Jp_inv = NULL, J = NULL))
    }

    alpha_tilde <- numeric(M)
    alpha_tilde[-alpha_0_index] <- alpha
    alpha_tilde[alpha_0_index] <- 0

    J <- Jp
    J[1:(M - 1), 1:(M - 1)] <- Jp[1:(M - 1), 1:(M - 1)] - (lambda * DtD)[-alpha_0_index, -alpha_0_index, drop = FALSE]

    edf <- sum(diag(J %*% Jp_inv))
    pen <- (lambda / 2) * as.numeric(t(alpha_tilde) %*% DtD %*% alpha_tilde)
    ic_val <- if (choose_lambda == "bic") -2 * (opt$value + pen) + log(n) * edf else -2 * (opt$value + pen) + 2 * edf

    return(list(ic = ic_val, alpha = alpha, beta = beta, edf = edf, Jp_inv = Jp_inv, J = J))
  })

  return(results_list)
}

calc_V <- function(yi, vi, X, group, mu_c, tau_c, alpha, beta, Jp_inv, J, V_type) {
  n <- length(yi)
  M <- length(mu_c)
  p <- ncol(X)
  alpha_0_index <- floor((M - 1) / 2) + 1

  if (V_type == "robust") {
    alpha_tilde <- numeric(M)
    alpha_tilde[-alpha_0_index] <- alpha
    alpha_tilde[alpha_0_index] <- 0

    max_alpha <- max(alpha_tilde)
    exp_alpha <- exp(alpha_tilde - max_alpha)
    w <- exp_alpha / sum(exp_alpha)

    mu_mat <- matrix(mu_c, nrow = n, ncol = M, byrow = TRUE)
    if (p > 0) mu_mat <- mu_mat + as.numeric(X %*% beta)

    sd2 <- vi + tau_c^2
    sd_mat <- matrix(sqrt(sd2), nrow = n, ncol = M)

    f_mat <- dnorm(yi, mean = mu_mat, sd = sd_mat)
    l_vec <- as.numeric(f_mat %*% w)
    P_mat <- t(t(f_mat) * w) / l_vec

    S_alpha <- P_mat[, -alpha_0_index, drop = FALSE] - matrix(w[-alpha_0_index], nrow = n, ncol = M - 1, byrow = TRUE)

    if (p > 0) {
      R_mat <- (yi - mu_mat) / sd2
      E_vec <- rowSums(P_mat * R_mat)
      S_beta <- X * E_vec
      S <- cbind(S_alpha, S_beta)
    } else {
      S <- S_alpha
    }

    S_group <- rowsum(S, group)
    B <- crossprod(S_group)

    return(Jp_inv %*% B %*% Jp_inv)
  } else if (V_type == "frequentist") {
    return(Jp_inv %*% J %*% Jp_inv)
  } else {
    return(Jp_inv)
  }
}

# --- Exported Functions & Methods ---

#' Fit a Penalized Gaussian Mixture Model
#'
#' Fits a penalized Gaussian mixture model for meta-analysis.
#'
#' @param yi Vector of effect sizes.
#' @param vi Vector of sampling variances.
#' @param X Optional design matrix for meta-regression.
#' @param group Optional grouping variable.
#' @param M Number of mixture components. Defaults to 40.
#' @param penalty_order Order of the penalty matrix. Defaults to 3.
#' @param lambda_grid Grid of penalty values to evaluate.
#' @param V_type Type of variance estimator ("robust", "frequentist", or "bayesian").
#' @param choose_lambda Information criterion for lambda selection ("aic" or "bic").
#' @param alpha_init Optional initial values for alpha.
#' @param beta_init Optional initial values for beta.
#' @param mu_grid Optional grid of means for the mixture components.
#' @param ... Additional arguments passed to the trust optimizer.
#'
#' @return An object of class \code{res.pgm}.
#'
#' @importFrom trust trust
#' @importFrom stats dnorm pnorm qnorm median mad lm coef uniroot
#' @importFrom graphics abline axis contour image legend lines matplot points polygon rect par plot hist
#' @importFrom grDevices adjustcolor hcl.colors
#' @export
rma.pgm <- function(yi, vi, X = NULL, group = NULL, M = 40, penalty_order = 3,
                    lambda_grid = exp(seq(-5, 5, by = 1)),
                    V_type = "robust", choose_lambda = "aic",
                    alpha_init = NULL, beta_init = NULL, mu_grid = NULL, ...) {
  if (!(choose_lambda %in% c("aic", "bic"))) stop("choose_lambda must be 'aic' or 'bic'.")
  if (!(V_type %in% c("frequentist", "bayesian", "robust"))) stop("V_type must be 'frequentist', 'bayesian', or 'robust'.")

  if (!is.null(X)) {
    X_mat <- as.matrix(X)
    p <- ncol(X_mat)
    if (is.null(colnames(X_mat))) colnames(X_mat) <- paste0("X", 1:p)
  } else {
    X_mat <- matrix(0, nrow = length(yi), ncol = 0)
    p <- 0
  }

  if (is.null(group)) group <- seq_along(yi)
  mu_c <- if (is.null(mu_grid)) generate_mu_grid(yi, vi, M, X = X_mat) else mu_grid

  M <- length(mu_c)
  tau_c <- (2 / 3) * (mu_c[2] - mu_c[1])
  alpha_0_index <- floor((M - 1) / 2) + 1

  D <- diff(diag(M), differences = penalty_order)
  DtD <- t(D) %*% D

  if (is.null(alpha_init)) alpha_init <- rep(0, M - 1)

  if (p > 0 && is.null(beta_init)) {
    fit_init <- lm(yi ~ X_mat, weights = 1 / vi)
    beta_init <- unname(coef(fit_init)[-1])
    beta_init[is.na(beta_init)] <- 0
  }

  path_results <- fit_pgm_path(yi, vi, X_mat, mu_c, tau_c, alpha_init, beta_init, lambda_grid, choose_lambda, DtD, ...)
  ics <- sapply(path_results, function(x) x$ic)
  best_idx <- which.min(ics)
  best <- path_results[[best_idx]]

  if (is.null(best) || is.infinite(best$ic)) stop("Optimization failed for all lambda values.\n")

  V_eta <- calc_V(yi, vi, X_mat, group, mu_c, tau_c, best$alpha, best$beta, best$Jp_inv, best$J, V_type)

  alpha_tilde <- numeric(M)
  alpha_tilde[-alpha_0_index] <- best$alpha
  alpha_tilde[alpha_0_index] <- 0

  max_alpha <- max(alpha_tilde)
  w <- exp(alpha_tilde - max_alpha) / sum(exp(alpha_tilde - max_alpha))

  out <- list(
    best_model = list(alpha = best$alpha, beta = best$beta, w = w, V_eta = V_eta, ic = best$ic, edf = best$edf),
    path_data = list(lambda = lambda_grid, ic = ics, edf = sapply(path_results, function(x) x$edf), alphas = do.call(rbind, lapply(path_results, function(x) x$alpha))),
    data = list(yi = yi, vi = vi, X = X_mat, group = group, V_type = V_type, choose_lambda = choose_lambda, mu_c = mu_c, tau_c = tau_c, alpha_0_index = alpha_0_index),
    p = p
  )
  class(out) <- "res.pgm"

  return(out)
}

#' Summarize a PGM Model
#'
#' @param object An object of class \code{res.pgm}.
#' @param level Confidence level. Defaults to 0.95.
#' @param less_than Optional numeric vector to calculate P(theta < val).
#' @param ... Additional arguments.
#' @export
summary.res.pgm <- function(object, level = 0.95, less_than = NULL, ...) {
  mu_c <- object$data$mu_c
  tau_c <- object$data$tau_c
  alpha_0_index <- object$data$alpha_0_index
  M <- length(mu_c)
  p <- object$p

  w <- object$best_model$w
  V_alpha <- object$best_model$V_eta[1:(M - 1), 1:(M - 1), drop = FALSE]
  z_crit <- qnorm(1 - (1 - level) / 2)

  mu_est <- sum(w * mu_c)
  grad_mu <- as.numeric(w * (mu_c - mu_est))[-alpha_0_index]
  var_mu <- as.numeric(t(grad_mu) %*% V_alpha %*% grad_mu)

  se_mu <- sqrt(max(0, var_mu))
  zval_mu <- ifelse(se_mu > 0, mu_est / se_mu, NA)
  pval_mu <- ifelse(!is.na(zval_mu), 2 * (1 - pnorm(abs(zval_mu))), NA)
  ci_mu <- c(mu_est - z_crit * se_mu, mu_est + z_crit * se_mu)

  res_mat <- matrix(c(mu_est, se_mu, zval_mu, pval_mu, ci_mu[1], ci_mu[2]), nrow = 1)
  rownames(res_mat) <- "Mean (\u03bc)"

  if (!is.null(less_than)) {
    prob_list <- lapply(less_than, function(val) {
      prob_val <- sum(w * pnorm(val, mean = mu_c, sd = tau_c))
      grad_p <- ((pnorm(val, mu_c, tau_c) - prob_val) * w)[-alpha_0_index]
      var_p <- as.numeric(t(grad_p) %*% V_alpha %*% grad_p)

      se_p <- sqrt(max(0, var_p))
      zval_p <- ifelse(se_p > 0, prob_val / se_p, NA)
      pval_p <- ifelse(!is.na(zval_p), 2 * (1 - pnorm(abs(zval_p))), NA)
      ci_p <- c(max(0, prob_val - z_crit * se_p), min(1, prob_val + z_crit * se_p))

      c(prob_val, se_p, zval_p, pval_p, ci_p[1], ci_p[2])
    })

    prob_stats <- do.call(rbind, prob_list)
    rownames(prob_stats) <- paste0("P(\u03b8 < ", less_than, ")")
    res_mat <- rbind(res_mat, prob_stats)
  }

  colnames(res_mat) <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")

  if (p > 0) {
    beta_est <- object$best_model$beta
    V_beta <- object$best_model$V_eta[(M):(M - 1 + p), (M):(M - 1 + p), drop = FALSE]
    se_beta <- sqrt(pmax(0, diag(V_beta)))
    zval_beta <- ifelse(se_beta > 0, beta_est / se_beta, NA)
    pval_beta <- ifelse(!is.na(zval_beta), 2 * (1 - pnorm(abs(zval_beta))), NA)
    ci_lb_beta <- beta_est - z_crit * se_beta
    ci_ub_beta <- beta_est + z_crit * se_beta

    beta_mat <- cbind(beta_est, se_beta, zval_beta, pval_beta, ci_lb_beta, ci_ub_beta)
    rownames(beta_mat) <- colnames(object$data$X)
    colnames(beta_mat) <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
  } else {
    beta_mat <- NULL
  }

  out <- list(res_mat = res_mat, beta_mat = beta_mat, n_obs = length(object$data$yi), n_groups = length(unique(object$data$group)), ic_name = toupper(object$data$choose_lambda), ic_val = object$best_model$ic, edf = object$best_model$edf, p = p)
  class(out) <- "summary.res.pgm"

  return(out)
}

#' Print PGM Summary
#'
#' @param x An object of class \code{summary.res.pgm}.
#' @param ... Additional arguments.
#' @export
print.summary.res.pgm <- function(x, ...) {
  cat("\n======================================================\n")
  cat("          Penalized Gaussian Mixture Results          \n")
  cat("======================================================\n")
  cat(sprintf("Observations : %d \nGroups       : %d\n", x$n_obs, x$n_groups))
  cat(sprintf("Model Fit    : %s = %.4f | EDF = %.4f\n", x$ic_name, x$ic_val, x$edf))

  if (x$p == 0) {
    cat("\n--- Parameter Estimates ---\n")
  } else {
    zero_str <- if (x$p == 1) "0" else paste(rep("0", x$p), collapse = ", ")
    cat(sprintf("\n--- Parameter Estimates (X = (%s)) ---\n", zero_str))
  }

  print(round(x$res_mat, digits = 4))

  if (!is.null(x$beta_mat)) {
    cat("\n--- Covariate Coefficients ---\n")
    print(round(x$beta_mat, digits = 4))
  }

  cat("\n")
  invisible(x)
}

#' Predict Method for PGM Model
#'
#' @param object An object of class \code{res.pgm}.
#' @param level Prediction interval level. Defaults to 0.95.
#' @param newmods Optional new moderator values for prediction.
#' @param ... Additional arguments.
#' @export
predict.res.pgm <- function(object, level = 0.95, newmods = NULL, ...) {
  mu_c <- object$data$mu_c
  tau_c <- object$data$tau_c
  w <- object$best_model$w
  p <- object$p

  if (p > 0) {
    if (is.null(newmods)) newmods <- matrix(0, nrow = 1, ncol = p)
    if (is.null(dim(newmods))) newmods <- matrix(newmods, nrow = 1)
    if (ncol(newmods) != p) stop(paste("newmods must have", p, "columns."))
    shifts <- as.numeric(newmods %*% object$best_model$beta)
  } else {
    newmods <- matrix(NA, nrow = 1, ncol = 0)
    shifts <- c(0)
  }

  res_list <- list()
  for (k in seq_along(shifts)) {
    shift <- shifts[k]
    cdf_fun <- function(q) sum(w * pnorm(q, mu_c + shift, tau_c))
    search_range <- c(min(mu_c + shift) - 5 * tau_c, max(mu_c + shift) + 5 * tau_c)

    pi_l <- tryCatch(uniroot(function(q) cdf_fun(q) - ((1 - level) / 2), interval = search_range, extendInt = "yes")$root, error = function(e) NA)
    pi_u <- tryCatch(uniroot(function(q) cdf_fun(q) - (1 - (1 - level) / 2), interval = search_range, extendInt = "yes")$root, error = function(e) NA)

    res_list[[k]] <- c(pi.lb = pi_l, pi.ub = pi_u)
  }

  res_mat <- do.call(rbind, res_list)

  if (p == 0) {
    rownames(res_mat) <- "Mean"
  } else {
    rownames(res_mat) <- apply(newmods, 1, function(row) {
      if (p == 1) paste0("X = ", round(row[1], 4)) else paste0("X = (", paste(round(row, 4), collapse = ", "), ")")
    })
  }

  cat(sprintf("\n--- %g%% Prediction Intervals (PI) ---\n", level * 100))
  print(round(res_mat, digits = 4))
  cat("\n")

  invisible(res_mat)
}

#' Plot Method for PGM Model
#'
#' @param x An object of class \code{res.pgm}.
#' @param y Ignored.
#' @param level Confidence level. Defaults to 0.95.
#' @param newmods Optional new moderator values.
#' @param addhist Logical indicating whether to add a histogram of yi.
#' @param ... Additional arguments.
#' @export
plot.res.pgm <- function(x, y = NULL, level = 0.95, newmods = NULL, addhist = TRUE, ...) {
  mu_c <- x$data$mu_c
  tau_c <- x$data$tau_c
  alpha_0_index <- x$data$alpha_0_index
  M <- length(mu_c)
  p <- x$p

  if (p > 0) {
    if (is.null(newmods)) newmods <- matrix(0, nrow = 1, ncol = p)
    if (is.null(dim(newmods))) newmods <- matrix(newmods, nrow = 1)
    if (ncol(newmods) != p) stop(paste("newmods must have", p, "columns."))
    shifts <- as.numeric(newmods %*% x$best_model$beta)
  } else {
    newmods <- matrix(NA, nrow = 1, ncol = 0)
    shifts <- c(0)
  }

  w <- x$best_model$w
  V_eta <- x$best_model$V_eta
  V_alpha <- V_eta[1:(M - 1), 1:(M - 1), drop = FALSE]
  z_crit <- qnorm(1 - (1 - level) / 2)

  colors <- c("firebrick", "steelblue", "forestgreen", "darkorange", "purple", "cyan4")
  if (length(shifts) > length(colors)) colors <- rep(colors, ceiling(length(shifts) / length(colors)))

  plot_data <- list()
  global_x_min <- Inf
  global_x_max <- -Inf
  global_y_max <- -Inf

  for (k in seq_along(shifts)) {
    shift <- shifts[k]
    mu_shifted <- mu_c + shift
    x_grid <- seq(min(mu_shifted) - 5 * tau_c, max(mu_shifted) + 5 * tau_c, length.out = 500)
    n_grid <- length(x_grid)
    global_x_min <- min(global_x_min, min(x_grid))
    global_x_max <- max(global_x_max, max(x_grid))

    f_mat_grid <- matrix(dnorm(x_grid, mean = rep(mu_shifted, each = n_grid), sd = tau_c), nrow = n_grid, ncol = M)
    pdf_vals <- as.numeric(f_mat_grid %*% w)

    grad_f_alpha <- t(t(f_mat_grid - pdf_vals) * w)[, -alpha_0_index, drop = FALSE]
    if (p > 0) {
      R_mat_grid <- (x_grid - matrix(mu_shifted, nrow = n_grid, ncol = M, byrow = TRUE)) / (tau_c^2)
      grad_f_beta <- rowSums((f_mat_grid * matrix(w, nrow = n_grid, ncol = M, byrow = TRUE)) * R_mat_grid) %*% t(newmods[k, ])
      grad_f_full <- cbind(grad_f_alpha, grad_f_beta)
      var_f <- rowSums((grad_f_full %*% V_eta) * grad_f_full)
    } else {
      var_f <- rowSums((grad_f_alpha %*% V_alpha) * grad_f_alpha)
    }

    se_pdf <- sqrt(pmax(0, var_f))
    y_pgm_lower <- pmax(0, pdf_vals - z_crit * se_pdf)
    y_pgm_upper <- pdf_vals + z_crit * se_pdf

    global_y_max <- max(global_y_max, max(y_pgm_upper, na.rm = TRUE))
    plot_data[[k]] <- list(x = x_grid, pdf = pdf_vals, lower = y_pgm_lower, upper = y_pgm_upper)
  }

  yi <- x$data$yi
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(1, 1), family = "serif", mar = c(5, 5, 2, 2), cex.axis = 1, cex.lab = 1.1)

  if (addhist) {
    h <- hist(yi, plot = FALSE, breaks = 30)
    y_max_final <- max(c(h$density, global_y_max), na.rm = TRUE) * 1.25
  } else {
    y_max_final <- global_y_max * 1.25
  }

  plot(NULL, xlim = c(global_x_min, global_x_max), ylim = c(0, y_max_final), main = "PGM Estimated Density", xlab = "Effect Size", ylab = "Density", axes = FALSE, frame.plot = FALSE)
  axis(1, tick = TRUE)
  axis(2, las = 1, tick = TRUE)

  if (addhist) rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$density, col = "gray92", border = "gray75", lwd = 1)

  for (k in seq_along(shifts)) {
    polygon(c(plot_data[[k]]$x, rev(plot_data[[k]]$x)), c(plot_data[[k]]$lower, rev(plot_data[[k]]$upper)), col = adjustcolor(colors[k], alpha.f = 0.2), border = NA)
    lines(plot_data[[k]]$x, plot_data[[k]]$pdf, lwd = 2, col = colors[k])
  }

  if (p > 0) {
    legend_labels <- apply(newmods, 1, function(row) {
      if (p == 1) paste0("X = ", round(row[1], 4)) else paste0("X = (", paste(round(row, 4), collapse = ", "), ")")
    })
    legend("topright", legend = legend_labels, col = colors[1:length(shifts)], lwd = 2, bty = "n")
  }
}

#' Diagnostic Plots for PGM Model
#'
#' @param x An object of class \code{res.pgm}.
#' @param ... Additional arguments.
#' @export
diagplot.res.pgm <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(3, 1), mar = c(4, 5, 2, 2), family = "serif")
  log_lam <- log(x$path_data$lambda)
  best_idx <- which.min(x$path_data$ic)
  best_log_lam <- log_lam[best_idx]
  ic_label <- toupper(x$data$choose_lambda)

  plot(log_lam, x$path_data$ic, type = "b", pch = 16, col = "darkblue", xlab = expression(log(lambda)), ylab = ic_label, main = paste(ic_label, "Profile"))
  abline(v = best_log_lam, col = "firebrick", lty = 2, lwd = 2)

  plot(log_lam, x$path_data$edf, type = "b", pch = 16, col = "forestgreen", xlab = expression(log(lambda)), ylab = "EDF", main = "EDF Profile")
  abline(v = best_log_lam, col = "firebrick", lty = 2, lwd = 2)

  matplot(log_lam, x$path_data$alphas, type = "l", lty = 1, xlab = expression(log(lambda)), ylab = expression(alpha), main = "Profile Alpha Path")
  abline(v = best_log_lam, col = "firebrick", lty = 2, lwd = 2)
}

# --- Internal Conditional Helper Functions (Not Exported) ---

pgm_cond_objective <- function(eta, yi, vi, zi, mu_c, tau_c, lambda_alpha, lambda_gamma, DtD) {
  n <- length(yi)
  M <- length(mu_c)
  alpha_0_index <- floor((M - 1) / 2) + 1

  alpha <- eta[1:(M - 1)]
  gamma <- eta[M:(2 * M - 2)]

  alpha_tilde <- numeric(M)
  alpha_tilde[-alpha_0_index] <- alpha
  gamma_tilde <- numeric(M)
  gamma_tilde[-alpha_0_index] <- gamma

  eta_mat <- matrix(alpha_tilde, nrow = n, ncol = M, byrow = TRUE) + as.numeric(zi) %*% t(gamma_tilde)

  max_eta <- apply(eta_mat, 1, max)
  exp_eta <- exp(eta_mat - max_eta)
  W_mat <- exp_eta / rowSums(exp_eta)

  sd_mat <- matrix(sqrt(vi + tau_c^2), nrow = n, ncol = M)
  mu_mat <- matrix(mu_c, nrow = n, ncol = M, byrow = TRUE)
  f_mat <- dnorm(yi, mean = mu_mat, sd = sd_mat)

  l_vec <- rowSums(f_mat * W_mat)
  l <- sum(log(l_vec))

  pen_alpha <- (lambda_alpha / 2) * as.numeric(t(alpha_tilde) %*% DtD %*% alpha_tilde)
  pen_gamma <- (lambda_gamma / 2) * as.numeric(t(gamma_tilde) %*% DtD %*% gamma_tilde)
  pen <- pen_alpha + pen_gamma

  P_mat <- (f_mat * W_mat) / l_vec
  diff_mat <- P_mat - W_mat

  s_tilde_alpha <- colSums(diff_mat)
  s_tilde_gamma <- colSums(diff_mat * as.numeric(zi))

  grad_pen_alpha_tilde <- as.numeric(lambda_alpha * (DtD %*% alpha_tilde))
  grad_pen_gamma_tilde <- as.numeric(lambda_gamma * (DtD %*% gamma_tilde))

  grad_alpha <- (s_tilde_alpha - grad_pen_alpha_tilde)[-alpha_0_index]
  grad_gamma <- (s_tilde_gamma - grad_pen_gamma_tilde)[-alpha_0_index]
  grad_obj <- c(grad_alpha, grad_gamma)

  P_z <- P_mat * as.numeric(zi)
  W_z <- W_mat * as.numeric(zi)

  H_alphaalpha <- diag(s_tilde_alpha) - (crossprod(P_mat) - crossprod(W_mat)) - lambda_alpha * DtD
  H_alphagamma <- diag(s_tilde_gamma) - (crossprod(P_mat, P_z) - crossprod(W_mat, W_z))
  H_gammagamma <- diag(colSums(diff_mat * as.numeric(zi)^2)) - (crossprod(P_z) - crossprod(W_z)) - lambda_gamma * DtD

  idx <- -alpha_0_index
  hess_obj <- rbind(
    cbind(H_alphaalpha[idx, idx, drop = FALSE], H_alphagamma[idx, idx, drop = FALSE]),
    cbind(t(H_alphagamma[idx, idx, drop = FALSE]), H_gammagamma[idx, idx, drop = FALSE])
  )

  return(list(value = l - pen, gradient = grad_obj, hessian = hess_obj))
}

fit_pgm_cond_path <- function(yi, vi, zi, mu_c, tau_c, eta_init, lambda_alpha_grid, lambda_gamma_grid, choose_lambda, DtD, ...) {
  n <- length(yi)
  M <- length(mu_c)
  alpha_0_index <- floor((M - 1) / 2) + 1

  trust_args <- list(...)
  if (is.null(trust_args$rinit)) trust_args$rinit <- 1
  if (is.null(trust_args$rmax)) trust_args$rmax <- 100

  grid <- expand.grid(lam_a = lambda_alpha_grid, lam_g = lambda_gamma_grid)

  results_list <- lapply(1:nrow(grid), function(i) {
    lam_a <- grid$lam_a[i]
    lam_g <- grid$lam_g[i]
    call_args <- c(list(objfun = pgm_cond_objective, parinit = eta_init, minimize = FALSE, yi = yi, vi = vi, zi = zi, mu_c = mu_c, tau_c = tau_c, lambda_alpha = lam_a, lambda_gamma = lam_g, DtD = DtD), trust_args)
    opt <- tryCatch(do.call(trust::trust, call_args), error = function(e) NULL)

    if (is.null(opt) || !opt$converged) {
      return(list(ic = Inf, alpha = rep(NA, M - 1), gamma = rep(NA, M - 1), edf = NA, Jp_inv = NULL, J = NULL, lambda_alpha = lam_a, lambda_gamma = lam_g))
    }

    alpha <- opt$argument[1:(M - 1)]
    gamma <- opt$argument[M:(2 * M - 2)]

    Jp <- -opt$hessian
    Jp_inv <- tryCatch(chol2inv(chol(Jp)), error = function(e) NULL)
    if (is.null(Jp_inv)) {
      return(list(ic = Inf, alpha = rep(NA, M - 1), gamma = rep(NA, M - 1), edf = NA, Jp_inv = NULL, J = NULL, lambda_alpha = lam_a, lambda_gamma = lam_g))
    }

    alpha_tilde <- numeric(M)
    alpha_tilde[-alpha_0_index] <- alpha
    gamma_tilde <- numeric(M)
    gamma_tilde[-alpha_0_index] <- gamma

    J <- Jp
    J[1:(M - 1), 1:(M - 1)] <- Jp[1:(M - 1), 1:(M - 1)] - (lam_a * DtD)[-alpha_0_index, -alpha_0_index, drop = FALSE]
    J[M:(2 * M - 2), M:(2 * M - 2)] <- Jp[M:(2 * M - 2), M:(2 * M - 2)] - (lam_g * DtD)[-alpha_0_index, -alpha_0_index, drop = FALSE]

    edf <- sum(diag(J %*% Jp_inv))
    pen <- (lam_a / 2) * as.numeric(t(alpha_tilde) %*% DtD %*% alpha_tilde) + (lam_g / 2) * as.numeric(t(gamma_tilde) %*% DtD %*% gamma_tilde)
    ic_val <- if (choose_lambda == "bic") -2 * (opt$value + pen) + log(n) * edf else -2 * (opt$value + pen) + 2 * edf

    return(list(ic = ic_val, alpha = alpha, gamma = gamma, edf = edf, Jp_inv = Jp_inv, J = J, lambda_alpha = lam_a, lambda_gamma = lam_g))
  })

  return(results_list)
}

calc_V_cond <- function(yi, vi, zi, group, mu_c, tau_c, alpha, gamma, Jp_inv, J, V_type) {
  n <- length(yi)
  M <- length(mu_c)
  alpha_0_index <- floor((M - 1) / 2) + 1

  if (V_type == "robust") {
    alpha_tilde <- numeric(M)
    alpha_tilde[-alpha_0_index] <- alpha
    gamma_tilde <- numeric(M)
    gamma_tilde[-alpha_0_index] <- gamma

    eta_mat <- matrix(alpha_tilde, nrow = n, ncol = M, byrow = TRUE) + as.numeric(zi) %*% t(gamma_tilde)

    max_eta <- apply(eta_mat, 1, max)
    exp_eta <- exp(eta_mat - max_eta)
    W_mat <- exp_eta / rowSums(exp_eta)

    sd_mat <- matrix(sqrt(vi + tau_c^2), nrow = n, ncol = M)
    mu_mat <- matrix(mu_c, nrow = n, ncol = M, byrow = TRUE)
    f_mat <- dnorm(yi, mean = mu_mat, sd = sd_mat)

    l_vec <- rowSums(f_mat * W_mat)
    P_mat <- (f_mat * W_mat) / l_vec
    diff_mat <- P_mat - W_mat

    S_alpha <- diff_mat[, -alpha_0_index, drop = FALSE]
    S_gamma <- (diff_mat * as.numeric(zi))[, -alpha_0_index, drop = FALSE]
    S <- cbind(S_alpha, S_gamma)

    S_group <- rowsum(S, group)
    B <- crossprod(S_group)

    return(Jp_inv %*% B %*% Jp_inv)
  } else if (V_type == "frequentist") {
    return(Jp_inv %*% J %*% Jp_inv)
  } else {
    return(Jp_inv)
  }
}

# --- Exported Conditional Functions & Methods ---

#' Fit a Conditional Penalized Gaussian Mixture Model
#'
#' Fits a conditional penalized Gaussian mixture model for meta-analysis.
#'
#' @param yi Vector of effect sizes.
#' @param vi Vector of sampling variances.
#' @param zi Vector of conditional variable values.
#' @param group Optional grouping variable.
#' @param M Number of mixture components. Defaults to 40.
#' @param penalty_order Order of the penalty matrix. Defaults to 3.
#' @param lambda_alpha_grid Grid of penalty values for alpha.
#' @param lambda_gamma_grid Grid of penalty values for gamma.
#' @param V_type Type of variance estimator ("robust", "frequentist", or "bayesian").
#' @param choose_lambda Information criterion for lambda selection ("aic" or "bic").
#' @param alpha_init Optional initial values for alpha.
#' @param gamma_init Optional initial values for gamma.
#' @param mu_grid Optional grid of means for the mixture components.
#' @param ... Additional arguments passed to the trust optimizer.
#'
#' @return An object of class \code{res.pgm.cond}.
#' @export
rma.pgm.cond <- function(yi, vi, zi, group = NULL, M = 40, penalty_order = 3,
                         lambda_alpha_grid = exp(seq(-5, 5, by = 1)),
                         lambda_gamma_grid = exp(seq(-5, 5, by = 1)),
                         V_type = "robust", choose_lambda = "aic",
                         alpha_init = NULL, gamma_init = NULL, mu_grid = NULL, ...) {
  if (!(choose_lambda %in% c("aic", "bic"))) stop("choose_lambda must be 'aic' or 'bic'.")
  if (!(V_type %in% c("frequentist", "bayesian", "robust"))) stop("V_type must be 'frequentist', 'bayesian', or 'robust'.")
  if (is.null(group)) group <- seq_along(yi)

  mu_c <- if (is.null(mu_grid)) generate_mu_grid(yi, vi, M) else mu_grid

  M <- length(mu_c)
  tau_c <- (2 / 3) * (mu_c[2] - mu_c[1])
  alpha_0_index <- floor((M - 1) / 2) + 1

  D <- diff(diag(M), differences = penalty_order)
  DtD <- t(D) %*% D

  if (is.null(alpha_init)) alpha_init <- rep(0, M - 1)
  if (is.null(gamma_init)) gamma_init <- rep(0, M - 1)
  eta_init <- c(alpha_init, gamma_init)

  path_results <- fit_pgm_cond_path(yi, vi, zi, mu_c, tau_c, eta_init, lambda_alpha_grid, lambda_gamma_grid, choose_lambda, DtD, ...)
  ics <- sapply(path_results, function(x) x$ic)
  best_idx <- which.min(ics)
  best <- path_results[[best_idx]]

  if (is.null(best) || is.infinite(best$ic)) stop("Optimization failed for all lambda combinations.\n")

  V_eta <- calc_V_cond(yi, vi, zi, group, mu_c, tau_c, best$alpha, best$gamma, best$Jp_inv, best$J, V_type)

  out <- list(
    best_model = list(alpha = best$alpha, gamma = best$gamma, V_eta = V_eta, ic = best$ic, edf = best$edf, lambda_alpha = best$lambda_alpha, lambda_gamma = best$lambda_gamma),
    path_data = list(lambda_alpha = sapply(path_results, function(x) x$lambda_alpha), lambda_gamma = sapply(path_results, function(x) x$lambda_gamma), ic = ics, edf = sapply(path_results, function(x) x$edf), alphas = do.call(rbind, lapply(path_results, function(x) x$alpha)), gammas = do.call(rbind, lapply(path_results, function(x) x$gamma))),
    data = list(yi = yi, vi = vi, zi = zi, group = group, V_type = V_type, choose_lambda = choose_lambda, mu_c = mu_c, tau_c = tau_c, alpha_0_index = alpha_0_index)
  )
  class(out) <- "res.pgm.cond"

  return(out)
}

#' Summarize a Conditional PGM Model
#'
#' @param object An object of class \code{res.pgm.cond}.
#' @param level Confidence level. Defaults to 0.95.
#' @param less_than Optional numeric vector to calculate P(theta < val).
#' @param at Optional numeric vector of z values to evaluate the summary at.
#' @param compare Optional vector or matrix to compare specific z values.
#' @param ... Additional arguments.
#' @export
summary.res.pgm.cond <- function(object, level = 0.95, less_than = NULL, at = NULL, compare = NULL, ...) {
  mu_c <- object$data$mu_c
  tau_c <- object$data$tau_c
  alpha_0_index <- object$data$alpha_0_index
  M <- length(mu_c)
  V_eta <- object$best_model$V_eta
  z_crit <- qnorm(1 - (1 - level) / 2)

  alpha_tilde <- numeric(M)
  alpha_tilde[-alpha_0_index] <- object$best_model$alpha
  gamma_tilde <- numeric(M)
  gamma_tilde[-alpha_0_index] <- object$best_model$gamma

  at_vals <- if (is.null(at)) 0 else at
  all_res <- list()

  for (z_val in at_vals) {
    eta_z <- alpha_tilde + gamma_tilde * z_val

    max_eta_z <- max(eta_z)
    w_z <- exp(eta_z - max_eta_z) / sum(exp(eta_z - max_eta_z))

    mu_est <- sum(w_z * mu_c)
    J_w <- diag(w_z) - tcrossprod(w_z)
    J_mu_alpha <- as.numeric(crossprod(J_w, mu_c))[-alpha_0_index]
    J_mu_gamma <- z_val * J_mu_alpha
    grad_mu <- c(J_mu_alpha, J_mu_gamma)

    se_mu <- sqrt(max(0, as.numeric(t(grad_mu) %*% V_eta %*% grad_mu)))
    zval_mu <- ifelse(se_mu > 0, mu_est / se_mu, NA)
    pval_mu <- ifelse(!is.na(zval_mu), 2 * (1 - pnorm(abs(zval_mu))), NA)
    ci_mu <- c(mu_est - z_crit * se_mu, mu_est + z_crit * se_mu)

    res_mat <- matrix(c(mu_est, se_mu, zval_mu, pval_mu, ci_mu[1], ci_mu[2]), nrow = 1)
    rownames(res_mat) <- "Mean (\u03bc)"

    if (!is.null(less_than)) {
      prob_list <- lapply(less_than, function(val) {
        p_c <- pnorm(val, mean = mu_c, sd = tau_c)
        prob_val <- sum(w_z * p_c)
        J_p_alpha <- as.numeric(crossprod(J_w, p_c))[-alpha_0_index]
        J_p_gamma <- z_val * J_p_alpha
        grad_p <- c(J_p_alpha, J_p_gamma)

        se_p <- sqrt(max(0, as.numeric(t(grad_p) %*% V_eta %*% grad_p)))
        zval_p <- ifelse(se_p > 0, prob_val / se_p, NA)
        pval_p <- ifelse(!is.na(zval_p), 2 * (1 - pnorm(abs(zval_p))), NA)
        ci_p <- c(max(0, prob_val - z_crit * se_p), min(1, prob_val + z_crit * se_p))

        c(prob_val, se_p, zval_p, pval_p, ci_p[1], ci_p[2])
      })
      prob_stats <- do.call(rbind, prob_list)
      rownames(prob_stats) <- paste0("P(\u03b8 < ", less_than, ")")
      res_mat <- rbind(res_mat, prob_stats)
    }

    colnames(res_mat) <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")

    header <- paste("z =", z_val)
    all_res[[header]] <- res_mat
  }

  compare_res <- NULL
  if (!is.null(compare)) {
    if (is.vector(compare)) compare <- matrix(compare, nrow = 1)
    if (ncol(compare) != 2) stop("compare must be a vector of length 2 or a matrix with 2 columns.")

    compare_res <- list()
    for (i in 1:nrow(compare)) {
      z1 <- compare[i, 1]
      z2 <- compare[i, 2]

      eta_z1 <- alpha_tilde + gamma_tilde * z1
      w_z1 <- exp(eta_z1 - max(eta_z1)) / sum(exp(eta_z1 - max(eta_z1)))
      J_w1 <- diag(w_z1) - tcrossprod(w_z1)

      eta_z2 <- alpha_tilde + gamma_tilde * z2
      w_z2 <- exp(eta_z2 - max(eta_z2)) / sum(exp(eta_z2 - max(eta_z2)))
      J_w2 <- diag(w_z2) - tcrossprod(w_z2)

      mu_diff <- sum(w_z2 * mu_c) - sum(w_z1 * mu_c)

      J_mu_alpha1 <- as.numeric(crossprod(J_w1, mu_c))[-alpha_0_index]
      grad_mu1 <- c(J_mu_alpha1, z1 * J_mu_alpha1)

      J_mu_alpha2 <- as.numeric(crossprod(J_w2, mu_c))[-alpha_0_index]
      grad_mu2 <- c(J_mu_alpha2, z2 * J_mu_alpha2)

      grad_mu_diff <- grad_mu2 - grad_mu1

      se_mu_diff <- sqrt(max(0, as.numeric(t(grad_mu_diff) %*% V_eta %*% grad_mu_diff)))
      zval_mu_diff <- ifelse(se_mu_diff > 0, mu_diff / se_mu_diff, NA)
      pval_mu_diff <- ifelse(!is.na(zval_mu_diff), 2 * (1 - pnorm(abs(zval_mu_diff))), NA)
      ci_mu_diff <- c(mu_diff - z_crit * se_mu_diff, mu_diff + z_crit * se_mu_diff)

      res_mat_diff <- matrix(c(mu_diff, se_mu_diff, zval_mu_diff, pval_mu_diff, ci_mu_diff[1], ci_mu_diff[2]), nrow = 1)
      rownames(res_mat_diff) <- "Mean Diff (\u03bc)"

      if (!is.null(less_than)) {
        prob_list_diff <- lapply(less_than, function(val) {
          p_c <- pnorm(val, mean = mu_c, sd = tau_c)
          prob_val1 <- sum(w_z1 * p_c)
          prob_val2 <- sum(w_z2 * p_c)
          prob_diff <- prob_val2 - prob_val1

          J_p_alpha1 <- as.numeric(crossprod(J_w1, p_c))[-alpha_0_index]
          grad_p1 <- c(J_p_alpha1, z1 * J_p_alpha1)

          J_p_alpha2 <- as.numeric(crossprod(J_w2, p_c))[-alpha_0_index]
          grad_p2 <- c(J_p_alpha2, z2 * J_p_alpha2)

          grad_p_diff <- grad_p2 - grad_p1

          se_p_diff <- sqrt(max(0, as.numeric(t(grad_p_diff) %*% V_eta %*% grad_p_diff)))
          zval_p_diff <- ifelse(se_p_diff > 0, prob_diff / se_p_diff, NA)
          pval_p_diff <- ifelse(!is.na(zval_p_diff), 2 * (1 - pnorm(abs(zval_p_diff))), NA)

          ci_p_diff <- c(prob_diff - z_crit * se_p_diff, prob_diff + z_crit * se_p_diff)

          c(prob_diff, se_p_diff, zval_p_diff, pval_p_diff, ci_p_diff[1], ci_p_diff[2])
        })

        prob_stats_diff <- do.call(rbind, prob_list_diff)
        rownames(prob_stats_diff) <- paste0("\u0394 P(\u03b8 < ", less_than, ")")
        res_mat_diff <- rbind(res_mat_diff, prob_stats_diff)
      }

      colnames(res_mat_diff) <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")

      header_diff <- paste("z =", z2, "vs z =", z1)
      compare_res[[header_diff]] <- res_mat_diff
    }
  }

  out <- list(all_res = all_res, compare_res = compare_res, n_obs = length(object$data$yi), n_groups = length(unique(object$data$group)), ic_name = toupper(object$data$choose_lambda), ic_val = object$best_model$ic, edf = object$best_model$edf)
  class(out) <- "summary.res.pgm.cond"

  return(out)
}

#' Print Conditional PGM Summary
#'
#' @param x An object of class \code{summary.res.pgm.cond}.
#' @param ... Additional arguments.
#' @export
print.summary.res.pgm.cond <- function(x, ...) {
  cat("\n======================================================\n")
  cat("    Conditional Penalized Gaussian Mixture Results    \n")
  cat("======================================================\n")
  cat(sprintf("Observations : %d \nGroups       : %d\n", x$n_obs, x$n_groups))
  cat(sprintf("Model Fit    : %s = %.4f | EDF = %.4f\n", x$ic_name, x$ic_val, x$edf))

  for (header in names(x$all_res)) {
    cat(sprintf("\n--- Parameter Estimates (%s) ---\n", header))
    print(round(x$all_res[[header]], digits = 4))
  }

  if (!is.null(x$compare_res)) {
    for (header in names(x$compare_res)) {
      cat(sprintf("\n--- Comparison (%s) ---\n", header))
      print(round(x$compare_res[[header]], digits = 4))
    }
  }

  cat("\n")
  invisible(x)
}

#' Predict Method for Conditional PGM Model
#'
#' @param object An object of class \code{res.pgm.cond}.
#' @param level Prediction interval level. Defaults to 0.95.
#' @param at Optional numeric vector of z values to predict at.
#' @param ... Additional arguments.
#' @export
predict.res.pgm.cond <- function(object, level = 0.95, at = NULL, ...) {
  mu_c <- object$data$mu_c
  tau_c <- object$data$tau_c
  alpha_0_index <- object$data$alpha_0_index
  M <- length(mu_c)

  alpha_tilde <- numeric(M)
  alpha_tilde[-alpha_0_index] <- object$best_model$alpha
  gamma_tilde <- numeric(M)
  gamma_tilde[-alpha_0_index] <- object$best_model$gamma

  at_vals <- if (is.null(at)) 0 else at
  res_list <- list()

  for (k in seq_along(at_vals)) {
    z_val <- at_vals[k]
    eta_z <- alpha_tilde + gamma_tilde * z_val

    max_eta_z <- max(eta_z)
    w_z <- exp(eta_z - max_eta_z) / sum(exp(eta_z - max_eta_z))

    cdf_fun <- function(q) sum(w_z * pnorm(q, mu_c, tau_c))
    search_range <- c(min(mu_c) - 5 * tau_c, max(mu_c) + 5 * tau_c)

    pi_l <- tryCatch(uniroot(function(q) cdf_fun(q) - ((1 - level) / 2), interval = search_range, extendInt = "yes")$root, error = function(e) NA)
    pi_u <- tryCatch(uniroot(function(q) cdf_fun(q) - (1 - (1 - level) / 2), interval = search_range, extendInt = "yes")$root, error = function(e) NA)

    res_list[[k]] <- c(pi.lb = pi_l, pi.ub = pi_u)
  }

  res_mat <- do.call(rbind, res_list)
  rownames(res_mat) <- paste("z =", at_vals)

  cat(sprintf("\n--- %g%% Prediction Intervals (PI) ---\n", level * 100))
  print(round(res_mat, digits = 4))
  cat("\n")

  invisible(res_mat)
}

#' Plot Method for Conditional PGM Model
#'
#' @param x An object of class \code{res.pgm.cond}.
#' @param y Ignored.
#' @param level Confidence level. Defaults to 0.95.
#' @param at Optional numeric vector of z values to plot.
#' @param addhist Logical indicating whether to add a histogram.
#' @param ... Additional arguments.
#' @export
plot.res.pgm.cond <- function(x, y = NULL, level = 0.95, at = NULL, addhist = TRUE, ...) {
  mu_c <- x$data$mu_c
  tau_c <- x$data$tau_c
  alpha_0_index <- x$data$alpha_0_index
  M <- length(mu_c)
  V_eta <- x$best_model$V_eta
  z_crit <- qnorm(1 - (1 - level) / 2)

  alpha_tilde <- numeric(M)
  alpha_tilde[-alpha_0_index] <- x$best_model$alpha
  gamma_tilde <- numeric(M)
  gamma_tilde[-alpha_0_index] <- x$best_model$gamma

  at_vals <- if (is.null(at)) 0 else at
  colors <- c("firebrick", "steelblue", "forestgreen", "darkorange", "purple", "cyan4")
  if (length(at_vals) > length(colors)) colors <- rep(colors, ceiling(length(at_vals) / length(colors)))

  plot_data <- list()
  global_y_max <- -Inf
  x_grid <- seq(min(mu_c) - 5 * tau_c, max(mu_c) + 5 * tau_c, length.out = 500)
  n_grid <- length(x_grid)
  f_mat_grid <- matrix(dnorm(x_grid, mean = rep(mu_c, each = n_grid), sd = tau_c), nrow = n_grid, ncol = M)

  for (k in seq_along(at_vals)) {
    z_val <- at_vals[k]
    eta_z <- alpha_tilde + gamma_tilde * z_val

    max_eta_z <- max(eta_z)
    w_z <- exp(eta_z - max_eta_z) / sum(exp(eta_z - max_eta_z))

    pdf_vals <- as.numeric(f_mat_grid %*% w_z)
    J_w <- diag(w_z) - tcrossprod(w_z)
    grad_f_alpha <- t(J_w %*% t(f_mat_grid))[, -alpha_0_index, drop = FALSE]
    grad_f_gamma <- z_val * grad_f_alpha
    grad_f_full <- cbind(grad_f_alpha, grad_f_gamma)

    var_f <- rowSums((grad_f_full %*% V_eta) * grad_f_full)
    se_pdf <- sqrt(pmax(0, var_f))
    y_lower <- pmax(0, pdf_vals - z_crit * se_pdf)
    y_upper <- pdf_vals + z_crit * se_pdf

    global_y_max <- max(global_y_max, max(y_upper, na.rm = TRUE))
    plot_data[[k]] <- list(x = x_grid, pdf = pdf_vals, lower = y_lower, upper = y_upper)
  }

  yi <- x$data$yi
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(1, 1), family = "serif", mar = c(5, 5, 2, 2), cex.axis = 1, cex.lab = 1.1)

  if (addhist) {
    h <- hist(yi, plot = FALSE, breaks = 30)
    y_max_final <- max(c(h$density, global_y_max), na.rm = TRUE) * 1.25
  } else {
    y_max_final <- global_y_max * 1.25
  }

  plot(NULL, xlim = range(x_grid), ylim = c(0, y_max_final), main = "PGM Estimated Density", xlab = "Effect Size", ylab = "Density", axes = FALSE, frame.plot = FALSE)
  axis(1, tick = TRUE)
  axis(2, las = 1, tick = TRUE)

  if (addhist) rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$density, col = "gray92", border = "gray75", lwd = 1)

  for (k in seq_along(at_vals)) {
    polygon(c(plot_data[[k]]$x, rev(plot_data[[k]]$x)), c(plot_data[[k]]$lower, rev(plot_data[[k]]$upper)), col = adjustcolor(colors[k], alpha.f = 0.2), border = NA)
    lines(plot_data[[k]]$x, plot_data[[k]]$pdf, lwd = 2, col = colors[k])
  }

  legend_labels <- paste("z =", at_vals)
  legend("topright", legend = legend_labels, col = colors[1:length(at_vals)], lwd = 2, bty = "n")
}

#' Diagnostic Plots for Conditional PGM Model
#'
#' @param x An object of class \code{res.pgm.cond}.
#' @param ... Additional arguments.
#' @export
diagplot.res.pgm.cond <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  lam_a_vals <- unique(x$path_data$lambda_alpha)
  lam_g_vals <- unique(x$path_data$lambda_gamma)
  log_lam_a <- log(lam_a_vals)
  log_lam_g <- log(lam_g_vals)

  ic_mat <- matrix(x$path_data$ic, nrow = length(lam_a_vals), ncol = length(lam_g_vals))
  edf_mat <- matrix(x$path_data$edf, nrow = length(lam_a_vals), ncol = length(lam_g_vals))

  best_lam_a <- log(x$best_model$lambda_alpha)
  best_lam_g <- log(x$best_model$lambda_gamma)
  ic_label <- toupper(x$data$choose_lambda)

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), family = "serif")

  image(log_lam_a, log_lam_g, ic_mat, col = hcl.colors(50, "YlOrRd", rev = TRUE), xlab = expression(log(lambda[alpha])), ylab = expression(log(lambda[gamma])), main = paste(ic_label, "Heatmap"))
  contour(log_lam_a, log_lam_g, ic_mat, add = TRUE, col = "black")
  points(best_lam_a, best_lam_g, pch = 4, col = "blue", lwd = 3, cex = 2)

  image(log_lam_a, log_lam_g, edf_mat, col = hcl.colors(50, "Viridis"), xlab = expression(log(lambda[alpha])), ylab = expression(log(lambda[gamma])), main = "EDF Heatmap")
  contour(log_lam_a, log_lam_g, edf_mat, add = TRUE, col = "white")
  points(best_lam_a, best_lam_g, pch = 4, col = "red", lwd = 3, cex = 2)

  idx_g <- which(x$path_data$lambda_gamma == x$best_model$lambda_gamma)
  matplot(log_lam_a, x$path_data$alphas[idx_g, ], type = "l", lty = 1, xlab = expression(log(lambda[alpha])), ylab = expression(alpha), main = "Alpha Path")
  abline(v = best_lam_a, col = "firebrick", lty = 2, lwd = 2)

  idx_a <- which(x$path_data$lambda_alpha == x$best_model$lambda_alpha)
  matplot(log_lam_g, x$path_data$gammas[idx_a, ], type = "l", lty = 1, xlab = expression(log(lambda[gamma])), ylab = expression(gamma), main = "Gamma Path")
  abline(v = best_lam_g, col = "firebrick", lty = 2, lwd = 2)
}
