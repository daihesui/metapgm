#' Diagnostic Plots
#'
#' Generic function for creating diagnostic plots for fitted models.
#'
#' @param x An object for which diagnostic plots should be produced.
#' @param ... Additional arguments passed to methods.
#' @export
diagplot <- function(x, ...) UseMethod("diagplot")

#' @noRd
pad_zero <- function(x, k0) {
  res <- numeric(length(x) + 1)
  res[-k0] <- x
  return(res)
}

#' @noRd
format_newX <- function(newX, p) {
  if (p == 0) {
    return(matrix(NA, nrow = 1, ncol = 0))
  }
  if (is.null(newX)) {
    return(matrix(0, nrow = 1, ncol = p))
  }
  if (!is.matrix(newX)) newX <- rbind(as.numeric(newX))
  if (ncol(newX) != p) stop(paste("newX must have", p, "columns."))
  return(newX)
}

#' Fit a Penalized Gaussian Mixture Model
#'
#' Fits a penalized Gaussian mixture model for meta-analysis, estimating the underlying distribution of true effect sizes.
#'
#' @param y Numeric vector of observed effect sizes.
#' @param v Numeric vector of sampling variances.
#' @param X Optional numeric matrix of covariates for meta-regression.
#' @param group Optional vector indicating cluster/study membership for robust variance estimation.
#' @param penalty_order Integer specifying the order of the difference penalty matrix. Defaults to 3.
#' @param K Integer specifying the number of mixture components. Defaults to 40.
#' @param lambda_alpha_grid Numeric vector of penalty values to evaluate.
#' @param choose_lambda Character string specifying the information criterion ("aic" or "bic").
#' @param alpha_init Optional numeric vector of initial values for the alpha parameters.
#' @param beta_init Optional numeric vector of initial values for the beta parameters.
#' @param mu_grid Optional numeric vector of component means.
#' @param get_ci Logical indicating whether to compute variance-covariance matrices for confidence intervals.
#' @param ... Additional arguments passed to \code{\link[trust]{trust}}.
#'
#' @return An object of class \code{res.pgm}.
#'
#' @importFrom trust trust
#' @importFrom stats dnorm pnorm qnorm median mad uniroot
#' @importFrom graphics abline axis contour image legend lines points polygon rect par plot hist
#' @importFrom grDevices adjustcolor hcl.colors
#' @export
rma.pgm <- function(y, v, X = NULL, group = NULL, penalty_order = 3, K = 30,
                    lambda_alpha_grid = exp(seq(-5, 5, by = 1)), choose_lambda = "aic",
                    alpha_init = NULL, beta_init = NULL, mu_grid = NULL, get_ci = TRUE, ...) {
  if (!(choose_lambda %in% c("aic", "bic"))) stop("choose_lambda must be 'aic' or 'bic'.")
  if (!is.numeric(y) || !is.numeric(v)) stop("'y' and 'v' must be numeric.")

  N <- length(y)
  if (length(v) != N) stop("Length of 'y' and 'v' must match.")
  if (any(is.na(y)) || any(is.na(v))) stop("Missing values (NA) are not allowed in 'y' or 'v'.")
  if (any(v <= 0)) stop("All sampling variances 'v' must be strictly positive.")

  if (!is.null(X)) {
    if (is.data.frame(X)) X <- as.matrix(X)
    if (!is.matrix(X) || !is.numeric(X)) stop("'X' must be a numeric matrix or data frame.")
    if (nrow(X) != N) stop("Number of rows in 'X' must match length of 'y'.")
    if (any(is.na(X))) stop("Missing values (NA) are not allowed in 'X'.")
    p <- ncol(X)
    if (is.null(colnames(X))) colnames(X) <- paste0("X", 1:p)
  } else {
    X <- matrix(0, nrow = N, ncol = 0)
    p <- 0
  }

  use_robust <- !is.null(group)
  if (is.null(group)) {
    group <- seq_len(N)
  } else {
    if (length(group) != N) stop("Length of 'group' must match 'y'.")
    if (any(is.na(group))) stop("Missing values (NA) are not allowed in 'group'.")
  }

  mu <- if (is.null(mu_grid)) generate_mu_grid(y, v, X, K) else mu_grid
  K <- length(mu)
  tau_c <- (2 / 3) * (mu[2] - mu[1])
  k0 <- floor((K - 1) / 2) + 1
  D <- diff(diag(K), differences = penalty_order)
  DtD <- crossprod(D)

  if (is.null(alpha_init)) alpha_init <- rep(0, K - 1)

  if (p > 0 && is.null(beta_init)) {
    fit_init <- stats::lm.wfit(x = cbind(1, X), y = y, w = 1 / v)
    beta_init <- unname(fit_init$coefficients[-1])
    beta_init[is.na(beta_init)] <- 0
  }

  path_results <- fit_pgm_path(y, v, X, lambda_alpha_grid, choose_lambda, alpha_init, beta_init, mu, tau_c, DtD, k0, ...)
  ics <- sapply(path_results, function(x) x$ic)
  best_idx <- which.min(ics)
  best <- path_results[[best_idx]]

  if (is.null(best) || is.infinite(best$ic)) stop("Optimization failed for all lambda values.\n")

  if (get_ci) {
    if (use_robust) {
      V_eta <- calc_V_robust(y, v, X, group, mu, tau_c, best$alpha, best$beta, best$Jp_inv, k0)
    } else {
      V_eta <- best$Jp_inv %*% best$J %*% best$Jp_inv
    }
  } else {
    V_eta <- NULL
  }

  alpha_star <- pad_zero(best$alpha, k0)
  max_alpha <- max(alpha_star)
  w <- exp(alpha_star - max_alpha) / sum(exp(alpha_star - max_alpha))

  out <- list(
    best_model = list(alpha = best$alpha, beta = best$beta, w = w, V_eta = V_eta, ic = best$ic, edf = best$edf),
    path_data = list(lambda = lambda_alpha_grid, ic = ics, edf = sapply(path_results, function(x) x$edf)),
    data = list(y = y, X = X, group = group, choose_lambda = choose_lambda, mu = mu, tau_c = tau_c, k0 = k0, get_ci = get_ci),
    p = p
  )
  class(out) <- "res.pgm"

  return(out)
}

#' Summarize a PGM Model
#'
#' @param object An object of class \code{res.pgm}.
#' @param level Numeric scalar between 0 and 1 indicating the confidence level. Defaults to 0.95.
#' @param less_than Optional numeric vector to calculate the cumulative probability P(theta < val).
#' @param ... Additional arguments passed to other methods.
#' @return An object of class \code{summary.res.pgm}.
#' @export
summary.res.pgm <- function(object, level = 0.95, less_than = NULL, ...) {
  get_ci <- isTRUE(object$data$get_ci) || is.null(object$data$get_ci)
  mu <- object$data$mu
  tau_c <- object$data$tau_c
  k0 <- object$data$k0
  K <- length(mu)
  p <- object$p
  w <- object$best_model$w

  if (get_ci) {
    V_alpha <- object$best_model$V_eta[1:(K - 1), 1:(K - 1), drop = FALSE]
  }
  z_crit <- qnorm(1 - (1 - level) / 2)

  mu_est <- sum(w * mu)
  if (get_ci) {
    grad_mu <- as.numeric(w * (mu - mu_est))[-k0]
    var_mu <- as.numeric(crossprod(grad_mu, V_alpha %*% grad_mu))
    se_mu <- sqrt(max(0, var_mu))
    zval_mu <- ifelse(se_mu > 0, mu_est / se_mu, NA)
    pval_mu <- ifelse(!is.na(zval_mu), 2 * (1 - pnorm(abs(zval_mu))), NA)
    ci_mu <- c(mu_est - z_crit * se_mu, mu_est + z_crit * se_mu)
  } else {
    se_mu <- NA
    zval_mu <- NA
    pval_mu <- NA
    ci_mu <- c(NA, NA)
  }

  res_mat <- matrix(c(mu_est, se_mu, zval_mu, pval_mu, ci_mu[1], ci_mu[2]), nrow = 1)
  rownames(res_mat) <- "Mean (\u03bc)"

  mu2_est <- sum(w * mu^2)
  sigma2_est <- tau_c^2 + mu2_est - mu_est^2

  if (get_ci) {
    grad_sigma2 <- as.numeric(w * ((mu^2 - mu2_est) - 2 * mu_est * (mu - mu_est)))[-k0]
    var_sigma2 <- as.numeric(crossprod(grad_sigma2, V_alpha %*% grad_sigma2))
    se_sigma2 <- sqrt(max(0, var_sigma2))
    zval_sigma2 <- ifelse(se_sigma2 > 0, sigma2_est / se_sigma2, NA)
    pval_sigma2 <- ifelse(!is.na(zval_sigma2), 2 * (1 - pnorm(abs(zval_sigma2))), NA)
    ci_lb_sigma2 <- max(0, sigma2_est - z_crit * se_sigma2)
    ci_ub_sigma2 <- sigma2_est + z_crit * se_sigma2
  } else {
    se_sigma2 <- NA
    zval_sigma2 <- NA
    pval_sigma2 <- NA
    ci_lb_sigma2 <- NA
    ci_ub_sigma2 <- NA
  }

  res_mat <- rbind(res_mat, c(sigma2_est, se_sigma2, zval_sigma2, pval_sigma2, ci_lb_sigma2, ci_ub_sigma2))
  rownames(res_mat)[2] <- "Variance (\u03c3\u00b2)"

  if (!is.null(less_than)) {
    prob_list <- lapply(less_than, function(val) {
      p_c <- pnorm(val, mean = mu, sd = tau_c)
      prob_val <- sum(w * p_c)
      if (get_ci) {
        grad_p <- ((p_c - prob_val) * w)[-k0]
        var_p <- as.numeric(crossprod(grad_p, V_alpha %*% grad_p))
        se_p <- sqrt(max(0, var_p))
        zval_p <- ifelse(se_p > 0, prob_val / se_p, NA)
        pval_p <- ifelse(!is.na(zval_p), 2 * (1 - pnorm(abs(zval_p))), NA)
        ci_p <- c(max(0, prob_val - z_crit * se_p), min(1, prob_val + z_crit * se_p))
      } else {
        se_p <- NA
        zval_p <- NA
        pval_p <- NA
        ci_p <- c(NA, NA)
      }
      c(prob_val, se_p, zval_p, pval_p, ci_p[1], ci_p[2])
    })

    prob_stats <- do.call(rbind, prob_list)
    rownames(prob_stats) <- paste0("P(\u03b8 < ", less_than, ")")
    res_mat <- rbind(res_mat, prob_stats)
  }

  colnames(res_mat) <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")

  if (p > 0) {
    beta_est <- object$best_model$beta
    if (get_ci) {
      idx_beta <- (K):(K - 1 + p)
      V_beta <- object$best_model$V_eta[idx_beta, idx_beta, drop = FALSE]
      se_beta <- sqrt(pmax(0, diag(V_beta)))
      zval_beta <- ifelse(se_beta > 0, beta_est / se_beta, NA)
      pval_beta <- ifelse(!is.na(zval_beta), 2 * (1 - pnorm(abs(zval_beta))), NA)
      ci_lb_beta <- beta_est - z_crit * se_beta
      ci_ub_beta <- beta_est + z_crit * se_beta
    } else {
      se_beta <- rep(NA, p)
      zval_beta <- rep(NA, p)
      pval_beta <- rep(NA, p)
      ci_lb_beta <- rep(NA, p)
      ci_ub_beta <- rep(NA, p)
    }
    beta_mat <- cbind(beta_est, se_beta, zval_beta, pval_beta, ci_lb_beta, ci_ub_beta)
    rownames(beta_mat) <- colnames(object$data$X)
    colnames(beta_mat) <- c("estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
  } else {
    beta_mat <- NULL
  }

  out <- list(
    res_mat = res_mat,
    beta_mat = beta_mat,
    n_obs = length(object$data$y),
    n_groups = length(unique(object$data$group)),
    ic_name = toupper(object$data$choose_lambda),
    ic_val = object$best_model$ic,
    edf = object$best_model$edf,
    p = p
  )
  class(out) <- "summary.res.pgm"

  return(out)
}

#' Print PGM Summary
#'
#' @param x An object of class \code{summary.res.pgm}.
#' @param digits Integer indicating the number of decimal places to print. Defaults to 4.
#' @param ... Additional arguments passed to other methods.
#' @export
print.summary.res.pgm <- function(x, digits = 4, ...) {
  cat("\n======================================================\n")
  cat("          Penalized Gaussian Mixture Results          \n")
  cat("======================================================\n")
  cat(sprintf("Observations : %d \nGroups       : %d\n", x$n_obs, x$n_groups))
  cat(sprintf("Model Fit    : %s = %.*f | EDF = %.*f\n", x$ic_name, digits, x$ic_val, digits, x$edf))

  if (x$p == 0) {
    cat("\n--- Parameter Estimates ---\n")
  } else {
    x_str <- if (x$p == 1) "0" else sprintf("(%s)", paste(rep("0", x$p), collapse = ", "))
    cat(sprintf("\n--- Parameter Estimates (X = %s) ---\n", x_str))
  }

  print(round(x$res_mat, digits = digits))

  if (!is.null(x$beta_mat)) {
    cat("\n--- Covariate Coefficients ---\n")
    print(round(x$beta_mat, digits = digits))
  }

  cat("\n")
  invisible(x)
}

#' Predict Method for PGM Model
#'
#' @param object An object of class \code{res.pgm}.
#' @param level Numeric scalar indicating the prediction interval level. Defaults to 0.95.
#' @param newX Optional numeric matrix of new moderator values for prediction.
#' @param digits Integer indicating the number of decimal places to print. Defaults to 4.
#' @param ... Additional arguments passed to other methods.
#' @return A matrix containing the lower and upper prediction intervals.
#' @export
predict.res.pgm <- function(object, level = 0.95, newX = NULL, digits = 4, ...) {
  mu <- object$data$mu
  tau_c <- object$data$tau_c
  w <- object$best_model$w
  p <- object$p
  K <- length(mu)

  newX <- format_newX(newX, p)
  shifts <- if (p > 0) as.numeric(newX %*% object$best_model$beta) else c(0)

  res_list <- vector("list", length(shifts))

  for (k in seq_along(shifts)) {
    shift <- shifts[k]
    cdf_fun <- function(q) sum(w * pnorm(q, mean = mu + shift, sd = tau_c))

    mu_shifted <- mu + shift
    search_range <- c(mu_shifted[1] - 5 * tau_c, mu_shifted[K] + 5 * tau_c)

    pi_l <- tryCatch(uniroot(function(q) cdf_fun(q) - ((1 - level) / 2), interval = search_range, extendInt = "yes")$root, error = function(e) NA)
    pi_u <- tryCatch(uniroot(function(q) cdf_fun(q) - (1 - (1 - level) / 2), interval = search_range, extendInt = "yes")$root, error = function(e) NA)
    res_list[[k]] <- c(pi.lb = pi_l, pi.ub = pi_u)
  }

  res_mat <- do.call(rbind, res_list)

  if (p == 0) {
    rownames(res_mat) <- "Mean"
  } else {
    rownames(res_mat) <- apply(newX, 1, function(row) {
      if (p == 1) paste0("X = ", round(row[1], digits)) else paste0("X = (", paste(round(row, digits), collapse = ", "), ")")
    })
  }

  cat(sprintf("\n--- %g%% Prediction Intervals (PI) ---\n", level * 100))
  print(round(res_mat, digits = digits))
  cat("\n")

  invisible(res_mat)
}

#' Plot Method for PGM Model
#'
#' @param x An object of class \code{res.pgm}.
#' @param level Confidence level for the uncertainty bands. Defaults to 0.95.
#' @param newX Optional numeric matrix of new moderator values.
#' @param add_hist Logical indicating whether to overlay a histogram of the observed effect sizes. Defaults to FALSE.
#' @param show_ci Logical indicating whether to display confidence bands around the density. Defaults to FALSE.
#' @param ... Additional graphic arguments.
#' @export
plot.res.pgm <- function(x, level = 0.95, newX = NULL, add_hist = FALSE, show_ci = FALSE, ...) {
  get_ci <- isTRUE(x$data$get_ci) || is.null(x$data$get_ci)
  if (!get_ci) show_ci <- FALSE

  mu <- x$data$mu
  tau_c <- x$data$tau_c
  k0 <- x$data$k0
  K <- length(mu)
  p <- x$p

  newX <- format_newX(newX, p)
  shifts <- if (p > 0) as.numeric(newX %*% x$best_model$beta) else c(0)

  w <- x$best_model$w

  if (get_ci) {
    V_eta <- x$best_model$V_eta
    V_alpha <- V_eta[1:(K - 1), 1:(K - 1), drop = FALSE]
  }

  z_crit <- qnorm(1 - (1 - level) / 2)
  cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  colors <- rep(cb_palette, length.out = length(shifts))

  plot_data <- list()
  global_x_min <- Inf
  global_x_max <- -Inf
  global_y_max <- -Inf

  for (k in seq_along(shifts)) {
    shift <- shifts[k]
    mu_shifted <- mu + shift

    x_grid <- seq(mu_shifted[1] - 5 * tau_c, mu_shifted[K] + 5 * tau_c, length.out = 500)
    n_grid <- length(x_grid)
    global_x_min <- min(global_x_min, x_grid[1])
    global_x_max <- max(global_x_max, x_grid[n_grid])

    Phi_grid <- matrix(dnorm(x_grid, mean = rep(mu_shifted, each = n_grid), sd = tau_c), nrow = n_grid, ncol = K)
    pdf_vals <- as.numeric(Phi_grid %*% w)

    y_pgm_lower <- NULL
    y_pgm_upper <- NULL

    if (show_ci && get_ci) {
      grad_f_alpha <- ((Phi_grid - pdf_vals) * rep(w, each = n_grid))[, -k0, drop = FALSE]
      if (p > 0) {
        R_mat_grid <- (x_grid - rep(mu_shifted, each = n_grid)) / (tau_c^2)
        grad_f_beta <- outer(rowSums((Phi_grid * rep(w, each = n_grid)) * R_mat_grid), newX[k, ])
        grad_f_full <- cbind(grad_f_alpha, grad_f_beta)
        var_f <- rowSums((grad_f_full %*% V_eta) * grad_f_full)
      } else {
        var_f <- rowSums((grad_f_alpha %*% V_alpha) * grad_f_alpha)
      }
      se_pdf <- sqrt(pmax(0, var_f))
      y_pgm_lower <- pmax(0, pdf_vals - z_crit * se_pdf)
      y_pgm_upper <- pdf_vals + z_crit * se_pdf
      global_y_max <- max(global_y_max, max(y_pgm_upper, na.rm = TRUE))
    } else {
      global_y_max <- max(global_y_max, max(pdf_vals, na.rm = TRUE))
    }

    plot_data[[k]] <- list(x = x_grid, pdf = pdf_vals, lower = y_pgm_lower, upper = y_pgm_upper)
  }

  y <- x$data$y
  par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.axis = 1, cex.lab = 1.1)

  if (add_hist) {
    h <- hist(y, plot = FALSE, breaks = 30)
    y_max_final <- max(c(h$density, global_y_max), na.rm = TRUE) * 1.25
  } else {
    y_max_final <- global_y_max * 1.25
  }

  plot(NULL, xlim = c(global_x_min, global_x_max), ylim = c(0, y_max_final), main = "PGM Estimated Density", xlab = "Effect Size", ylab = "Density", axes = FALSE, frame.plot = FALSE)
  axis(1, tick = TRUE)
  axis(2, las = 1, tick = TRUE)

  if (add_hist) rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$density, col = "gray92", border = "gray75", lwd = 1)

  for (k in seq_along(shifts)) {
    if (show_ci && get_ci) {
      polygon(c(plot_data[[k]]$x, rev(plot_data[[k]]$x)), c(plot_data[[k]]$lower, rev(plot_data[[k]]$upper)), col = adjustcolor(colors[k], alpha.f = 0.2), border = NA)
    }
    lines(plot_data[[k]]$x, plot_data[[k]]$pdf, lwd = 2, col = colors[k])
  }

  if (p > 0) {
    legend_labels <- apply(newX, 1, function(row) {
      if (p == 1) paste0("X = ", round(row[1], 4)) else paste0("X = (", paste(round(row, 4), collapse = ", "), ")")
    })
    legend("topright", legend = legend_labels, col = colors[1:length(shifts)], lwd = 2, bty = "n")
  }
}

#' Diagnostic Plots for PGM Model
#'
#' @param x An object of class \code{res.pgm}.
#' @param ... Additional graphic arguments.
#' @export
diagplot.res.pgm <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(2, 1), mar = c(4, 5, 2, 2), cex.axis = 1, cex.lab = 1.1)

  log_lambda_alpha <- log(x$path_data$lambda)
  best_idx <- which.min(x$path_data$ic)
  best_log_lambda_alpha <- log_lambda_alpha[best_idx]
  ic_label <- toupper(x$data$choose_lambda)

  plot_ic <- x$path_data$ic
  plot_ic[is.infinite(plot_ic)] <- NA

  plot(log_lambda_alpha, plot_ic, type = "b", pch = 16, col = "#0072B2", xlab = expression(log(lambda[alpha])), ylab = ic_label, main = paste(ic_label, "Profile"))
  abline(v = best_log_lambda_alpha, col = "#D55E00", lty = 2, lwd = 2)

  plot(log_lambda_alpha, x$path_data$edf, type = "b", pch = 16, col = "#009E73", xlab = expression(log(lambda[alpha])), ylab = "EDF", main = "EDF Profile")
  abline(v = best_log_lambda_alpha, col = "#D55E00", lty = 2, lwd = 2)
}

#' Fit a Conditional Penalized Gaussian Mixture Model
#'
#' Fits a conditional PGM model where mixture weights vary as a function of a moderator.
#'
#' @param y Numeric vector of observed effect sizes.
#' @param v Numeric vector of sampling variances.
#' @param z Numeric vector of the conditional continuous moderator variable.
#' @param group Optional vector indicating cluster/study membership for robust variance estimation.
#' @param penalty_order Integer specifying the order of the difference penalty matrix. Defaults to 3.
#' @param K Integer specifying the number of mixture components. Defaults to 40.
#' @param lambda_alpha_grid Numeric vector of penalty values for the intercept weights (alpha).
#' @param lambda_gamma_grid Numeric vector of penalty values for the slope weights (gamma).
#' @param choose_lambda Character string specifying the information criterion ("aic" or "bic").
#' @param alpha_init Optional numeric vector of initial values for alpha parameters.
#' @param gamma_init Optional numeric vector of initial values for gamma parameters.
#' @param mu_grid Optional numeric vector of component means.
#' @param get_ci Logical indicating whether to compute variance-covariance matrices for confidence intervals.
#' @param ... Additional arguments passed to \code{\link[trust]{trust}}.
#'
#' @return An object of class \code{res.pgm.cond}.
#' @export
rma.pgm.cond <- function(y, v, z, group = NULL, penalty_order = 3, K = 30,
                         lambda_alpha_grid = exp(seq(-5, 5, by = 1)),
                         lambda_gamma_grid = exp(seq(-5, 5, by = 1)),
                         choose_lambda = "aic", alpha_init = NULL, gamma_init = NULL, mu_grid = NULL, get_ci = TRUE, ...) {
  if (!(choose_lambda %in% c("aic", "bic"))) stop("choose_lambda must be 'aic' or 'bic'.")
  if (!is.numeric(y) || !is.numeric(v)) stop("'y' and 'v' must be numeric.")

  N <- length(y)
  if (length(v) != N) stop("Length of 'y' and 'v' must match.")
  if (any(is.na(y)) || any(is.na(v))) stop("Missing values (NA) are not allowed in 'y' or 'v'.")
  if (any(v <= 0)) stop("All sampling variances 'v' must be strictly positive.")
  if (is.factor(z)) stop("'z' cannot be a factor. Please provide dummy variables for categorical moderators.")

  z <- as.numeric(z)
  if (length(z) != N) stop("Length of 'z' must match 'y'.")
  if (any(is.na(z))) stop("Missing values (NA) are not allowed in 'z'.")

  use_robust <- !is.null(group)
  if (is.null(group)) {
    group <- seq_len(N)
  } else {
    if (length(group) != N) stop("Length of 'group' must match 'y'.")
    if (any(is.na(group))) stop("Missing values (NA) are not allowed in 'group'.")
  }

  mu <- if (is.null(mu_grid)) generate_mu_grid(y, v, matrix(0, nrow = N, ncol = 0), K) else mu_grid
  K <- length(mu)
  tau_c <- (2 / 3) * (mu[2] - mu[1])
  k0 <- floor((K - 1) / 2) + 1
  D <- diff(diag(K), differences = penalty_order)
  DtD <- crossprod(D)

  if (is.null(alpha_init)) alpha_init <- rep(0, K - 1)
  if (is.null(gamma_init)) gamma_init <- rep(0, K - 1)
  eta_init <- c(alpha_init, gamma_init)

  path_results <- fit_pgm_cond_path(y, v, z, lambda_alpha_grid, lambda_gamma_grid, choose_lambda, eta_init, mu, tau_c, DtD, k0, ...)
  ics <- sapply(path_results, function(x) x$ic)
  best_idx <- which.min(ics)
  best <- path_results[[best_idx]]

  if (is.null(best) || is.infinite(best$ic)) stop("Optimization failed for all lambda combinations.\n")

  if (get_ci) {
    if (use_robust) {
      V_eta <- calc_V_cond_robust(y, v, z, group, mu, tau_c, best$alpha, best$gamma, best$Jp_inv, k0)
    } else {
      V_eta <- best$Jp_inv %*% best$J %*% best$Jp_inv
    }
  } else {
    V_eta <- NULL
  }

  out <- list(
    best_model = list(alpha = best$alpha, gamma = best$gamma, V_eta = V_eta, ic = best$ic, edf = best$edf, lambda_alpha = best$lambda_alpha, lambda_gamma = best$lambda_gamma),
    path_data = list(
      lambda_alpha = sapply(path_results, function(x) x$lambda_alpha),
      lambda_gamma = sapply(path_results, function(x) x$lambda_gamma),
      ic = ics,
      edf = sapply(path_results, function(x) x$edf)
    ),
    data = list(y = y, group = group, choose_lambda = choose_lambda, mu = mu, tau_c = tau_c, k0 = k0, get_ci = get_ci)
  )
  class(out) <- "res.pgm.cond"

  return(out)
}

#' Summarize a Conditional PGM Model
#'
#' @param object An object of class \code{res.pgm.cond}.
#' @param level Numeric scalar between 0 and 1 indicating the confidence level. Defaults to 0.95.
#' @param less_than Optional numeric vector to calculate P(theta < val).
#' @param newz Optional numeric vector of conditional z values to evaluate the summary at.
#' @param compare Optional numeric vector of length 2 or a two-column matrix specifying pairs of z values to compare.
#' @param ... Additional arguments passed to other methods.
#' @return An object of class \code{summary.res.pgm.cond}.
#' @export
summary.res.pgm.cond <- function(object, level = 0.95, less_than = NULL, newz = NULL, compare = NULL, ...) {
  get_ci <- isTRUE(object$data$get_ci) || is.null(object$data$get_ci)
  mu <- object$data$mu
  tau_c <- object$data$tau_c
  k0 <- object$data$k0

  if (get_ci) {
    V_eta <- object$best_model$V_eta
  }
  z_crit <- qnorm(1 - (1 - level) / 2)

  alpha_star <- pad_zero(object$best_model$alpha, k0)
  gamma_star <- pad_zero(object$best_model$gamma, k0)

  newz_vals <- if (is.null(newz)) 0 else as.numeric(newz)
  all_res <- list()

  for (z_val in newz_vals) {
    eta_z <- alpha_star + gamma_star * z_val
    max_eta_z <- max(eta_z)
    w <- exp(eta_z - max_eta_z) / sum(exp(eta_z - max_eta_z))

    mu_est <- sum(w * mu)

    if (get_ci) {
      J_w <- diag(w) - tcrossprod(w)
      J_mu_alpha <- as.numeric(crossprod(J_w, mu))[-k0]
      J_mu_gamma <- z_val * J_mu_alpha
      grad_mu <- c(J_mu_alpha, J_mu_gamma)
      var_mu <- as.numeric(crossprod(grad_mu, V_eta %*% grad_mu))
      se_mu <- sqrt(max(0, var_mu))
      zval_mu <- ifelse(se_mu > 0, mu_est / se_mu, NA)
      pval_mu <- ifelse(!is.na(zval_mu), 2 * (1 - pnorm(abs(zval_mu))), NA)
      ci_mu <- c(mu_est - z_crit * se_mu, mu_est + z_crit * se_mu)
    } else {
      se_mu <- NA
      zval_mu <- NA
      pval_mu <- NA
      ci_mu <- c(NA, NA)
    }

    res_mat <- matrix(c(mu_est, se_mu, zval_mu, pval_mu, ci_mu[1], ci_mu[2]), nrow = 1)
    rownames(res_mat) <- "Mean (\u03bc)"

    mu2_est <- sum(w * mu^2)
    sigma2_est <- tau_c^2 + mu2_est - mu_est^2

    if (get_ci) {
      J_mu2_alpha <- as.numeric(crossprod(J_w, mu^2))[-k0]
      J_sigma2_alpha <- J_mu2_alpha - 2 * mu_est * J_mu_alpha
      J_sigma2_gamma <- z_val * J_sigma2_alpha

      grad_sigma2 <- c(J_sigma2_alpha, J_sigma2_gamma)
      var_sigma2 <- as.numeric(crossprod(grad_sigma2, V_eta %*% grad_sigma2))
      se_sigma2 <- sqrt(max(0, var_sigma2))
      zval_sigma2 <- ifelse(se_sigma2 > 0, sigma2_est / se_sigma2, NA)
      pval_sigma2 <- ifelse(!is.na(zval_sigma2), 2 * (1 - pnorm(abs(zval_sigma2))), NA)
      ci_lb_sigma2 <- max(0, sigma2_est - z_crit * se_sigma2)
      ci_ub_sigma2 <- sigma2_est + z_crit * se_sigma2
    } else {
      se_sigma2 <- NA
      zval_sigma2 <- NA
      pval_sigma2 <- NA
      ci_lb_sigma2 <- NA
      ci_ub_sigma2 <- NA
    }

    res_mat <- rbind(res_mat, c(sigma2_est, se_sigma2, zval_sigma2, pval_sigma2, ci_lb_sigma2, ci_ub_sigma2))
    rownames(res_mat)[2] <- "Variance (\u03c3\u00b2)"

    if (!is.null(less_than)) {
      prob_list <- lapply(less_than, function(val) {
        p_c <- pnorm(val, mean = mu, sd = tau_c)
        prob_val <- sum(w * p_c)
        if (get_ci) {
          J_p_alpha <- as.numeric(crossprod(J_w, p_c))[-k0]
          J_p_gamma <- z_val * J_p_alpha
          grad_p <- c(J_p_alpha, J_p_gamma)
          var_p <- as.numeric(crossprod(grad_p, V_eta %*% grad_p))
          se_p <- sqrt(max(0, var_p))
          zval_p <- ifelse(se_p > 0, prob_val / se_p, NA)
          pval_p <- ifelse(!is.na(zval_p), 2 * (1 - pnorm(abs(zval_p))), NA)
          ci_p <- c(max(0, prob_val - z_crit * se_p), min(1, prob_val + z_crit * se_p))
        } else {
          se_p <- NA
          zval_p <- NA
          pval_p <- NA
          ci_p <- c(NA, NA)
        }
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
    for (i in seq_len(nrow(compare))) {
      z1 <- compare[i, 1]
      z2 <- compare[i, 2]

      eta_z1 <- alpha_star + gamma_star * z1
      w_z1 <- exp(eta_z1 - max(eta_z1)) / sum(exp(eta_z1 - max(eta_z1)))

      eta_z2 <- alpha_star + gamma_star * z2
      w_z2 <- exp(eta_z2 - max(eta_z2)) / sum(exp(eta_z2 - max(eta_z2)))

      mu_diff <- sum(w_z2 * mu) - sum(w_z1 * mu)

      if (get_ci) {
        J_w1 <- diag(w_z1) - tcrossprod(w_z1)
        J_w2 <- diag(w_z2) - tcrossprod(w_z2)

        J_mu_alpha1 <- as.numeric(crossprod(J_w1, mu))[-k0]
        grad_mu1 <- c(J_mu_alpha1, z1 * J_mu_alpha1)

        J_mu_alpha2 <- as.numeric(crossprod(J_w2, mu))[-k0]
        grad_mu2 <- c(J_mu_alpha2, z2 * J_mu_alpha2)

        grad_mu_diff <- grad_mu2 - grad_mu1
        var_mu_diff <- as.numeric(crossprod(grad_mu_diff, V_eta %*% grad_mu_diff))

        se_mu_diff <- sqrt(max(0, var_mu_diff))
        zval_mu_diff <- ifelse(se_mu_diff > 0, mu_diff / se_mu_diff, NA)
        pval_mu_diff <- ifelse(!is.na(zval_mu_diff), 2 * (1 - pnorm(abs(zval_mu_diff))), NA)
        ci_mu_diff <- c(mu_diff - z_crit * se_mu_diff, mu_diff + z_crit * se_mu_diff)
      } else {
        se_mu_diff <- NA
        zval_mu_diff <- NA
        pval_mu_diff <- NA
        ci_mu_diff <- c(NA, NA)
      }

      res_mat_diff <- matrix(c(mu_diff, se_mu_diff, zval_mu_diff, pval_mu_diff, ci_mu_diff[1], ci_mu_diff[2]), nrow = 1)
      rownames(res_mat_diff) <- "Mean Diff (\u03bc)"

      mu2_z1 <- sum(w_z1 * mu^2)
      sigma2_z1 <- tau_c^2 + mu2_z1 - sum(w_z1 * mu)^2
      mu2_z2 <- sum(w_z2 * mu^2)
      sigma2_z2 <- tau_c^2 + mu2_z2 - sum(w_z2 * mu)^2
      sigma2_diff <- sigma2_z2 - sigma2_z1

      if (get_ci) {
        J_mu2_alpha1 <- as.numeric(crossprod(J_w1, mu^2))[-k0]
        J_sigma2_alpha1 <- J_mu2_alpha1 - 2 * sum(w_z1 * mu) * J_mu_alpha1
        grad_sigma2_1 <- c(J_sigma2_alpha1, z1 * J_sigma2_alpha1)

        J_mu2_alpha2 <- as.numeric(crossprod(J_w2, mu^2))[-k0]
        J_sigma2_alpha2 <- J_mu2_alpha2 - 2 * sum(w_z2 * mu) * J_mu_alpha2
        grad_sigma2_2 <- c(J_sigma2_alpha2, z2 * J_sigma2_alpha2)

        grad_sigma2_diff <- grad_sigma2_2 - grad_sigma2_1
        var_sigma2_diff <- as.numeric(crossprod(grad_sigma2_diff, V_eta %*% grad_sigma2_diff))

        se_sigma2_diff <- sqrt(max(0, var_sigma2_diff))
        zval_sigma2_diff <- ifelse(se_sigma2_diff > 0, sigma2_diff / se_sigma2_diff, NA)
        pval_sigma2_diff <- ifelse(!is.na(zval_sigma2_diff), 2 * (1 - pnorm(abs(zval_sigma2_diff))), NA)
        ci_lb_sigma2_diff <- sigma2_diff - z_crit * se_sigma2_diff
        ci_ub_sigma2_diff <- sigma2_diff + z_crit * se_sigma2_diff
      } else {
        se_sigma2_diff <- NA
        zval_sigma2_diff <- NA
        pval_sigma2_diff <- NA
        ci_lb_sigma2_diff <- NA
        ci_ub_sigma2_diff <- NA
      }

      res_mat_diff <- rbind(res_mat_diff, c(sigma2_diff, se_sigma2_diff, zval_sigma2_diff, pval_sigma2_diff, ci_lb_sigma2_diff, ci_ub_sigma2_diff))
      rownames(res_mat_diff)[2] <- "\u0394 Variance (\u03c3\u00b2)"

      if (!is.null(less_than)) {
        prob_list_diff <- lapply(less_than, function(val) {
          p_c <- pnorm(val, mean = mu, sd = tau_c)
          prob_val1 <- sum(w_z1 * p_c)
          prob_val2 <- sum(w_z2 * p_c)
          prob_diff <- prob_val2 - prob_val1

          if (get_ci) {
            J_p_alpha1 <- as.numeric(crossprod(J_w1, p_c))[-k0]
            grad_p1 <- c(J_p_alpha1, z1 * J_p_alpha1)

            J_p_alpha2 <- as.numeric(crossprod(J_w2, p_c))[-k0]
            grad_p2 <- c(J_p_alpha2, z2 * J_p_alpha2)

            grad_p_diff <- grad_p2 - grad_p1
            var_p_diff <- as.numeric(crossprod(grad_p_diff, V_eta %*% grad_p_diff))

            se_p_diff <- sqrt(max(0, var_p_diff))
            zval_p_diff <- ifelse(se_p_diff > 0, prob_diff / se_p_diff, NA)
            pval_p_diff <- ifelse(!is.na(zval_p_diff), 2 * (1 - pnorm(abs(zval_p_diff))), NA)
            ci_p_diff <- c(prob_diff - z_crit * se_p_diff, prob_diff + z_crit * se_p_diff)
          } else {
            se_p_diff <- NA
            zval_p_diff <- NA
            pval_p_diff <- NA
            ci_p_diff <- c(NA, NA)
          }

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

  out <- list(all_res = all_res, compare_res = compare_res, n_obs = length(object$data$y), n_groups = length(unique(object$data$group)), ic_name = toupper(object$data$choose_lambda), ic_val = object$best_model$ic, edf = object$best_model$edf)
  class(out) <- "summary.res.pgm.cond"

  return(out)
}

#' Print Conditional PGM Summary
#'
#' @param x An object of class \code{summary.res.pgm.cond}.
#' @param digits Integer indicating the number of decimal places to print. Defaults to 4.
#' @param ... Additional arguments passed to other methods.
#' @export
print.summary.res.pgm.cond <- function(x, digits = 4, ...) {
  cat("\n======================================================\n")
  cat("    Conditional Penalized Gaussian Mixture Results    \n")
  cat("======================================================\n")
  cat(sprintf("Observations : %d \nGroups       : %d\n", x$n_obs, x$n_groups))
  cat(sprintf("Model Fit    : %s = %.*f | EDF = %.*f\n", x$ic_name, digits, x$ic_val, digits, x$edf))

  for (header in names(x$all_res)) {
    cat(sprintf("\n--- Parameter Estimates (%s) ---\n", header))
    print(round(x$all_res[[header]], digits = digits))
  }

  if (!is.null(x$compare_res)) {
    for (header in names(x$compare_res)) {
      cat(sprintf("\n--- Comparison (%s) ---\n", header))
      print(round(x$compare_res[[header]], digits = digits))
    }
  }

  cat("\n")
  invisible(x)
}

#' Predict Method for Conditional PGM Model
#'
#' @param object An object of class \code{res.pgm.cond}.
#' @param level Numeric scalar indicating prediction interval level. Defaults to 0.95.
#' @param newz Optional numeric vector of conditional z values to predict at.
#' @param digits Integer indicating the number of decimal places to print. Defaults to 4.
#' @param ... Additional arguments passed to other methods.
#' @export
predict.res.pgm.cond <- function(object, level = 0.95, newz = NULL, digits = 4, ...) {
  mu <- object$data$mu
  tau_c <- object$data$tau_c
  k0 <- object$data$k0
  K <- length(mu)

  alpha_star <- pad_zero(object$best_model$alpha, k0)
  gamma_star <- pad_zero(object$best_model$gamma, k0)

  newz_vals <- if (is.null(newz)) 0 else as.numeric(newz)
  res_list <- vector("list", length(newz_vals))

  search_range <- c(mu[1] - 5 * tau_c, mu[K] + 5 * tau_c)

  for (k in seq_along(newz_vals)) {
    z_val <- newz_vals[k]
    eta_z <- alpha_star + gamma_star * z_val
    max_eta_z <- max(eta_z)
    w <- exp(eta_z - max_eta_z) / sum(exp(eta_z - max_eta_z))

    cdf_fun <- function(q) sum(w * pnorm(q, mean = mu, sd = tau_c))
    pi_l <- tryCatch(uniroot(function(q) cdf_fun(q) - ((1 - level) / 2), interval = search_range, extendInt = "yes")$root, error = function(e) NA)
    pi_u <- tryCatch(uniroot(function(q) cdf_fun(q) - (1 - (1 - level) / 2), interval = search_range, extendInt = "yes")$root, error = function(e) NA)
    res_list[[k]] <- c(pi.lb = pi_l, pi.ub = pi_u)
  }

  res_mat <- do.call(rbind, res_list)
  rownames(res_mat) <- paste("z =", round(newz_vals, digits))

  cat(sprintf("\n--- %g%% Prediction Intervals (PI) ---\n", level * 100))
  print(round(res_mat, digits = digits))
  cat("\n")

  invisible(res_mat)
}

#' Plot Method for Conditional PGM Model
#'
#' @param x An object of class \code{res.pgm.cond}.
#' @param level Confidence level for uncertainty bands. Defaults to 0.95.
#' @param newz Optional numeric vector of z values to plot curves for.
#' @param add_hist Logical indicating whether to overlay a histogram. Defaults to FALSE.
#' @param show_ci Logical indicating whether to show confidence intervals. Defaults to FALSE.
#' @param ... Additional graphic arguments.
#' @export
plot.res.pgm.cond <- function(x, level = 0.95, newz = NULL, add_hist = FALSE, show_ci = FALSE, ...) {
  get_ci <- isTRUE(x$data$get_ci) || is.null(x$data$get_ci)
  if (!get_ci) show_ci <- FALSE

  mu <- x$data$mu
  tau_c <- x$data$tau_c
  k0 <- x$data$k0
  K <- length(mu)

  if (get_ci) {
    V_eta <- x$best_model$V_eta
  }

  z_crit <- qnorm(1 - (1 - level) / 2)

  alpha_star <- pad_zero(x$best_model$alpha, k0)
  gamma_star <- pad_zero(x$best_model$gamma, k0)

  newz_vals <- if (is.null(newz)) 0 else as.numeric(newz)
  cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  colors <- rep(cb_palette, length.out = length(newz_vals))

  plot_data <- list()
  global_y_max <- -Inf

  x_grid <- seq(mu[1] - 5 * tau_c, mu[K] + 5 * tau_c, length.out = 500)
  n_grid <- length(x_grid)
  Phi_grid <- matrix(dnorm(x_grid, mean = rep(mu, each = n_grid), sd = tau_c), nrow = n_grid, ncol = K)

  for (k in seq_along(newz_vals)) {
    z_val <- newz_vals[k]
    eta_z <- alpha_star + gamma_star * z_val
    max_eta_z <- max(eta_z)
    w <- exp(eta_z - max_eta_z) / sum(exp(eta_z - max_eta_z))
    pdf_vals <- as.numeric(Phi_grid %*% w)

    y_lower <- NULL
    y_upper <- NULL

    if (show_ci && get_ci) {
      J_w <- diag(w) - tcrossprod(w)
      grad_f_alpha <- (Phi_grid %*% J_w)[, -k0, drop = FALSE]
      grad_f_gamma <- z_val * grad_f_alpha
      grad_f_full <- cbind(grad_f_alpha, grad_f_gamma)
      var_f <- rowSums((grad_f_full %*% V_eta) * grad_f_full)
      se_pdf <- sqrt(pmax(0, var_f))
      y_lower <- pmax(0, pdf_vals - z_crit * se_pdf)
      y_upper <- pdf_vals + z_crit * se_pdf
      global_y_max <- max(global_y_max, max(y_upper, na.rm = TRUE))
    } else {
      global_y_max <- max(global_y_max, max(pdf_vals, na.rm = TRUE))
    }

    plot_data[[k]] <- list(pdf = pdf_vals, lower = y_lower, upper = y_upper)
  }

  y <- x$data$y
  par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.axis = 1, cex.lab = 1.1)

  if (add_hist) {
    h <- hist(y, plot = FALSE, breaks = 30)
    y_max_final <- max(c(h$density, global_y_max), na.rm = TRUE) * 1.25
  } else {
    y_max_final <- global_y_max * 1.25
  }

  plot(NULL, xlim = range(x_grid), ylim = c(0, y_max_final), main = "PGM Estimated Density", xlab = "Effect Size", ylab = "Density", axes = FALSE, frame.plot = FALSE)
  axis(1, tick = TRUE)
  axis(2, las = 1, tick = TRUE)

  if (add_hist) rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$density, col = "gray92", border = "gray75", lwd = 1)

  for (k in seq_along(newz_vals)) {
    if (show_ci && get_ci) {
      polygon(c(x_grid, rev(x_grid)), c(plot_data[[k]]$lower, rev(plot_data[[k]]$upper)), col = adjustcolor(colors[k], alpha.f = 0.2), border = NA)
    }
    lines(x_grid, plot_data[[k]]$pdf, lwd = 2, col = colors[k])
  }

  legend_labels <- paste("z =", newz_vals)
  legend("topright", legend = legend_labels, col = colors[1:length(newz_vals)], lwd = 2, bty = "n")
}

#' Diagnostic Plots for Conditional PGM Model
#'
#' @param x An object of class \code{res.pgm.cond}.
#' @param ... Additional graphic arguments.
#' @export
diagplot.res.pgm.cond <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  lambda_alpha_vals <- unique(x$path_data$lambda_alpha)
  lambda_gamma_vals <- unique(x$path_data$lambda_gamma)
  log_lambda_alpha <- log(lambda_alpha_vals)
  log_lambda_gamma <- log(lambda_gamma_vals)

  safe_ic <- x$path_data$ic
  safe_ic[is.infinite(safe_ic)] <- NA
  ic_mat <- matrix(safe_ic, nrow = length(lambda_alpha_vals), ncol = length(lambda_gamma_vals))
  edf_mat <- matrix(x$path_data$edf, nrow = length(lambda_alpha_vals), ncol = length(lambda_gamma_vals))

  best_lambda_alpha <- log(x$best_model$lambda_alpha)
  best_lambda_gamma <- log(x$best_model$lambda_gamma)
  ic_label <- toupper(x$data$choose_lambda)

  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1), cex.axis = 1, cex.lab = 1.1)

  image(log_lambda_alpha, log_lambda_gamma, ic_mat, col = hcl.colors(50, "Plasma", rev = TRUE), xlab = expression(log(lambda[alpha])), ylab = expression(log(lambda[gamma])), main = paste(ic_label, "Heatmap"))
  contour(log_lambda_alpha, log_lambda_gamma, ic_mat, add = TRUE, col = "black")
  points(best_lambda_alpha, best_lambda_gamma, pch = 4, col = "#009E73", lwd = 3, cex = 2)

  image(log_lambda_alpha, log_lambda_gamma, edf_mat, col = hcl.colors(50, "Viridis"), xlab = expression(log(lambda[alpha])), ylab = expression(log(lambda[gamma])), main = "EDF Heatmap")
  contour(log_lambda_alpha, log_lambda_gamma, edf_mat, add = TRUE, col = "white")
  points(best_lambda_alpha, best_lambda_gamma, pch = 4, col = "#D55E00", lwd = 3, cex = 2)
}

#' @noRd
generate_mu_grid <- function(y, v, X, K) {
  v_typ <- median(v)

  if (ncol(X) == 0) {
    res_base <- y
  } else {
    fit <- stats::lm.wfit(x = cbind(1, X), y = y, w = 1 / v)
    beta_hat <- fit$coefficients[-1]
    beta_hat[is.na(beta_hat)] <- 0
    res_base <- y - as.numeric(X %*% beta_hat)
  }

  mu_rob <- median(res_base)
  res <- res_base - mu_rob
  var_total_rob <- mad(res_base)^2

  tau2_rob <- max(0, var_total_rob - v_typ)
  shrinkage <- tau2_rob / (tau2_rob + v)
  theta_tilde <- mu_rob + shrinkage * res

  pad <- sqrt(tau2_rob)
  if (pad == 0) pad <- sqrt(v_typ)

  lower_bound <- min(theta_tilde) - 1.5 * pad
  upper_bound <- max(theta_tilde) + 1.5 * pad

  return(seq(lower_bound, upper_bound, length.out = K))
}

#' @noRd
pgm_objective <- function(eta, y, v, X, mu, tau_c, lambda_alpha, DtD, k0) {
  N <- length(y)
  K <- length(mu)
  p <- ncol(X)

  alpha <- eta[1:(K - 1)]
  beta <- if (p > 0) eta[K:(K - 1 + p)] else NULL

  alpha_star <- pad_zero(alpha, k0)

  max_alpha <- max(alpha_star)
  exp_alpha <- exp(alpha_star - max_alpha)
  w <- exp_alpha / sum(exp_alpha)

  mu_mat <- matrix(mu, nrow = N, ncol = K, byrow = TRUE)
  if (p > 0) mu_mat <- mu_mat + as.numeric(X %*% beta)

  v_marg <- v + tau_c^2
  sd_mat <- matrix(sqrt(v_marg), nrow = N, ncol = K)
  Phi <- dnorm(y, mean = mu_mat, sd = sd_mat)
  f_y <- as.numeric(Phi %*% w)
  l <- sum(log(f_y))

  pen <- (lambda_alpha / 2) * as.numeric(crossprod(alpha_star, DtD %*% alpha_star))
  Q <- (Phi / f_y) * rep(w, each = N)
  s_tilde <- colSums(Q) - N * w

  grad_pen_tilde <- as.numeric(lambda_alpha * (DtD %*% alpha_star))
  grad_alpha <- (s_tilde - grad_pen_tilde)[-k0]

  if (p > 0) {
    resid_scaled <- (y - mu_mat) / v_marg
    expected_resid <- rowSums(Q * resid_scaled)
    grad_beta <- as.numeric(crossprod(X, expected_resid))
    grad_obj <- c(grad_alpha, grad_beta)
  } else {
    grad_obj <- grad_alpha
  }

  J_tilde <- crossprod(Q) - (N * tcrossprod(w)) - diag(s_tilde)
  H_pen_tilde <- lambda_alpha * DtD
  H_alphaalpha <- -(J_tilde + H_pen_tilde)[-k0, -k0, drop = FALSE]

  if (p > 0) {
    W_vec <- rowSums(Q / v_marg) - rowSums(Q * resid_scaled^2) + expected_resid^2
    H_betabeta <- -crossprod(X, W_vec * X)
    H_alphabeta <- crossprod(Q * (resid_scaled - expected_resid), X)
    H_alphabeta <- H_alphabeta[-k0, , drop = FALSE]
    hess_obj <- rbind(cbind(H_alphaalpha, H_alphabeta), cbind(t(H_alphabeta), H_betabeta))
  } else {
    hess_obj <- H_alphaalpha
  }

  return(list(value = l - pen, gradient = grad_obj, hessian = hess_obj))
}

#' @noRd
fit_pgm_path <- function(y, v, X, lambda_alpha_grid, choose_lambda, alpha_init, beta_init, mu, tau_c, DtD, k0, ...) {
  N <- length(y)
  K <- length(mu)
  p <- ncol(X)

  trust_args <- list(...)
  if (is.null(trust_args$rinit)) trust_args$rinit <- 1
  if (is.null(trust_args$rmax)) trust_args$rmax <- 100

  current_parinit <- if (p > 0) c(alpha_init, beta_init) else alpha_init
  results_list <- vector("list", length(lambda_alpha_grid))

  for (i in seq_along(lambda_alpha_grid)) {
    lambda_alpha <- lambda_alpha_grid[i]
    call_args <- c(
      list(
        objfun = pgm_objective, parinit = current_parinit, minimize = FALSE,
        y = y, v = v, X = X, mu = mu, tau_c = tau_c, lambda_alpha = lambda_alpha, DtD = DtD, k0 = k0
      ),
      trust_args
    )

    opt <- tryCatch(do.call(trust::trust, call_args), error = function(e) NULL)

    if (is.null(opt) || !opt$converged) {
      results_list[[i]] <- list(ic = Inf, alpha = rep(NA, length(alpha_init)), beta = rep(NA, p), edf = NA, Jp_inv = NULL, J = NULL)
      next
    }

    alpha <- opt$argument[1:(K - 1)]
    beta <- if (p > 0) opt$argument[K:(K - 1 + p)] else NULL
    current_parinit <- opt$argument

    Jp <- -opt$hessian
    Jp_inv <- tryCatch(chol2inv(chol(Jp)), error = function(e) NULL)

    if (is.null(Jp_inv)) {
      results_list[[i]] <- list(ic = Inf, alpha = rep(NA, length(alpha_init)), beta = rep(NA, p), edf = NA, Jp_inv = NULL, J = NULL)
      next
    }

    alpha_star <- pad_zero(alpha, k0)

    J <- Jp
    idx <- -k0
    J[1:(K - 1), 1:(K - 1)] <- Jp[1:(K - 1), 1:(K - 1)] - (lambda_alpha * DtD)[idx, idx, drop = FALSE]

    edf <- sum(diag(J %*% Jp_inv))
    pen <- (lambda_alpha / 2) * as.numeric(crossprod(alpha_star, DtD %*% alpha_star))

    ic_val <- if (choose_lambda == "bic") {
      -2 * (opt$value + pen) + log(N) * edf
    } else {
      -2 * (opt$value + pen) + 2 * edf
    }

    results_list[[i]] <- list(ic = ic_val, alpha = alpha, beta = beta, edf = edf, Jp_inv = Jp_inv, J = J)
  }

  return(results_list)
}

#' @noRd
calc_V_robust <- function(y, v, X, group, mu, tau_c, alpha, beta, Jp_inv, k0) {
  N <- length(y)
  K <- length(mu)
  p <- ncol(X)

  alpha_star <- pad_zero(alpha, k0)
  max_alpha <- max(alpha_star)
  exp_alpha <- exp(alpha_star - max_alpha)
  w <- exp_alpha / sum(exp_alpha)

  mu_mat <- matrix(mu, nrow = N, ncol = K, byrow = TRUE)
  if (p > 0) mu_mat <- mu_mat + as.numeric(X %*% beta)

  v_marg <- v + tau_c^2
  sd_mat <- matrix(sqrt(v_marg), nrow = N, ncol = K)
  Phi <- dnorm(y, mean = mu_mat, sd = sd_mat)
  f_y <- as.numeric(Phi %*% w)
  Q <- (Phi / f_y) * rep(w, each = N)

  idx <- -k0
  S_alpha <- Q[, idx, drop = FALSE] - rep(w[idx], each = N)

  if (p > 0) {
    resid_scaled <- (y - mu_mat) / v_marg
    expected_resid <- rowSums(Q * resid_scaled)
    S_beta <- X * expected_resid
    S <- cbind(S_alpha, S_beta)
  } else {
    S <- S_alpha
  }

  S_group <- rowsum(S, group)
  B <- crossprod(S_group)

  return(Jp_inv %*% B %*% Jp_inv)
}

#' @noRd
pgm_cond_objective <- function(eta, y, v, z, mu, tau_c, lambda_alpha, lambda_gamma, DtD, k0) {
  N <- length(y)
  K <- length(mu)

  alpha <- eta[1:(K - 1)]
  gamma <- eta[K:(2 * K - 2)]

  alpha_star <- pad_zero(alpha, k0)
  gamma_star <- pad_zero(gamma, k0)

  eta_mat <- matrix(alpha_star, nrow = N, ncol = K, byrow = TRUE) + z %*% t(gamma_star)
  max_eta <- apply(eta_mat, 1, max)
  exp_eta <- exp(eta_mat - max_eta)
  W <- exp_eta / rowSums(exp_eta)

  v_marg <- v + tau_c^2
  sd_mat <- matrix(sqrt(v_marg), nrow = N, ncol = K)
  mu_mat <- matrix(mu, nrow = N, ncol = K, byrow = TRUE)
  Phi <- dnorm(y, mean = mu_mat, sd = sd_mat)

  f_y <- rowSums(Phi * W)
  l <- sum(log(f_y))

  pen_alpha <- (lambda_alpha / 2) * as.numeric(crossprod(alpha_star, DtD %*% alpha_star))
  pen_gamma <- (lambda_gamma / 2) * as.numeric(crossprod(gamma_star, DtD %*% gamma_star))
  pen <- pen_alpha + pen_gamma

  Q <- (Phi / f_y) * W
  diff_mat <- Q - W

  s_tilde_alpha <- colSums(diff_mat)
  s_tilde_gamma <- as.numeric(crossprod(diff_mat, z))

  grad_pen_alpha_tilde <- as.numeric(lambda_alpha * (DtD %*% alpha_star))
  grad_pen_gamma_tilde <- as.numeric(lambda_gamma * (DtD %*% gamma_star))

  grad_alpha <- (s_tilde_alpha - grad_pen_alpha_tilde)[-k0]
  grad_gamma <- (s_tilde_gamma - grad_pen_gamma_tilde)[-k0]
  grad_obj <- c(grad_alpha, grad_gamma)

  Q_z <- Q * z
  W_z <- W * z

  H_alphaalpha <- diag(s_tilde_alpha) - (crossprod(Q) - crossprod(W)) - lambda_alpha * DtD
  H_alphagamma <- diag(s_tilde_gamma) - (crossprod(Q, Q_z) - crossprod(W, W_z))
  H_gammagamma <- diag(as.numeric(crossprod(diff_mat, z^2))) - (crossprod(Q_z) - crossprod(W_z)) - lambda_gamma * DtD

  idx <- -k0
  hess_obj <- rbind(
    cbind(H_alphaalpha[idx, idx, drop = FALSE], H_alphagamma[idx, idx, drop = FALSE]),
    cbind(t(H_alphagamma[idx, idx, drop = FALSE]), H_gammagamma[idx, idx, drop = FALSE])
  )

  return(list(value = l - pen, gradient = grad_obj, hessian = hess_obj))
}

#' @noRd
fit_pgm_cond_path <- function(y, v, z, lambda_alpha_grid, lambda_gamma_grid, choose_lambda, eta_init, mu, tau_c, DtD, k0, ...) {
  N <- length(y)
  K <- length(mu)

  trust_args <- list(...)
  if (is.null(trust_args$rinit)) trust_args$rinit <- 1
  if (is.null(trust_args$rmax)) trust_args$rmax <- 100

  grid <- expand.grid(lambda_alpha = lambda_alpha_grid, lambda_gamma = lambda_gamma_grid)
  n_models <- nrow(grid)
  results_list <- vector("list", n_models)
  current_parinit <- eta_init

  for (i in seq_len(n_models)) {
    lambda_alpha <- grid$lambda_alpha[i]
    lambda_gamma <- grid$lambda_gamma[i]

    call_args <- c(
      list(
        objfun = pgm_cond_objective, parinit = current_parinit, minimize = FALSE,
        y = y, v = v, z = z, mu = mu, tau_c = tau_c,
        lambda_alpha = lambda_alpha, lambda_gamma = lambda_gamma, DtD = DtD, k0 = k0
      ),
      trust_args
    )

    opt <- tryCatch(do.call(trust::trust, call_args), error = function(e) NULL)

    if (is.null(opt) || !opt$converged) {
      results_list[[i]] <- list(ic = Inf, alpha = rep(NA, K - 1), gamma = rep(NA, K - 1), edf = NA, Jp_inv = NULL, J = NULL, lambda_alpha = lambda_alpha, lambda_gamma = lambda_gamma)
      next
    }

    alpha <- opt$argument[1:(K - 1)]
    gamma <- opt$argument[K:(2 * K - 2)]
    current_parinit <- opt$argument

    Jp <- -opt$hessian
    Jp_inv <- tryCatch(chol2inv(chol(Jp)), error = function(e) NULL)

    if (is.null(Jp_inv)) {
      results_list[[i]] <- list(ic = Inf, alpha = rep(NA, K - 1), gamma = rep(NA, K - 1), edf = NA, Jp_inv = NULL, J = NULL, lambda_alpha = lambda_alpha, lambda_gamma = lambda_gamma)
      next
    }

    alpha_star <- pad_zero(alpha, k0)
    gamma_star <- pad_zero(gamma, k0)

    J <- Jp
    idx <- -k0
    J[1:(K - 1), 1:(K - 1)] <- Jp[1:(K - 1), 1:(K - 1)] - (lambda_alpha * DtD)[idx, idx, drop = FALSE]
    J[K:(2 * K - 2), K:(2 * K - 2)] <- Jp[K:(2 * K - 2), K:(2 * K - 2)] - (lambda_gamma * DtD)[idx, idx, drop = FALSE]

    edf <- sum(diag(J %*% Jp_inv))
    pen_alpha <- (lambda_alpha / 2) * as.numeric(crossprod(alpha_star, DtD %*% alpha_star))
    pen_gamma <- (lambda_gamma / 2) * as.numeric(crossprod(gamma_star, DtD %*% gamma_star))
    pen <- pen_alpha + pen_gamma

    ic_val <- if (choose_lambda == "bic") {
      -2 * (opt$value + pen) + log(N) * edf
    } else {
      -2 * (opt$value + pen) + 2 * edf
    }

    results_list[[i]] <- list(ic = ic_val, alpha = alpha, gamma = gamma, edf = edf, Jp_inv = Jp_inv, J = J, lambda_alpha = lambda_alpha, lambda_gamma = lambda_gamma)
  }

  return(results_list)
}

#' @noRd
calc_V_cond_robust <- function(y, v, z, group, mu, tau_c, alpha, gamma, Jp_inv, k0) {
  N <- length(y)
  K <- length(mu)

  alpha_star <- pad_zero(alpha, k0)
  gamma_star <- pad_zero(gamma, k0)

  eta_mat <- matrix(alpha_star, nrow = N, ncol = K, byrow = TRUE) + z %*% t(gamma_star)
  max_eta <- apply(eta_mat, 1, max)
  exp_eta <- exp(eta_mat - max_eta)
  W <- exp_eta / rowSums(exp_eta)

  v_marg <- v + tau_c^2
  sd_mat <- matrix(sqrt(v_marg), nrow = N, ncol = K)
  mu_mat <- matrix(mu, nrow = N, ncol = K, byrow = TRUE)
  Phi <- dnorm(y, mean = mu_mat, sd = sd_mat)

  f_y <- rowSums(Phi * W)
  Q <- (Phi / f_y) * W
  diff_mat <- Q - W

  idx <- -k0
  S_alpha <- diff_mat[, idx, drop = FALSE]
  S_gamma <- (diff_mat * z)[, idx, drop = FALSE]
  S <- cbind(S_alpha, S_gamma)

  S_group <- rowsum(S, group)
  B <- crossprod(S_group)

  return(Jp_inv %*% B %*% Jp_inv)
}
