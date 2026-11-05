# ----------------
# Define functions
# ----------------

# Function for calculating the angle between two vectors
func_angulardiff <- function(x,y) {
  dot.prod <- dot(x,y)
  norm.x <- Norm(x,p = 2)
  norm.y <- Norm(y,p = 2)
  # clamp to [-1, 1] to avoid NaN from rounding errors on near-colinear vectors
  theta <- acos(max(min(dot.prod / (norm.x * norm.y), 1), -1))
  return(theta)
}

# Function for calculating the bootstrapped statistics
reg.boot <- function(regdat,index) {
  xv <- regdat[,1][index]
  yv <- regdat[,2][index]
  model <- lm(yv~xv)
  PI.boot <-
    -summary(model)[[4]][1] / summary(model)[[4]][2]
  #coef(model)
}

reg.boot_anch <- function(regdat,index) {
  xv <- regdat[,1][index]
  yv <- regdat[,2][index]
  wv <- regdat[,3][index]
  model <- lm(yv~xv - 1,weights=wv)
  PI.boot <- 1 / coef(model)
}


mytheme <- theme_bw() +
  theme(text = element_text(colour="black"),
        #panel.border = element_blank(),
        axis.title = element_text(size = rel(1.3)),
        axis.text = element_text(colour="red",size = rel(1.1)))

###########################################################################
# Graphic style
###########################################################################
status_shapes <- c("Good" = 16, "Bad (selection)" = 15, "Outlier (Cook)" = 17)
status_colors <- c("Good" = "steelblue", "Bad (selection)" = "darkorange", "Outlier (Cook)" = "red")

# same order of the legend entries for shape and colour (only levels present are shown)
# larger legend (title and entries) for the figures that have one
legend_big <- theme(legend.title = element_text(size = 16),
                    legend.text  = element_text(size = 14))

status_limits <- function(x) intersect(names(status_shapes), x)

msp_status <- function(quality) {
  factor(unname(c(Good = "Good", Bad = "Bad (selection)", Outlier = "Outlier (Cook)")[as.character(quality)]),
         levels = names(status_shapes))
}

pi_label <- function(res, note = "") {
  sprintf("PI = %.1f [%.1f - %.1f] (\u00b5T)%s", res$pi_val, res$l95, res$u95, note)
}

# Q vs field plot. pts: data.frame with x, y, quality (and dQ if errbars = TRUE)
msp_plot <- function(pts, line, band, res, ylab, anchor = TRUE, errbars = FALSE, note = "") {
  pts$status <- msp_status(pts$quality)
  p <- ggplot() +
    # 95% confidence interval of the PI, same colour as the bootstrap density
    annotate("rect", xmin = res$l95, xmax = res$u95, ymin = -Inf, ymax = Inf,
             fill = "steelblue", alpha = 0.4) +
    geom_hline(yintercept = 0, linewidth = 1, colour = "black") +
    geom_polygon(data = band, aes(x = Field, y = y_err), fill = "pink", alpha = 0.3) +
    geom_line(data = line, aes(x = Field, y = fit), colour = "darkred")
  if (errbars) {
    p <- p + geom_errorbar(data = subset(pts, quality != "Bad"),
                           aes(x = x, ymin = y - dQ, ymax = y + dQ, colour = status),
                           width = 0, linewidth = 0.5, alpha = 0.6)
  }
  p <- p +
    geom_point(data = pts, aes(x = x, y = y, shape = status, colour = status), size = 2.5) +
    scale_shape_manual(values = status_shapes, limits = status_limits) +
    scale_colour_manual(values = status_colors, limits = status_limits)
  if (anchor) {
    p <- p + geom_point(aes(x = 0, y = -1), colour = "black", shape = 4, size = 3)
  }
  p +
    # paleointensity (Q = 0): solid line (the 95% CI is the coloured band)
    geom_vline(xintercept = res$pi_val, colour = "darkred", linewidth = 1) +
    # annotate("text", x = -Inf, y = Inf, label = pi_label(res, note),
    #          hjust = -0.05, vjust = 1.5, size = 5, colour = "grey20") +
    labs(x = expression(mu[0] ~ H[lab] ~ "(" * mu * "T)"), y = ylab, shape = "Status", colour = "Status") +
    coord_cartesian(xlim = c(0, max(pts$x) + 3)) +
    theme_minimal() +
    theme(axis.title = element_text(size = 20),
          axis.text  = element_text(size = 18)) +
    legend_big
}

# Cook's distance plot (points above the threshold are the outliers)
cook_plot <- function(q, th) {
  q <- q[!is.na(q$Cook), ]
  q$status <- msp_status(q$quality)
  ggplot(q, aes(x = x, y = Cook)) +
    geom_hline(yintercept = th, linetype = "dashed", colour = "grey40") +
    geom_point(aes(shape = status, colour = status), size = 2.5) +
    scale_shape_manual(values = status_shapes, limits = status_limits) +
    scale_colour_manual(values = status_colors, limits = status_limits) +
    labs(x = "Laboratory Field (\u00b5T)", y = "Cook's distance", shape = "Status", colour = "Status") +
    theme_minimal() +
    theme(axis.title = element_text(size = 20),
          axis.text  = element_text(size = 18)) +
    legend_big
}

# Bootstrap distribution of the paleointensity, with the retained PI (solid)
# and its 95% confidence interval (dashed).
boot_plot <- function(dat, res, note = "", xr = NULL) {
  b   <- as.numeric(dat$b.val)
  b   <- b[is.finite(b)]
  dns <- density(b, from = min(b), to = max(b), n = 512)   # same estimate as geom_density
  curve <- data.frame(x = dns$x, y = dns$y)
  lo <- max(res$l95, min(b))
  hi <- min(res$u95, max(b))
  p <- ggplot() +
    geom_line(data = curve, aes(x = x, y = y), colour = "steelblue4", linewidth = 0.8)
  if (lo < hi) {   # only the 95% confidence interval is coloured
    xi <- seq(lo, hi, length.out = 400)
    ci <- data.frame(x = xi, y = approx(dns$x, dns$y, xout = xi)$y)
    p <- p + geom_ribbon(data = ci, aes(x = x, ymin = 0, ymax = y), fill = "steelblue", alpha = 0.4)
  }
  p <- p +
    geom_vline(xintercept = res$pi_val, colour = "darkred", linewidth = 1) +
    geom_vline(xintercept = c(res$l95, res$u95), colour = "grey40", linetype = "dashed") +
    # annotate("text", x = -Inf, y = Inf, label = pi_label(res, note),
    #          hjust = -0.05, vjust = 1.5, size = 5, colour = "grey20") +
    labs(x = expression(mu[0] ~ H[anc] ~ "(" * mu * "T)"), y = "Density") +
    theme_minimal() +
    theme(axis.title = element_text(size = 20),
          axis.text  = element_text(size = 18))
  # x-range chosen with the slider: exactly [xmin, xmax] (density computed on all values)
  if (!is.null(xr)) p <- p + scale_x_continuous(expand = c(0, 0)) + coord_cartesian(xlim = xr)
  p
}

conf_level = c(0.68,0.95); # confidence interval
#m4m3threshold <- -0.1; # (m4-m3) / m4 (m4 should be greater than m3) 
#m1m3threshold <- -0.05; # (m1-m3) / m1 (m3 should be lower than m1)
m4m3threshold <- -1; # (m4-m3) / m4 (m4 should be greater than m3) 
m1m3threshold <- -1; # (m1-m3) / m1 (m3 should be lower than m1)

###########################################################################
# MSP-DSC : strategies for the alpha parameter (cases 1 to 4)
#   case 1 : constant alpha (slider, default 0.5)
#   case 2 : single alpha common to all specimens, optimized over [0, 0.8]
#   case 3 : per-specimen alpha, optimized, same bounds [0, 0.8] for all
#   case 4 : per-specimen alpha, optimized, individual bounds [a_min, a_max]
#            (last two columns of the input file)
#   case 5 : per-specimen alpha, optimized in [0, 0.8], with a Gaussian prior
#            centred on the middle of [a_min, a_max] of each specimen and
#            of standard deviation sigma_i = (a_max - a_min) / k.
#            lambda_i = 1 / sigma_i^2 is fixed by the data file (no LOOCV).
###########################################################################

alpha_lo  <- 0     # lower bound of alpha (cases 2 and 3)
alpha_hi  <- 0.8   # upper bound of alpha (cases 2 and 3)
alpha_case1 <- 0.5 # constant alpha of case 1
alpha_pop <- 0.5   # ridge shrinkage target (case 3) and fallback value
lambda_grid <- round(seq(0.5, 2, by = 0.05), 2)   # ridge strengths searched by leave-one-out CV
alpha_dQ  <- 0.5   # alpha used to compute the error bars dQ in cases 2-5
sigma_k_default <- 4     # case 5: the interval [a_min, a_max] is k standard deviations wide (k = 4: +/- 2 sigma)
sigma_min       <- 0.05  # case 5: lower limit of sigma_i (avoids an infinite ridge strength for very narrow intervals)

# Case 5 : ridge centre and strength of each specimen from its interval [a_min, a_max]
case5_prior <- function(a_min, a_max, k = sigma_k_default) {
  sigma <- pmax((a_max - a_min) / k, sigma_min)
  list(ctr = (a_min + a_max) / 2, lam = 1 / sigma^2)
}

# Error bar on Qdsc (Eq. 19-22, Fabian & Leonhardt 2010, Appendix A.2)
calc_dQ <- function(m0, m1, m3, m4, m2c, alpha) {
  epsi   <- abs((m1 - m4) / m1)                                  # Eq. 19
  D      <- 2 * m0 - m1 - m2c
  dQalt2 <- epsi^2 * 4 / D^4 * (                                 # Eq. 20
    alpha^2 * m3^2 * D^2 +
      m2c^2 * (m0 - (1 + alpha) * m1 + alpha * m3)^2 +
      m1^2 * ((1 + 2 * alpha) * m0 - (1 + alpha) * m2c - alpha * m3)^2
  )
  dQds2  <- (1 / 3) * ((m3 - m1) / D)^2                          # Eq. 21
  sqrt(dQalt2 + dQds2)                                           # Eq. 22
}

# Qdsc as a function of alpha (a per-specimen vector is allowed)
compute_Qdsc <- function(m0, m1, m2, m3, alpha) {
  2 * ((1 + alpha) * m1 - m0 - alpha * m3) / (2 * m0 - m1 - m2)
}

# Case 2 : for a given global alpha, profiled weighted-least-squares slope
# (regression anchored at (0,-1)) and weighted residual sum of squares.
# data_sub needs the columns Field, m0, m1, m3, m2_corr, Qfc, dQ.
profile_slope <- function(alpha, data_sub) {
  c_i    <- 2 * (data_sub$m1 - data_sub$m3) / (2 * data_sub$m0 - data_sub$m1 - data_sub$m2_corr)
  Qdsc_i <- data_sub$Qfc + alpha * c_i
  w      <- 1 / data_sub$dQ^2
  
  slope <- sum(w * data_sub$Field * (Qdsc_i + 1)) / sum(w * data_sub$Field^2)
  resid <- Qdsc_i - (slope * data_sub$Field - 1)
  list(slope = slope, SSR = sum(w * resid^2))
}

# Ridge target of each alpha_i : a scalar (case 3), "mid" = midpoint of the
# specimen's own interval [a_min, a_max] (case 4), or a vector.
resolve_alpha_center <- function(alpha_center, a_min, a_max, a_ctr = NULL) {
  n <- length(a_min)
  if (is.character(alpha_center)) {
    if (!alpha_center %in% c("mid", "col")) stop("alpha_center must be numeric, \"mid\" or \"col\"")
    # "mid": middle of [a_min, a_max]; "col": column a_ctr of the data (case 5)
    ctr <- if (alpha_center == "col") a_ctr else (a_min + a_max) / 2
  } else if (length(alpha_center) == 1) {
    ctr <- rep(alpha_center, n)
  } else {
    ctr <- alpha_center
  }
  ctr[is.na(ctr)] <- alpha_pop
  ctr
}

# Cases 3, 4 and 5 : slope + per-specimen alpha by block coordinate descent
# (ridge shrinkage toward the centre, alpha_i clipped to [a_min_i, a_max_i]).
# data_sub needs the columns Field, m0, m1, m3, m2_corr, Qfc, dQ, a_min, a_max.
# Case 5 : data_sub also carries a_ctr (ridge centre) and lam_i (ridge strength of
# each specimen) and alpha_center = "col". They travel with the rows, so they stay
# aligned in the bootstrap samples.
fit_block_alpha <- function(data_sub, lambda, alpha_center = alpha_pop,
                            max_iter = 50, tol = 1e-8) {
  c_i   <- 2 * (data_sub$m1 - data_sub$m3) / (2 * data_sub$m0 - data_sub$m1 - data_sub$m2_corr)
  Qfc_i <- data_sub$Qfc
  w     <- 1 / data_sub$dQ^2
  a_min <- data_sub$a_min
  a_max <- data_sub$a_max
  if ("lam_i" %in% names(data_sub)) lambda <- data_sub$lam_i   # case 5: one lambda per specimen
  
  center <- resolve_alpha_center(alpha_center, a_min, a_max, data_sub$a_ctr)
  alpha  <- pmin(pmax(center, a_min), a_max)
  
  for (iter in seq_len(max_iter)) {
    Qdsc_i <- Qfc_i + alpha * c_i
    slope  <- sum(w * data_sub$Field * (Qdsc_i + 1)) / sum(w * data_sub$Field^2)
    
    k_i       <- Qfc_i - slope * data_sub$Field + 1
    alpha_new <- (center * lambda - w * c_i * k_i) / (w * c_i^2 + lambda)
    alpha_new <- pmin(pmax(alpha_new, a_min), a_max)
    
    if (max(abs(alpha_new - alpha)) < tol) break
    alpha <- alpha_new
  }
  list(slope = slope, alpha = alpha, c_i = c_i)
}

# Leave-one-out cross-validation of the ridge strength lambda (cases 3 and 4)
# (as in paleointensity_analysis_V3.R, section 6bis).
# d = selected specimens (with a_min, a_max); alpha_center = alpha_pop (case 3)
# or "mid" (case 4). Returns lambda and its LOOCV error (weighted squared error).
select_lambda_loocv <- function(d, lambda_grid, alpha_center = alpha_pop) {
  n      <- nrow(d)
  center <- resolve_alpha_center(alpha_center, d$a_min, d$a_max)
  c_all  <- 2 * (d$m1 - d$m3) / (2 * d$m0 - d$m1 - d$m2_corr)
  w_all  <- 1 / d$dQ^2
  
  errors <- numeric(length(lambda_grid))
  for (li in seq_along(lambda_grid)) {
    lam    <- lambda_grid[li]
    sq_err <- numeric(n)
    for (i in seq_len(n)) {
      fit_i <- fit_block_alpha(d[-i, , drop = FALSE], lam, alpha_center = center[-i])
      
      k_i   <- d$Qfc[i] - fit_i$slope * d$Field[i] + 1
      a_i   <- (center[i] * lam - w_all[i] * c_all[i] * k_i) / (w_all[i] * c_all[i]^2 + lam)
      a_i   <- min(max(a_i, d$a_min[i]), d$a_max[i])
      
      resid <- d$Qfc[i] + a_i * c_all[i] - (fit_i$slope * d$Field[i] - 1)
      sq_err[i] <- w_all[i] * resid^2
    }
    errors[li] <- sum(sq_err)
  }
  data.frame(lambda = lambda_grid, loocv_error = errors)
}

# Case-resampling bootstrap. fit_fun(data_resampled) returns the fitted slope.
# Returns the bootstrapped paleointensities (1/slope), their 68% and 95%
# percentile intervals, and a pointwise 95% band for the line slope*Field - 1.
bootstrap_pi <- function(d, fit_fun, field_grid, B, seed = 123, progress = NULL) {
  n <- nrow(d)
  boot_slope <- rep(NA_real_, B)
  
  withr::with_seed(seed, {
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      boot_slope[b] <- tryCatch(fit_fun(d[idx, , drop = FALSE]),
                                error = function(e) NA_real_)
      if (!is.null(progress) && b %% max(1, B %/% 100) == 0) progress(b / B)
    }
  })
  
  boot_pi  <- 1 / boot_slope
  ok       <- is.finite(boot_pi)
  fit_grid <- outer(boot_slope, field_grid, function(s, f) s * f - 1)
  
  list(
    pi         = boot_pi[ok],
    ci68       = quantile(boot_pi[ok], c(0.16, 0.84), names = FALSE),
    ci95       = quantile(boot_pi[ok], c(0.025, 0.975), names = FALSE),
    band_lower = apply(fit_grid, 2, quantile, probs = 0.025, na.rm = TRUE, names = FALSE),
    band_upper = apply(fit_grid, 2, quantile, probs = 0.975, na.rm = TRUE, names = FALSE),
    n_failed   = sum(!ok)
  )
}

# Stage 1 : fit on ALL selected specimens (including those flagged as outliers
#           under the constant-alpha fit).
# Stage 2 : new Cook's-distance outlier detection on the optimized fit
#           (threshold 3 x mean), then final fit on the remaining specimens.
# Both functions take d = data.frame of the selected specimens (with an `id`
# column = row index in the input file) and return the same list structure.

# Case 2 : single alpha common to all specimens, optimized over [lo, hi]
run_global_alpha <- function(d, lo, hi, B, field_grid, progress = NULL) {
  opt1 <- optimize(function(a) profile_slope(a, d)$SSR, interval = c(lo, hi))
  d$Q1 <- compute_Qdsc(d$m0, d$m1, d$m2_corr, d$m3, opt1$minimum)
  
  mod1 <- lm(I(Q1 + 1) ~ Field - 1, data = d, weights = 1 / dQ^2)
  d$cook    <- cooks.distance(mod1)
  d$outlier <- d$cook > 3 * mean(d$cook, na.rm = TRUE)
  d$outlier[is.na(d$outlier)] <- FALSE
  
  clean <- d[!d$outlier, , drop = FALSE]
  opt2  <- optimize(function(a) profile_slope(a, clean)$SSR, interval = c(lo, hi))
  alpha_final <- opt2$minimum
  slope       <- profile_slope(alpha_final, clean)$slope
  
  fit_fun <- function(db) {
    ob <- optimize(function(a) profile_slope(a, db)$SSR, interval = c(lo, hi))
    profile_slope(ob$minimum, db)$slope
  }
  boot <- bootstrap_pi(clean, fit_fun, field_grid, B, progress = progress)
  
  clean$Qf <- compute_Qdsc(clean$m0, clean$m1, clean$m2_corr, clean$m3, alpha_final)
  r2 <- summary(lm(I(Qf + 1) ~ Field - 1, data = clean, weights = 1 / dQ^2))$r.squared
  
  list(slope = slope, alpha_global = alpha_final,
       alpha_clean = rep(alpha_final, nrow(clean)), clean_ids = clean$id,
       cook = d$cook, outlier_ids = d$id[d$outlier], boot = boot, r2 = r2)
}

# Cases 3, 4 and 5 : one alpha per specimen (ridge-shrunk), bounds taken from the
# columns a_min / a_max of d (case 5: d also has a_ctr and lam_i, lambda is ignored).
run_perspecimen_alpha <- function(d, lambda, alpha_center, B, field_grid, progress = NULL) {
  f1   <- fit_block_alpha(d, lambda, alpha_center)
  d$Q1 <- d$Qfc + f1$alpha * f1$c_i
  
  mod1 <- lm(I(Q1 + 1) ~ Field - 1, data = d, weights = 1 / dQ^2)
  d$cook    <- cooks.distance(mod1)
  d$outlier <- d$cook > 3 * mean(d$cook, na.rm = TRUE)
  d$outlier[is.na(d$outlier)] <- FALSE
  
  clean <- d[!d$outlier, , drop = FALSE]
  f2    <- fit_block_alpha(clean, lambda, alpha_center)
  
  fit_fun <- function(db) fit_block_alpha(db, lambda, alpha_center)$slope
  boot <- bootstrap_pi(clean, fit_fun, field_grid, B, progress = progress)
  
  clean$Qf <- clean$Qfc + f2$alpha * f2$c_i
  r2 <- summary(lm(I(Qf + 1) ~ Field - 1, data = clean, weights = 1 / dQ^2))$r.squared
  
  list(slope = f2$slope, alpha_global = NA_real_,
       alpha_clean = f2$alpha, clean_ids = clean$id,
       cook = d$cook, outlier_ids = d$id[d$outlier], boot = boot, r2 = r2)
}



# ---------------------------------------------------------------------------
# Comparison of the five MSP-DSC cases : paired (case-resampling) bootstrap
# ---------------------------------------------------------------------------
# The five cases are computed on the same specimens, so their paleointensities
# are strongly correlated and cannot be compared with a test for independent
# estimates. Here the SAME resampled specimens are used for the five cases at
# each replicate, so that the differences PI_k - PI_j keep the correlation.
#
# n     : number of specimens of the common set (specimens passing the selection)
# cases : list (one element per case) of list(data, keep, fit):
#           data = data.frame with n rows (all the columns the case needs)
#           keep = logical(n), specimens retained by the case in its final fit
#                  (the exclusion of the influential points is fixed, as in the
#                  bootstrap of each case)
#           fit  = function(db) returning the slope of the anchored regression
# Returns a B x K matrix of bootstrapped paleointensities (NA if a fit fails or
# if the slope is not positive).
paired_bootstrap_cases <- function(n, cases, B, seed = 321, progress = NULL,
                                   min_n = 4) {
  K <- length(cases)
  slopes <- matrix(NA_real_, nrow = B, ncol = K)
  withr::with_seed(seed, {
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      for (j in seq_len(K)) {
        ij <- idx[cases[[j]]$keep[idx]]
        if (length(ij) >= min_n) {
          slopes[b, j] <- tryCatch(cases[[j]]$fit(cases[[j]]$data[ij, , drop = FALSE]),
                                   error = function(e) NA_real_)
        }
      }
      if (!is.null(progress) && b %% max(1, B %/% 100) == 0) progress(b / B)
    }
  })
  slopes[!is.finite(slopes) | slopes <= 0] <- NA_real_
  1 / slopes
}

# Differences between cases (later case minus earlier case) with percentile
# intervals and bootstrap p-values (two-sided), Holm-adjusted for the number of
# comparisons. pi_hat : paleointensity of each case (full data set).
compare_cases <- function(pi_boot, pi_hat, labels = paste("Case", seq_along(pi_hat)),
                          conf = 0.95) {
  K <- ncol(pi_boot)
  rows <- list()
  for (j in seq_len(K - 1)) for (k in (j + 1):K) {
    d  <- pi_boot[, k] - pi_boot[, j]
    d  <- d[is.finite(d)]
    nb <- length(d)
    ci <- quantile(d, c((1 - conf) / 2, 1 - (1 - conf) / 2), names = FALSE)
    p  <- min(1, 2 * (1 + min(sum(d <= 0), sum(d >= 0))) / (nb + 1))
    rows[[length(rows) + 1]] <- data.frame(
      Comparison = paste(labels[k], "-", labels[j]),
      Difference = pi_hat[k] - pi_hat[j],
      CI_low = ci[1], CI_up = ci[2], p_value = p, n_boot = nb,
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$p_Holm <- p.adjust(out$p_value, method = "holm")
  out
}


# ---------------------------------------------------------------------------
# Plot of the alpha values of two cases, specimen by specimen (1:1 line, colour =
# class of laboratory field). cmp : data.frame with the columns B (field), Sample
# and the two alpha columns named by xcol and ycol.
# ---------------------------------------------------------------------------
alpha_compare_plot <- function(cmp, xcol, ycol, xlab, ylab) {
  # classes of laboratory field, step 10 uT
  lo  <- floor(min(cmp$B) / 10) * 10
  hi  <- max(ceiling(max(cmp$B) / 10) * 10, lo + 10)
  brk <- seq(lo, hi, by = 10)
  cmp$Field_class <- cut(cmp$B, breaks = brk, include.lowest = TRUE, right = FALSE,
                         labels = paste0(head(brk, -1), "-", tail(brk, -1)))

  ggplot(cmp, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = Field_class), size = 3) +
    ggrepel::geom_text_repel(aes(label = Sample), size = 3.5,
                             box.padding = 0.5, point.padding = 0.3,
                             min.segment.length = 0, segment.colour = "grey60",
                             segment.size = 0.3, max.overlaps = Inf, seed = 42) +
    scale_colour_viridis_d(option = "viridis", name = expression(H[lab]), drop = FALSE) +
    coord_equal(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
    labs(x = xlab, y = ylab) +
    theme_minimal() +
    theme(axis.title   = element_text(size = 15),
          axis.text    = element_text(size = 13),
          legend.title = element_text(size = 16),
          legend.text  = element_text(size = 14))
}
