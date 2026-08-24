#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

set.seed(25025)

output_dir <- file.path("images", "generated")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set FIGURE_QA_DIR to emit full-resolution PNG previews without changing the
# committed SVG-only book output. This is useful for visual regression checks.
qa_dir <- Sys.getenv("FIGURE_QA_DIR", unset = "")
if (nzchar(qa_dir)) {
  dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
}

palette <- list(
  blue = "#2364AA",
  orange = "#E6862A",
  green = "#278A64",
  gray = "#68717C",
  light_gray = "#D8DEE6",
  purple = "#7656A5",
  ink = "#18212B"
)

# A finite-dimensional check of the Hilbert-space projection identity used in
# Chapters 8 and 25: the residual must be orthogonal to the nuisance span.
nuisance_basis <- cbind(c(1, 0, 1, -1), c(0, 1, 1, 1))
ordinary_score <- c(2, -1, 3, 1)
nuisance_projection <- nuisance_basis %*%
  solve(crossprod(nuisance_basis), crossprod(nuisance_basis, ordinary_score))
efficient_residual <- ordinary_score - nuisance_projection
stopifnot(max(abs(crossprod(nuisance_basis, efficient_residual))) < 1e-12)

book_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", colour = palette$ink),
      plot.subtitle = element_text(colour = palette$gray),
      plot.caption = element_text(colour = palette$gray, hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E8ECF1", linewidth = 0.35),
      strip.text = element_text(face = "bold", colour = palette$ink),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.title = element_text(colour = palette$ink),
      axis.text = element_text(colour = palette$ink)
    )
}

save_svg <- function(plot, filename, width = 9, height = 5.4) {
  svglite(file.path(output_dir, filename), width = width, height = height,
          bg = "white", system_fonts = list(sans = "Arial", serif = "Times New Roman"))
  print(plot)
  invisible(dev.off())

  if (nzchar(qa_dir)) {
    ggsave(
      file.path(qa_dir, sub("\\.svg$", ".png", filename)),
      plot = plot, device = ragg::agg_png, width = width, height = height,
      units = "in", dpi = 160, bg = "white"
    )
  }
}

# Chapter 6: exponential tilting of a joint Gaussian --------------------------
c_tilt <- 0.65
sigma <- matrix(c(1, c_tilt, c_tilt, 1), 2, 2)
grid <- expand.grid(
  statistic = seq(-3.2, 3.8, length.out = 180),
  log_lr = seq(-3.2, 4.0, length.out = 180)
)

bvn_density <- function(x, mean, sigma) {
  coordinates <- as.matrix(x[, c("statistic", "log_lr")])
  centered <- sweep(coordinates, 2, mean)
  const <- 1 / (2 * pi * sqrt(det(sigma)))
  const * exp(-0.5 * rowSums((centered %*% solve(sigma)) * centered))
}

grid$reference <- bvn_density(grid, c(0, 0), sigma)
grid$tilted <- bvn_density(grid, c(c_tilt, 1), sigma)

dx <- diff(unique(grid$statistic))[1]
dy <- diff(unique(grid$log_lr))[1]
stopifnot(abs(sum(grid$reference) * dx * dy - 1) < 0.02)
stopifnot(abs(sum(grid$tilted) * dx * dy - 1) < 0.02)
stopifnot(abs(c_tilt - sigma[1, 2]) < 1e-12)
tilted_from_reference <- grid$reference * exp(grid$log_lr - 0.5)
tilt_normalization <- sum(tilted_from_reference) * dx * dy
tilt_mean <- sum(grid$statistic * tilted_from_reference) * dx * dy /
  tilt_normalization
stopifnot(abs(tilt_normalization - 1) < 0.02)
stopifnot(abs(tilt_mean - c_tilt) < 0.02)

p_tilt_joint <- ggplot(grid, aes(statistic, log_lr)) +
  geom_contour(aes(z = reference, colour = "Reference law"), bins = 7, linewidth = 0.75) +
  geom_contour(aes(z = tilted, colour = "After likelihood tilt"), bins = 7,
               linewidth = 0.75, linetype = 2) +
  geom_point(data = data.frame(statistic = c(0, c_tilt), log_lr = c(0, 1)),
             aes(statistic, log_lr, colour = c("Reference law", "After likelihood tilt")),
             size = 2.5, inherit.aes = FALSE) +
  scale_colour_manual(values = c("Reference law" = palette$blue,
                                 "After likelihood tilt" = palette$orange)) +
  coord_equal() +
  labs(title = "A likelihood tilt moves the joint Gaussian mean",
       subtitle = "The covariance ellipse keeps its shape; the statistic shifts by Cov(T, Z)",
       x = "Statistic T", y = "Centered likelihood coordinate Z") +
  book_theme()

marginal <- rbind(
  data.frame(t = seq(-3.2, 3.8, length.out = 400), law = "Reference law"),
  data.frame(t = seq(-3.2, 3.8, length.out = 400), law = "After likelihood tilt")
)
marginal$density <- ifelse(
  marginal$law == "Reference law",
  dnorm(marginal$t),
  dnorm(marginal$t, mean = c_tilt)
)

p_tilt_marginal <- ggplot(marginal, aes(t, density, colour = law, linetype = law)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(0, c_tilt), colour = c(palette$blue, palette$orange),
             linewidth = 0.55) +
  scale_colour_manual(values = c("Reference law" = palette$blue,
                                 "After likelihood tilt" = palette$orange)) +
  scale_linetype_manual(values = c("Reference law" = 1, "After likelihood tilt" = 2)) +
  labs(title = "Marginal consequence", subtitle = "T: N(0,1) becomes N(c,1)",
       x = "T", y = "Density") +
  book_theme()

save_svg(p_tilt_joint + p_tilt_marginal + plot_layout(widths = c(1.35, 1)),
         "chapter-6-exponential-tilt.svg", width = 11, height = 5.2)

# Chapter 7: exact Poisson likelihood and its LAN approximation ---------------
lambda0 <- 2
h_grid <- seq(-2.6, 3, length.out = 300)
n_values <- c(12, 60, 300)
poisson_panels <- lapply(n_values, function(n) {
  sum_x <- rpois(1, n * lambda0)
  lambda_h <- lambda0 + h_grid / sqrt(n)
  exact <- sum_x * log(lambda_h / lambda0) - n * (lambda_h - lambda0)
  delta_n <- (sum_x - n * lambda0) / (lambda0 * sqrt(n))
  lan <- h_grid * delta_n - 0.5 * h_grid^2 / lambda0
  rbind(
    data.frame(h = h_grid, value = exact, curve = "Exact local log likelihood", n = n),
    data.frame(h = h_grid, value = lan, curve = "LAN quadratic", n = n)
  )
})
poisson_data <- do.call(rbind, poisson_panels)

p_lan <- ggplot(poisson_data, aes(h, value, colour = curve, linetype = curve)) +
  geom_hline(yintercept = 0, colour = palette$light_gray, linewidth = 0.5) +
  geom_line(linewidth = 0.95) +
  facet_wrap(~ n, nrow = 1, labeller = label_bquote(n == .(n))) +
  scale_colour_manual(values = c("Exact local log likelihood" = palette$blue,
                                 "LAN quadratic" = palette$orange)) +
  scale_linetype_manual(values = c("Exact local log likelihood" = 1,
                                   "LAN quadratic" = 2)) +
  labs(title = "The local Poisson likelihood becomes quadratic",
       subtitle = expression(lambda == lambda[0] + h/sqrt(n) ~~~ (lambda[0] == 2)),
       x = "Local parameter h", y = "Log likelihood ratio") +
  book_theme()

save_svg(p_lan, "chapter-7-lan-quadratic.svg", width = 11, height = 4.8)

# Chapter 8: concentration, Hodges risk, and bowl-shaped loss -----------------
x <- seq(-4, 4, length.out = 600)
normal_shift <- rbind(
  data.frame(x = x, density = dnorm(x), distribution = "N(0, 1)"),
  data.frame(x = x, density = dnorm(x, 1.1), distribution = "N(1.1, 1)")
)
normal_scale <- rbind(
  data.frame(x = x, density = dnorm(x), distribution = "N(0, 1)"),
  data.frame(x = x, density = dnorm(x, sd = 0.55), distribution = "N(0, 0.55^2)")
)

p_shift <- ggplot(normal_shift, aes(x, density, colour = distribution, linetype = distribution)) +
  annotate("rect", xmin = -1, xmax = 1, ymin = 0, ymax = Inf,
           fill = palette$blue, alpha = 0.08) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = 3, colour = palette$gray) +
  scale_colour_manual(values = c("N(0, 1)" = palette$blue, "N(1.1, 1)" = palette$orange)) +
  scale_linetype_manual(values = c("N(0, 1)" = 1, "N(1.1, 1)" = 2)) +
  labs(title = "Fixed variance", subtitle = "Centering puts more mass near zero",
       x = "Error", y = "Density") + book_theme()

p_scale <- ggplot(normal_scale, aes(x, density, colour = distribution, linetype = distribution)) +
  annotate("rect", xmin = -1, xmax = 1, ymin = 0, ymax = Inf,
           fill = palette$green, alpha = 0.08) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = 3, colour = palette$gray) +
  scale_colour_manual(values = c("N(0, 0.55^2)" = palette$green, "N(0, 1)" = palette$blue)) +
  scale_linetype_manual(values = c("N(0, 0.55^2)" = 1, "N(0, 1)" = 2)) +
  labs(title = "Fixed center", subtitle = "Smaller variance concentrates mass",
       x = "Error", y = "Density") + book_theme()

save_svg(p_shift + p_scale, "chapter-8-normal-concentration.svg", width = 11, height = 4.8)

hodges_risk <- function(theta, n) {
  threshold <- n^(-1 / 4)
  lower <- sqrt(n) * (-threshold - theta)
  upper <- sqrt(n) * (threshold - theta)
  middle_prob <- pnorm(upper) - pnorm(lower)
  lower_second <- pnorm(lower) - lower * dnorm(lower)
  upper_second <- 1 - pnorm(upper) + upper * dnorm(upper)
  n * theta^2 * middle_prob + lower_second + upper_second
}

theta_grid <- seq(-1.25, 1.25, length.out = 900)
hodges_data <- do.call(rbind, lapply(c(10, 100, 1000), function(n) {
  data.frame(theta = theta_grid, risk = hodges_risk(theta_grid, n), n = factor(n))
}))

p_hodges <- ggplot(hodges_data, aes(theta, risk, colour = n, linetype = n)) +
  geom_hline(yintercept = 1, linewidth = 0.55, colour = palette$gray) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("10" = palette$gray, "100" = palette$orange,
                                 "1000" = palette$blue)) +
  scale_linetype_manual(values = c("10" = 3, "100" = 2, "1000" = 1)) +
  coord_cartesian(ylim = c(0, 20)) +
  labs(title = "Hodges' pointwise gain hides moving risk peaks",
       subtitle = bquote(n * E[theta] * (S[n] - theta)^2 ~
                           "; horizontal line: sample-mean risk"),
       x = expression(theta), y = "Scaled quadratic risk", colour = "Sample size",
       linetype = "Sample size") +
  book_theme()

save_svg(p_hodges, "chapter-8-hodges-risk.svg", width = 9.5, height = 5.2)

circle <- data.frame(angle = seq(0, 2 * pi, length.out = 300))
circle$x <- cos(circle$angle)
circle$y <- sin(circle$angle)
diamond <- data.frame(x = c(0, 1, 0, -1, 0), y = c(1, 0, -1, 0, 1))

p_sublevel <- ggplot() +
  geom_polygon(data = circle, aes(x, y), fill = palette$blue, alpha = 0.12,
               colour = palette$blue, linewidth = 0.9) +
  geom_polygon(data = diamond, aes(x, y), fill = palette$green, alpha = 0.12,
               colour = palette$green, linewidth = 0.9) +
  coord_equal() +
  annotate("text", x = 0.25, y = 0.72, label = "squared / L2", colour = palette$blue,
           hjust = 0, size = 3.7) +
  annotate("text", x = 0.18, y = -0.72, label = "L1", colour = palette$green,
           hjust = 0, size = 3.7) +
  labs(title = "Convex, symmetric sublevel sets", x = expression(x[1]), y = expression(x[2])) +
  book_theme() + theme(panel.grid = element_blank())

loss_x <- seq(-2, 2, length.out = 400)
loss_data <- rbind(
  data.frame(x = loss_x, loss = loss_x^2, kind = "Squared loss"),
  data.frame(x = loss_x, loss = abs(loss_x), kind = "Absolute loss")
)
p_losses <- ggplot(loss_data, aes(x, loss, colour = kind, linetype = kind)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("Squared loss" = palette$blue,
                                 "Absolute loss" = palette$green)) +
  scale_linetype_manual(values = c("Squared loss" = 1, "Absolute loss" = 2)) +
  labs(title = "Bowl-shaped losses", subtitle = "Increasing away from zero",
       x = "Estimation error", y = "Loss") + book_theme()

save_svg(p_sublevel + p_losses, "chapter-8-bowl-shaped-loss.svg", width = 10.5, height = 4.6)

# Chapter 15: local asymptotic power ------------------------------------------
alpha <- 0.05
delta <- seq(0, 4, length.out = 350)
power_data <- do.call(rbind, lapply(c(1, 0.75, 0.4), function(rho) {
  data.frame(delta = delta,
             power = 1 - pnorm(qnorm(1 - alpha) - rho * delta),
             alignment = factor(rho, levels = c(1, 0.75, 0.4)))
}))
stopifnot(all(diff(subset(power_data, alignment == 1)$power) >= 0))
stopifnot(all(subset(power_data, alignment == 1)$power + 1e-12 >=
              subset(power_data, alignment == 0.75)$power))
stopifnot(all(subset(power_data, alignment == 0.75)$power + 1e-12 >=
              subset(power_data, alignment == 0.4)$power))

p_power <- ggplot(power_data, aes(delta, power, colour = alignment, linetype = alignment)) +
  geom_hline(yintercept = alpha, colour = palette$gray, linewidth = 0.55, linetype = 3) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("1" = palette$green, "0.75" = palette$blue,
                                 "0.4" = palette$orange),
                      labels = c("Efficient direction (rho = 1)", "Partial alignment (rho = 0.75)",
                                 "Weak alignment (rho = 0.4)")) +
  scale_linetype_manual(values = c("1" = 1, "0.75" = 2, "0.4" = 3),
                        labels = c("Efficient direction (rho = 1)", "Partial alignment (rho = 0.75)",
                                   "Weak alignment (rho = 0.4)")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Local power is controlled by alignment",
       subtitle = "Equal-norm local displacements; one-sided test at alpha = 0.05",
       x = "Local-displacement norm delta", y = "Power") + book_theme()

save_svg(p_power, "chapter-15-local-power.svg", width = 9.5, height = 5.2)

# Chapter 18: pointwise versus uniform control --------------------------------
x01 <- seq(0, 1, length.out = 1200)
moving <- do.call(rbind, lapply(c(5, 10, 20), function(n) {
  data.frame(x = x01, value = exp(-((x01 - 1 / n) / (0.22 / n^2))^2),
             n = factor(n), panel = "Pointwise, but not uniform")
}))
shrinking <- do.call(rbind, lapply(c(5, 10, 20), function(n) {
  data.frame(x = x01, value = sin(2 * pi * x01) / sqrt(n),
             n = factor(n), panel = "Uniformly shrinking")
}))
function_data <- rbind(moving, shrinking)

p_uniform <- ggplot(function_data, aes(x, value, colour = n, linetype = n)) +
  geom_hline(yintercept = 0, colour = palette$light_gray) +
  geom_line(linewidth = 0.85) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("5" = palette$orange, "10" = palette$blue,
                                 "20" = palette$green)) +
  scale_linetype_manual(values = c("5" = 3, "10" = 2, "20" = 1)) +
  labs(title = "Pointwise convergence does not control an entire function",
       subtitle = "The moving spike has sup norm 1 for every n; the right-hand sequence shrinks uniformly",
       x = "Index t", y = expression(f[n](t)), colour = "n", linetype = "n") +
  book_theme()

save_svg(p_uniform, "chapter-18-pointwise-uniform.svg", width = 11, height = 4.8)

# Chapter 19: empirical process and Z-estimation ------------------------------
population <- faithful$eruptions
set.seed(19019)
sample_x <- sample(population, 65, replace = TRUE)
grid_x <- seq(min(population) - 0.15, max(population) + 0.15, length.out = 450)
F_population <- ecdf(population)(grid_x)
F_sample <- ecdf(sample_x)(grid_x)
ecdf_data <- data.frame(
  x = rep(grid_x, 2),
  value = c(F_population, F_sample),
  curve = rep(c("Reference empirical distribution", "Sample empirical distribution"),
              each = length(grid_x))
)
process_data <- data.frame(x = grid_x,
                           value = sqrt(length(sample_x)) * (F_sample - F_population))

p_ecdf <- ggplot(ecdf_data, aes(x, value, colour = curve, linetype = curve)) +
  geom_step(linewidth = 0.85) +
  scale_colour_manual(values = c("Reference empirical distribution" = palette$gray,
                                 "Sample empirical distribution" = palette$blue)) +
  scale_linetype_manual(values = c("Reference empirical distribution" = 2,
                                   "Sample empirical distribution" = 1)) +
  labs(title = "Empirical distributions are random functions",
       x = "Old Faithful eruption duration (minutes)", y = "Cumulative probability") +
  book_theme(base_size = 11)

p_process <- ggplot(process_data, aes(x, value)) +
  geom_hline(yintercept = 0, colour = palette$gray, linewidth = 0.5) +
  geom_step(colour = palette$orange, linewidth = 0.85) +
  labs(title = "Centered empirical process", x = "Eruption duration (minutes)",
       y = "Centered process") + book_theme(base_size = 11)

theta <- seq(1.6, 4.6, length.out = 300)
z_data <- rbind(
  data.frame(theta = theta, psi = mean(population) - theta,
             equation = "Reference equation"),
  data.frame(theta = theta, psi = mean(sample_x) - theta,
             equation = "Sample equation")
)
p_z <- ggplot(z_data, aes(theta, psi, colour = equation, linetype = equation)) +
  geom_hline(yintercept = 0, colour = palette$gray, linewidth = 0.5) +
  geom_line(linewidth = 0.95) +
  geom_vline(xintercept = c(mean(population), mean(sample_x)),
             colour = c(palette$gray, palette$green), linetype = c(2, 1)) +
  scale_colour_manual(values = c("Reference equation" = palette$gray,
                                 "Sample equation" = palette$green)) +
  scale_linetype_manual(values = c("Reference equation" = 2, "Sample equation" = 1)) +
  labs(title = "Z-estimation moves the zero crossing",
       x = expression(theta), y = expression(Psi(theta))) + book_theme(base_size = 11)

save_svg(p_ecdf | p_process | p_z,
         "chapter-19-empirical-z.svg", width = 13, height = 4.6)

# Chapter 25: a partially linear model ----------------------------------------
set.seed(2511)
n_pl <- 700
w <- runif(n_pl, -2, 2)
m_w <- 0.8 * sin(w) + 0.2 * w
v <- rnorm(n_pl, sd = 0.75)
d <- m_w + v
theta0 <- 1.2
g_w <- cos(1.5 * w) + 0.3 * w^2
eps <- rnorm(n_pl, sd = 0.8)
y <- theta0 * d + g_w + eps
ell_w <- theta0 * m_w + g_w
res_d <- d - m_w
res_y <- y - ell_w
orthogonal_moment <- mean(res_d * (res_y - theta0 * res_d))
stopifnot(abs(orthogonal_moment) < 0.08)

pl_data <- data.frame(w = w, d = d, y = y, res_d = res_d, res_y = res_y)

p_raw <- ggplot(pl_data, aes(d, y, colour = w)) +
  geom_point(alpha = 0.45, size = 1.25) +
  geom_smooth(method = "lm", se = FALSE, colour = palette$blue, linewidth = 0.9) +
  scale_colour_gradient2(low = palette$orange, mid = palette$gray, high = palette$purple,
                         midpoint = 0) +
  labs(title = "Before removing nuisance structure",
       subtitle = "Both treatment and outcome depend nonlinearly on W",
       x = "D", y = "Y", colour = "W") + book_theme()

p_res <- ggplot(pl_data, aes(res_d, res_y)) +
  geom_point(alpha = 0.45, size = 1.25, colour = palette$green) +
  geom_abline(slope = theta0, intercept = 0, colour = palette$blue, linewidth = 1) +
  labs(title = "After residualization",
       subtitle = "The remaining slope is the target after orthogonalization",
       x = expression(D - E[D * "|" * W]), y = expression(Y - E[Y * "|" * W])) +
  book_theme()

save_svg(p_raw + p_res, "chapter-25-partially-linear.svg", width = 11, height = 5)

message("Rendered reproducible SVG figures to ", normalizePath(output_dir))
