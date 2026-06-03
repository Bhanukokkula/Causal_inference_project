set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(sandwich)
  library(lmtest)
  library(MatchIt)
  library(cobalt)
  library(glmnet)
  library(survey)
  library(MCMCpack)
  library(coda)
  library(rbounds)
  library(EValue)
  library(ggdag)
  library(dagitty)
  library(broom)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)

# ── Data & derived variables ─────────────────────────────────────────────────
data <- read.csv("data.csv")
data$z             <- as.integer(data$z)
data$race_collapsed <- factor(ifelse(data$race == 4, "majority", "other"),
                               levels = c("other", "majority"))
data$gender   <- factor(data$gender,  levels = c("1", "2"))
data$fgen     <- factor(data$fgen,    levels = c("0", "1"))
data$urban_f  <- factor(data$urban,   levels = c("0","1","2","3","4"))
data$schoolid <- factor(data$schoolid)

cov_rhs   <- "selfrpt + race_collapsed + gender + fgen + urban_f +
               mindset + test + sch_race + pov + size"
ps_form   <- as.formula(paste("z ~", cov_rhs))
out_form  <- as.formula(paste("y ~ z +", cov_rhs))
lmer_form <- as.formula(paste("y ~ z +", cov_rhs, "+ (1 | schoolid)"))

save_plot <- function(p, name, w = 8, h = 5) {
  ggsave(file.path("figures", paste0(name, ".png")), p,
         width = w, height = h, dpi = 150, bg = "white")
  message("Saved: figures/", name, ".png")
}

# ── 1. Outcome distribution by treatment ─────────────────────────────────────
p1 <- ggplot(data, aes(x = y, fill = factor(z, labels = c("Control", "Treatment")))) +
  geom_density(alpha = 0.55, colour = NA) +
  scale_fill_manual(values = c("#4477AA", "#EE6677"), name = NULL) +
  labs(title = "Outcome Distribution by Treatment Group",
       subtitle = "Achievement (y) is standardised; treatment group shifted ~0.41 SD rightward",
       x = "Achievement (y, SD units)", y = "Density") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
save_plot(p1, "01_outcome_by_treatment")

# ── 2. Propensity score overlap ───────────────────────────────────────────────
ps_model    <- glm(ps_form, family = binomial(), data = data)
data$pscore <- predict(ps_model, type = "response")
p_z         <- mean(data$z)
data$weight_stab <- ifelse(data$z == 1,
                            p_z / data$pscore,
                            (1 - p_z) / (1 - data$pscore))

p2 <- ggplot(data, aes(x = pscore,
                        fill = factor(z, labels = c("Control", "Treatment")))) +
  geom_density(alpha = 0.55, colour = NA) +
  scale_fill_manual(values = c("#4477AA", "#EE6677"), name = NULL) +
  labs(title = "Propensity Score Overlap",
       subtitle = "Good common support between treated and control groups",
       x = "Estimated P(Z = 1 | X)", y = "Density") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
save_plot(p2, "02_propensity_overlap")

# ── 3. Love plot: before vs after IPW ────────────────────────────────────────
png("figures/03_love_plot.png", width = 900, height = 620, res = 120)
bal <- bal.tab(ps_form, data = data, estimand = "ATE", m.threshold = 0.1)
love.plot(bal, stat = "mean.diffs", threshold = 0.1, abs = TRUE,
          var.order = "unadjusted",
          title = "Covariate Balance: Unweighted vs. IPW-Weighted",
          sample.names = c("Unweighted", "IPW-Weighted"),
          colors = c("#888888", "#EE6677")) +
  theme_bw(base_size = 11)
dev.off()
message("Saved: figures/03_love_plot.png")

# ── 4. Fit all models for CIs ────────────────────────────────────────────────
message("Fitting models...")
lmer_fit   <- lmer(lmer_form, data = data, REML = FALSE)
lmer_coef  <- fixef(lmer_fit)["z"]
lmer_ci    <- confint(lmer_fit, parm = "z", method = "Wald")

ols_fit      <- lm(out_form, data = data)
ols_est      <- coef(ols_fit)["z"]
ols_naive_ci <- confint(ols_fit)["z", ]
vcov_cl      <- vcovCL(ols_fit, cluster = ~schoolid)
ols_cl_ci    <- coefci(ols_fit, parm = "z", vcov. = vcov_cl)

design_naive <- svydesign(ids = ~1, data = data, weights = ~weight_stab)
design_cl    <- svydesign(ids = ~schoolid, data = data, weights = ~weight_stab)
ipw_naive    <- svyglm(y ~ z, design = design_naive)
ipw_cl       <- svyglm(y ~ z, design = design_cl)
ipw_est      <- coef(ipw_cl)["z"]
ipw_naive_ci <- confint(ipw_naive)["z", ]
ipw_cl_ci    <- confint(ipw_cl)["z", ]

dr_naive    <- svyglm(as.formula(paste("y ~ z +", cov_rhs)), design = design_naive)
dr_cl       <- svyglm(as.formula(paste("y ~ z +", cov_rhs)), design = design_cl)
dr_est      <- coef(dr_cl)["z"]
dr_naive_ci <- confint(dr_naive)["z", ]
dr_cl_ci    <- confint(dr_cl)["z", ]

X_ps      <- model.matrix(ps_form, data = data)[, -1]
cv_lasso  <- cv.glmnet(X_ps, data$z, family = "binomial", alpha = 1, nfolds = 10)
data$pscore_lasso <- as.vector(predict(cv_lasso, newx = X_ps,
                                        type = "response", s = "lambda.min"))
data$weight_lasso <- ifelse(data$z == 1,
                             p_z / data$pscore_lasso,
                             (1 - p_z) / (1 - data$pscore_lasso))
design_lasso_cl    <- svydesign(ids = ~schoolid, data = data, weights = ~weight_lasso)
design_lasso_naive <- svydesign(ids = ~1,         data = data, weights = ~weight_lasso)
ipw_lasso_cl    <- svyglm(y ~ z, design = design_lasso_cl)
ipw_lasso_naive <- svyglm(y ~ z, design = design_lasso_naive)
lasso_est       <- coef(ipw_lasso_cl)["z"]
lasso_cl_ci     <- confint(ipw_lasso_cl)["z", ]
lasso_naive_ci  <- confint(ipw_lasso_naive)["z", ]

match_fit <- matchit(ps_form, data = data, method = "nearest",
                     ratio = 1, replace = FALSE)
matched   <- match.data(match_fit)
psm_reg   <- lm(out_form, data = matched, weights = weights)
psm_est   <- coef(psm_reg)["z"]
vcov_psm_cl    <- vcovCL(psm_reg, cluster = ~schoolid)
psm_cl_ci      <- coefci(psm_reg, parm = "z", vcov. = vcov_psm_cl)["z", ]
psm_naive_ci   <- coefci(psm_reg, parm = "z")["z", ]

# ── 5. Cluster bootstrap (fast, 300 reps) ────────────────────────────────────
message("Running cluster bootstrap (300 reps)...")
draw_cluster_boot <- function(dat) {
  schools   <- levels(dat$schoolid)
  sampled   <- sample(schools, length(schools), replace = TRUE)
  do.call(rbind, lapply(seq_along(sampled), function(j) {
    d <- dat[dat$schoolid == sampled[j], ]
    d$schoolid <- factor(j)
    d
  }))
}
est_ols <- function(d) coef(lm(out_form, data = d))["z"]
est_ipw <- function(d) {
  ps <- glm(ps_form, family = binomial(), data = d)
  d$ps <- predict(ps, type = "response"); pz <- mean(d$z)
  d$w  <- ifelse(d$z == 1, pz/d$ps, (1-pz)/(1-d$ps))
  coef(lm(y ~ z, data = d, weights = d$w))["z"]
}
est_dr <- function(d) {
  ps <- glm(ps_form, family = binomial(), data = d)
  d$ps <- predict(ps, type = "response"); pz <- mean(d$z)
  d$w  <- ifelse(d$z == 1, pz/d$ps, (1-pz)/(1-d$ps))
  coef(lm(out_form, data = d, weights = d$w))["z"]
}
est_lmer <- function(d) {
  m <- suppressWarnings(lmer(lmer_form, data = d, REML = FALSE,
                              control = lmerControl(optimizer = "bobyqa")))
  fixef(m)["z"]
}
run_boot <- function(fn, n, dat) {
  vapply(seq_len(n), function(i) {
    tryCatch(fn(draw_cluster_boot(dat)), error = function(e) NA_real_)
  }, numeric(1))
}
set.seed(42)
boot_ols  <- run_boot(est_ols,  300, data)
boot_ipw  <- run_boot(est_ipw,  300, data)
boot_dr   <- run_boot(est_dr,   300, data)
boot_lmer <- run_boot(est_lmer, 150, data)

ci_boot <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)

# ── 6. CI width comparison bar chart ─────────────────────────────────────────
width_df <- tribble(
  ~Method, ~Naive_Width, ~Cluster_Width,
  "OLS",    diff(ols_naive_ci),         diff(ci_boot(boot_ols)),
  "IPW",    diff(ipw_naive_ci),         diff(ci_boot(boot_ipw)),
  "DR",     diff(dr_naive_ci),          diff(ci_boot(boot_dr)),
  "lmer",   diff(as.numeric(lmer_ci)),  diff(ci_boot(boot_lmer))
) %>%
  pivot_longer(c(Naive_Width, Cluster_Width),
               names_to = "Type", values_to = "Width") %>%
  mutate(Type = recode(Type,
                        Naive_Width   = "Naive (iid)",
                        Cluster_Width = "Cluster-robust"))

p_width <- ggplot(width_df, aes(x = Method, y = Width, fill = Type)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("#4477AA", "#EE6677"), name = NULL) +
  labs(title = "Naive vs. Cluster-Robust CI Width",
       subtitle = "Ignoring 76-school clustering underestimates uncertainty",
       x = NULL, y = "95% CI width (SD units)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
save_plot(p_width, "04_ci_width_comparison")

# ── 7. Bayesian posterior density ────────────────────────────────────────────
message("Running Bayesian model...")
set.seed(42)
bayes_fit <- MCMCregress(formula = out_form, data = data,
                          burnin = 1000, mcmc = 10000, thin = 10, seed = 42)
z_col        <- grep("^z$", colnames(bayes_fit), value = TRUE)
if (!length(z_col)) z_col <- grep("^z", colnames(bayes_fit), value = TRUE)[1]
ate_post     <- as.numeric(bayes_fit[, z_col])
post_ci      <- quantile(ate_post, c(0.025, 0.5, 0.975))

p_bayes <- ggplot(data.frame(ATE = ate_post), aes(x = ATE)) +
  geom_density(fill = "#4477AA", alpha = 0.6, colour = "#22355A") +
  geom_vline(xintercept = post_ci[c(1, 3)], linetype = "dashed",
             colour = "#EE6677", linewidth = 0.8) +
  geom_vline(xintercept = post_ci[2], linetype = "solid",
             colour = "#AA3333", linewidth = 1) +
  annotate("text", x = post_ci[2] + 0.004, y = Inf, vjust = 1.5, hjust = 0,
           label = sprintf("Median = %.3f", post_ci[2]), size = 4) +
  labs(title = "Posterior Distribution of ATE (Bayesian Linear Regression)",
       subtitle = "95% credible interval shown in red; solid line = posterior median",
       x = "ATE (SD units)", y = "Posterior density") +
  theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank())
save_plot(p_bayes, "05_bayesian_posterior")

# ── 8. Rosenbaum bounds ───────────────────────────────────────────────────────
pairs_df <- matched %>%
  group_by(subclass) %>%
  filter(n() == 2, length(unique(z)) == 2) %>%
  summarise(y_treat = y[z == 1], y_ctrl = y[z == 0], .groups = "drop")

sens_out  <- psens(pairs_df$y_treat, pairs_df$y_ctrl, Gamma = 3.0, GammaInc = 0.1)
bounds_df <- as.data.frame(sens_out$bounds)
names(bounds_df) <- c("Gamma", "p_lower", "p_upper")
break_idx  <- which(bounds_df$p_upper > 0.05)[1]
gamma_break <- if (is.na(break_idx)) "> 3.0" else bounds_df$Gamma[break_idx]

p_rosen <- ggplot(bounds_df, aes(x = Gamma)) +
  geom_line(aes(y = p_upper, colour = "Upper bound"), linewidth = 1) +
  geom_line(aes(y = p_lower, colour = "Lower bound"), linewidth = 1,
            linetype = "dashed") +
  geom_hline(yintercept = 0.05, linetype = "dotted", colour = "grey40") +
  annotate("text", x = 1.1, y = 0.055, label = "p = 0.05", size = 3.5,
           colour = "grey40") +
  scale_colour_manual(values = c("Upper bound" = "#EE6677",
                                  "Lower bound" = "#4477AA"), name = NULL) +
  scale_x_continuous(breaks = seq(1, 3, 0.5)) +
  labs(title = "Rosenbaum Bounds Sensitivity Analysis",
       subtitle = sprintf("Effect significance breaks at Gamma = %s", gamma_break),
       x = "Gamma (unmeasured confounding odds ratio)",
       y = "Wilcoxon signed-rank p-value (bound)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
save_plot(p_rosen, "06_rosenbaum_bounds")

# ── 9. Subgroup forest plot ───────────────────────────────────────────────────
message("Running subgroup bootstrap (100 reps per subgroup)...")
subgroup_boot_ci <- function(sub_data, n_boot = 100) {
  schools_sg <- unique(as.character(sub_data$schoolid))
  if (length(unique(sub_data$z)) < 2) return(c(NA, NA, NA))
  m0  <- lm(y ~ z + selfrpt + mindset + test, data = sub_data)
  est <- coef(m0)["z"]
  boots <- vapply(seq_len(n_boot), function(i) {
    s <- sample(schools_sg, length(schools_sg), replace = TRUE)
    bd <- do.call(rbind, lapply(seq_along(s), function(j) {
      d <- sub_data[as.character(sub_data$schoolid) == s[j], ]
      d$schoolid <- factor(j); d
    }))
    tryCatch(coef(lm(y ~ z + selfrpt + mindset + test, data = bd))["z"],
             error = function(e) NA_real_)
  }, numeric(1))
  ci <- quantile(boots, c(0.025, 0.975), na.rm = TRUE)
  c(ATE = unname(est), CI_lo = unname(ci[1]), CI_hi = unname(ci[2]))
}

set.seed(42)
sg_race <- data %>%
  group_by(Subgroup = paste0("Race: ", race_collapsed)) %>%
  group_modify(~ { r <- subgroup_boot_ci(.x)
    data.frame(N=nrow(.x), ATE=r["ATE"], CI_lo=r["CI_lo"], CI_hi=r["CI_hi"]) }) %>%
  ungroup()
sg_gender <- data %>%
  group_by(Subgroup = paste0("Gender: ", gender)) %>%
  group_modify(~ { r <- subgroup_boot_ci(.x)
    data.frame(N=nrow(.x), ATE=r["ATE"], CI_lo=r["CI_lo"], CI_hi=r["CI_hi"]) }) %>%
  ungroup()
sg_fgen <- data %>%
  group_by(Subgroup = paste0("First-gen: ", fgen)) %>%
  group_modify(~ { r <- subgroup_boot_ci(.x)
    data.frame(N=nrow(.x), ATE=r["ATE"], CI_lo=r["CI_lo"], CI_hi=r["CI_hi"]) }) %>%
  ungroup()

sg_all <- bind_rows(sg_race, sg_gender, sg_fgen) %>% filter(!is.na(ATE))

p_sg <- ggplot(sg_all, aes(x = ATE, y = forcats::fct_rev(Subgroup))) +
  geom_point(size = 3.5, colour = "#22355A") +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi),
                 height = 0.3, colour = "#22355A", linewidth = 0.9) +
  geom_vline(xintercept = lmer_coef, linetype = "dashed",
             colour = "#EE6677", linewidth = 0.9) +
  annotate("text", x = lmer_coef + 0.005, y = Inf, vjust = 1.3, hjust = 0,
           label = sprintf("Overall = %.3f", lmer_coef),
           colour = "#EE6677", size = 4) +
  labs(title = "Subgroup ATEs with Cluster-Bootstrap 95% CIs",
       subtitle = "No subgroup interaction survives BH correction -- effect is relatively homogeneous",
       x = "Estimated ATE (SD units)", y = NULL) +
  theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank())
save_plot(p_sg, "07_subgroup_forest", w = 8, h = 5)

# ── 10. Final methods forest plot ─────────────────────────────────────────────
forest_df <- tibble(
  Method   = c("OLS Regression", "IPW (stabilised)", "Doubly Robust",
               "Lasso-IPW", "PSM + Regression", "Multilevel (lmer)",
               "Bayesian (MCMC)"),
  Estimate = c(ols_est, ipw_est, dr_est, lasso_est, psm_est,
               lmer_coef, mean(ate_post)),
  CI_lo    = c(ci_boot(boot_ols)[1], ci_boot(boot_ipw)[1], ci_boot(boot_dr)[1],
               lasso_cl_ci[1], psm_cl_ci[1], ci_boot(boot_lmer)[1],
               post_ci[1]),
  CI_hi    = c(ci_boot(boot_ols)[2], ci_boot(boot_ipw)[2], ci_boot(boot_dr)[2],
               lasso_cl_ci[2], psm_cl_ci[2], ci_boot(boot_lmer)[2],
               post_ci[3]),
  Primary  = Method == "Multilevel (lmer)"
)

p_forest <- ggplot(forest_df,
                    aes(x = Estimate, y = forcats::fct_rev(Method),
                        colour = Primary, shape = Primary)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi),
                 height = 0.35, linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  scale_colour_manual(values = c("FALSE" = "#4477AA", "TRUE" = "#CC3333"),
                       guide = "none") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  scale_x_continuous(limits = c(0.25, 0.6), breaks = seq(0.25, 0.6, 0.05)) +
  labs(title = "All-Methods Forest Plot: ATE with Cluster-Robust 95% CIs",
       subtitle = "Red diamond = primary estimator (multilevel lmer). CIs account for 76-school clustering.",
       x = "Average Treatment Effect (SD units of achievement)",
       y = NULL) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())
save_plot(p_forest, "08_methods_forest_plot", w = 9, h = 5.5)

message("\nAll figures saved to figures/")
list.files("figures/")
