#' R script to reproduce the results in:
#'   Crippa A, Discacciati A, Bottai M, Spiegelman D, Orsini N.
#'   One-stage dose-response meta-analysis for aggregated data.
#'   BMC medical research methodology. 2016 Aug 2;16(1):91.

# Packages required ----
library(tidyverse)
# devtools::install_github("alecri/dosresmeta")
library(dosresmeta)
library(rms)
library(lmtest)
library(gridExtra)
library(readxl)
library(knitr)
theme_set(theme_minimal())

# 4 Practical examples ----  

# 4.1 Example 1: fictitious data example ----

# true beta coefficients generating the simulated data 'sim_os'
bi <- structure(c(-3.02773781522055, -2.71360284216677, -2.81718024023898, 
                  -2.98525410370846, -3.07807174460969, -3.08090019565107, -2.76123148090383, 
                  -2.90138487974716, -3.23406524104945, -0.35, -0.246168150617906, 
                  -0.230850306644441, -0.35, -0.274625537367731, -0.267349209860969, 
                  -0.35, -0.233102337406713, -0.252156341820433, 0.025, 0.0259126663907604, 
                  0.0244049682211859, 0.025, 0.032217031434239, 0.0343394212681985, 
                  0.025, 0.0300706894228835, 0.0300235286430353), dim = c(9L, 3L),
                dimnames = list(1:9, paste("beta", 0:2)))
apply(bi, 2, mean)[-1]
apply(bi, 2, sd)[-1]
kable(bi)

# load simulated data (Table 2)
data("sim_os")
kable(sim_os, digits = 2)

# number of non-referent log rr
group_by(sim_os, id) %>% 
  summarise(n_nonref = sum(se != 0))

# estimation: one stage and two-stage drma (on a subset)
mod_os <- dosresmeta(logrr ~ dose + I(dose^2), id = id, cases = cases, n = n, type = type,
                     se = se, data = sim_os, proc = "1stage")
summary(mod_os)

id_ts <- names(which(table(sim_os$id) > 2))
mod_ts <- dosresmeta(logrr ~ dose + I(dose^2), id = id, cases = cases, n = n, type = type,
                     se = se, data = subset(sim_os, id %in% id_ts), proc = "2stage")
summary(mod_ts)

# comparison coefficients
data.frame(one_stage = mod_os$coefficients, two_stage = c(mod_ts$coefficients)) %>% 
  kable(digits = 3)

# comparison standard errors
data.frame(one_stage = mod_os$vcov[c(1, 4)]^.5, two_stage = mod_ts$vcov[c(1, 4)]^.5) %>% 
  kable(digits = 3)

# comparison prediction
data.frame(dose = c(1.5, 1, 5, 9)) %>%
  bind_cols(predict(mod_os, ., expo = T)[, 3:5],
            predict(mod_ts, ., expo = T)[, 3:5]) %>%
  mutate(diff_os = -ci.lb + ci.ub,
         diff_ts = -ci.lb1 + ci.ub1) %>% 
  kable(digits = 2)

# comparison variance components
data.frame(one_stage = mod_os$Psi[-2], two_stage = mod_ts$Psi[-2])

# graphical comparison using 1.5 as referent (Figure 2)
xref <- 1.5
newd <- data.frame(dose = c(xref, seq(0, max(sim_os$dose), length.out = 50)))
pred_os <- predict(mod_os, newd, xref = xref, expo = T)
pred_ts <- predict(mod_ts, newd, xref = xref, expo = T)
pred_true <- tibble(
  dose = newd$dose,
  pred = c(exp(cbind(dose - xref, dose^2- xref^2) %*% c(-0.27, 0.027)))
  )

ggplot(pred_ts, aes(dose, y = pred)) + 
  geom_line(aes(linetype = "Two-stage")) + 
  geom_line(data = pred_os, aes(linetype = "One-stage")) + 
  geom_line(data = pred_true, aes(linetype = "True")) + 
  scale_y_continuous(trans = "log", breaks = c(.7, 1, 1.7, 2.5)) +
  labs(y = "Odds ratio", x = "Dose", linetype = "Curve") +
  scale_linetype_manual(values = c(`Two-stage` = "dotted", `One-stage` = "dashed", True = "solid"))

# individual curves (Table 3)
bi_blup <- t(apply(blup(mod_os), 1, function(x) x + rbind(coef(mod_os))))
bi_ts <- t(apply(blup(mod_ts), 1, function(x) x + rbind(coef(mod_ts)))) %>% 
  as.data.frame() %>% cbind(id = id_ts)
bi_comp <- data.frame(id = seq_along(unique(sim_os$id)),
                      xi_ref = sim_os$dose[sim_os$se == 0],
                      true = bi[, -1], blup = bi_blup) %>%
  merge(bi_ts, all = TRUE)
colnames(bi_comp)[-2] <- c("ID", paste0(rep(c("true", "os", "ts"), each = 2), c(".b1", ".b2")))
kable(bi_comp)

# auxiliary function for individual predictions
pred_sq <- function(coef, x, xref){
  y <- do.call("cbind", Map(function(b, r){
    exp(cbind(x - r, x^2- r^2) %*% t(b))
  }, split(coef, 1:nrow(coef)), split(xref, seq_along(xref))))
  colnames(y) <- paste0("pred", 1:ncol(y))
  data.frame(y)
}
predi <- data.frame(
  x = newd$dose,
  y_true = pred_sq(coef =  bi_comp[, c("true.b1", "true.b2")], newd$dose, bi_comp$xi_ref),
  y_os = pred_sq(coef =  bi_comp[, c("os.b1", "os.b2")], newd$dose, bi_comp$xi_ref),
  y_ts = pred_sq(coef =  bi_comp[, c("ts.b1", "ts.b2")], newd$dose, bi_comp$xi_ref)
) %>% 
  pivot_longer(
    cols = c("y_true.pred1":"y_ts.pred9"),
    names_to = c("pred_type", "id"),
    names_pattern = "y_(.*).pred(.*)"
  ) %>% 
  mutate(
    pred_type = factor(pred_type, levels = c("true", "os", "ts"), 
                       labels = c("True", "One-stage", "Two-stage")),
    study_lab = paste("Study ID", id)
  ) %>% 
  merge(sim_os, by = "id")

# Figure 3
ggplot(predi, aes(x = dose, y = rr)) + 
  geom_errorbar(data = select(predi, dose, rr, lrr, urr, study_lab) %>% distinct(),
                aes(ymin = lrr, ymax = urr), width = .3) +
  geom_line(aes(x = x, y = value, linetype = pred_type)) +
  scale_y_continuous(trans = "log", breaks = c(.3, .6, 1, 2, 4)) +
  facet_wrap(. ~ study_lab, scale = "free") +
  scale_linetype_manual(values =  c(`True` = "solid", `One-stage` = "dashed", `Two-stage` = "dotted")) +
  labs(x = "Dose", y = "Odds Ratio", linetype = "Curve") + 
  theme_classic()

# vpc plot(Figure 4)
p_vpc <- ggplot(subset(sim_os, se != 0), aes(dose, vpc(mod_os))) +
  geom_point() + geom_smooth(method = "loess", se = F, col = "black") +
  labs(y = "VPC", x = "Dose") + theme_classic()
p_vpc
# compare with a singular heterogeneity measure (I^2)
summary(mod_ts)

# Model comparison and goodness-of-fit
# more complex models: cubic, spline, fracpol
mod_cub <- dosresmeta(logrr ~ dose + I(dose^2) + I(dose^3), id = id, 
                      cases = cases, n = n, type = type, se = se, data = sim_os, 
                      proc = "1stage", method = "ml", 
                      control = list(maxiter = 1500))
k <- quantile(sim_os$dose, c(.05, .2, .5, .8, .95))
mod_spl <- dosresmeta(logrr ~ rcs(dose, k), id = id, cases = cases, n = n, 
                      type = type, se = se, data = sim_os, proc = "1stage",
                      method = "ml", control = list(maxiter = 5000))
p <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
grid <- subset(expand.grid(p1 = p, p2 = p), p1 <= p2)
rownames(grid) <- seq(nrow(grid))
shift <- 0
scale <- 1
modi_pi <- lapply(split(grid, seq(nrow(grid))), function(p)
  dosresmeta(logrr ~ fracpol(dose, p = p, shift = shift, scale = scale), id = id,
             type = type, se = se, cases = cases, n = n, data = sim_os, 
             proc = "1stage", method = "ml"))
grid[which.min(lapply(modi_pi, AIC)), ]
mod_fracpol <- modi_pi[[which.min(lapply(modi_pi, AIC))]]

pred_cub <- predict(mod_cub, newd, xref = xref, expo = T)
pred_spl <- predict(mod_spl, newd, xref = xref, expo = T)
pred_fracpol <- predict(mod_fracpol, subset(newd, dose >0), expo = T)

# Figure 7
ggplot(pred_cub, aes(newd$dose, pred, linetype = "Cubic")) + geom_line() +
  geom_line(data = pred_spl, aes(linetype = "Restricted cubic splines")) +
  geom_line(data = pred_fracpol, aes(x = newd$dose[newd$dose>0], linetype = "Fractional polynomials (0.5, 3)")) +
  scale_y_continuous(trans = "log") + labs(x = "Dose", y = "Odds ratio", linetype = "Curve") +
  scale_linetype_manual(values = c(Cubic = "dotted", `Restricted cubic splines` = "dashed", `Fractional polynomials (0.5, 3)` = "solid")) +
  theme_classic()

# comparison AICs
data.frame("cubic" = AIC(mod_cub), "spline" = AIC(mod_spl), 
           "frac pol" = AIC(mod_fracpol)) %>% 
  kable(digits = 3)

# gof plot
sq_gof <- gof(mod_os)
sq_gof

# Figure 8
ggplot(subset(sim_os, se != 0), aes(dose, sq_gof$tdata$tresiduals)) + 
  geom_point() + geom_smooth(method = "loess", se = F, col = "black") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_y_continuous(lim = c(-4, 4)) +
  labs(y = "Decorrelated residuals", x = "Dose") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))


# 4.2 Example 2: coffee consumption and all-cause mortality ----

# load data
# with the new version of dosresmeta just type `data("coffee_mort_add")`

load(url("https://github.com/alecri/one-stage-dosresmeta/raw/master/coffee_mort_add.rda"))

# reconstructing covariances
k <- c(0.5, 2.0, 4.5) 
Slist <- lapply(unique(coffee_mort_add$id), function(i)
  with(subset(coffee_mort_add, id == i), {
    if (any(is.na(cases) | is.na(n))){
      diag(se[se != 0 & !is.na(se)]^2, nrow = sum(se != 0 & !is.na(se)))
    }
    else {
      covar.logrr(y = logrr, v = I(se^2), cases = cases, n = n,
                  type = type)
    }
  }))
names(Slist) <- unique(coffee_mort_add$id)

# data for plot and prediction
newd <- data.frame(dose = seq(0, 8, length.out = 100))
newd_tab <- data.frame(dose = c(0, 2, 4))

# Splines without exclusion
spl <- dosresmeta(logrr ~ rcs(dose, k), id = id, type = type,
                  cases = cases, n = n, se = se, data = coffee_mort_add, 
                  covariance = "user", Slist = Slist, proc = "1stage", method = "reml")
summary(spl)
# Splines with exclusion
spl_exc <- dosresmeta(logrr ~ rcs(dose, k), id = id, type = type, cases = cases, 
                      n = n, se = se, data = subset(coffee_mort_add, !(id %in% c(28, 29))), method = "reml",
                      covariance = "user", Slist = Slist[!names(Slist) %in% c("28", "29")])
summary(spl_exc)

# comparison predicted RRS
round(predict(spl, newd_tab, expo = T), 2)
round(predict(spl_exc, newd_tab, expo = T), 2)

# comparison coefficients
round(rbind(spl = c(coef(spl), vcov(spl)[-2]),
            spl_exc = c(coef(spl_exc), vcov(spl_exc)[-2])), 5)

# Figure 5
pred_spl <- predict(spl, newd, expo = T)
pred_spl_excl <- predict(spl_exc, newd, expo = T)
ggplot(pred_spl, aes(newd$dose, pred, linetype = "One-stage")) + 
  geom_line() +
  geom_line(data = pred_spl_excl, aes(linetype = "Two-stage")) +
  scale_y_continuous(trans = "log", breaks = c(1, .75, .8, .85, .9, .95)) + 
  labs(x = "Coffee consumption (cups/day)", y = "Relative risk", linetype = "Curve") +
  scale_linetype_manual(values = c(`One-stage` = "solid", `Two-stage` = "dashed")) +
  theme_classic()


# alternative models (ML estimation for likelihood comparison)
newd <- data.frame(dose = c(1, seq(0, 8, length.out = 100)))
spl_ml <- dosresmeta(logrr ~ rcs(dose, k), id = id, type = type, cases = cases, n = n, 
                     se = se, data = coffee_mort_add, covariance = "user",  Slist = Slist,
                     proc = "1stage", method = "ml")
pred_spl_ml <- data.frame(dose = newd$dose, pred_spl = predict(spl_ml, newd, expo = T)$pred)

# spike at zero
spl_spike <- dosresmeta(logrr ~ I(1*(dose < 1)) + I(rcs(dose, k)*(dose >= 1)), 
                        id = id, type = type,
                        cases = cases, n = n, se = se, data = coffee_mort_add, 
                        covariance = "user", proc = "1stage", method = "ml",
                        Slist = Slist, control = list(maxiter = 5000))
summary(spl_spike)
pred_spl_spike <- predict(spl_spike, newd, expo = T)

# categorical
k2 <- c(0, 1, 3, 5, 7, 10)
categ <- dosresmeta(logrr ~ relevel(cut(dose, breaks = k2, include.lowest = T, right = F), 2), 
                    id = id, type = type,
                    cases = cases, n = n, se = se, data = coffee_mort_add, 
                    covariance = "user", proc = "1stage", method = "reml",
                    Slist = Slist, control = list(maxiter = 5000))
summary(categ)
pred_categ <- predict(categ, newd, expo = T)

# Figure 9
pred_modi <- cbind(pred_spl_ml, 
                   pred_spl_spike = pred_spl_spike$pred,
                   pred_categ = pred_categ$pred)

p <- ggplot(pred_modi, aes(dose, pred_spl, linetype = "Spline")) + 
  geom_line() +
  geom_line(data = subset(pred_modi, dose < 1), aes(y = pred_spl_spike, linetype = "Spike at 0")) +
  geom_line(data = subset(pred_modi, dose >= 1), aes(y = pred_spl_spike, linetype = "Spike at 0")) +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(trans = "log", breaks = c(.8, .85, .9, .95, 1, 1.05, 1.1, 1.15)) +
  labs(x = "Coffee consumption (cups/day)", y = "Relative risk", linetype = "Model") +
  scale_linetype_manual(values = c(`Spline` = "solid", `Spike at 0` = "dashed",
                                   `Categories` = "longdash")) +
  theme_classic()
for (i in seq_along(k2[-1])){
  p <- p + geom_line(data = subset(pred_modi, dose >= k2[i] & dose < k2[i+1]), 
                     aes(y = pred_categ, linetype = "Categories"))
}
p

# comparison AIC
lapply(list(`spline` = spl_ml, `Categories` = categ, `Spike at 0` = spl_spike), AIC)

# Meta-regression
spl_reg <- dosresmeta(logrr ~ rcs(dose, k), id = id, type = type, cases = cases, 
                      n = n, se = se, data = coffee_mort_add, covariance = "user", 
                      proc = "1stage", Slist = Slist, mod = ~ gender, method = "reml")
summary(spl_reg)
lmtest::lrtest(spl_reg, spl_ml)

# vpc plot
coffee_vpc <- coffee_mort_add %>%
  filter(se != 0) %>%
  mutate(
    `Spline` = vpc(spl),
    `Spline meta-regression` = vpc(spl_reg)
  ) %>%
  gather(curve, vpc, Spline:`Spline meta-regression`)

# Figure 6
ggplot(coffee_vpc, aes(dose, vpc, group = curve)) +
  geom_point(aes(shape = curve)) +
  geom_smooth(aes(linetype = curve), method = "loess", se = F, col = "black") +
  labs(y = "VPC", x = "Coffee consumption (cups/day)", shape = "Curve", linetype = "Curve") + 
  theme_classic() + theme(legend.position = "top")
