## ----eval = TRUE, echo = TRUE, message=FALSE, warning=FALSE--------------
library(Hmisc)
library(tidyverse)
library(ggplot2)
library(nlme)
library(contrast)

theme_set(theme_bw(base_size = 12))
cbbPalette <- c("#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
scale_colour_discrete <- function(...) scale_colour_manual(..., values=cbbPalette)
scale_fill_discrete <- function(...) scale_fill_manual(..., values=cbbPalette)

## ----eval = FALSE, echo = TRUE-------------------------------------------
## install.packages(c("tidyverse", "Hmisc", "ggplot2", "nlme", "contrast"))

## ----generate_data1, echo=TRUE-------------------------------------------
# fixed effects parameters
Beta <- c('(Intercept)'=31.6, # mean ADAS at baseline
  female=-0.63, age=0.01,     # weak effects for sex and age
  month=0.44,                 # increase in ADAS per month in controls
  'month:active'=-0.11)       # relative slowing in active group

# random effects variance parameters
sigma_random_intercept <- 7.3
sigma_random_slope <- 0.45
sigma_residual <- 3.4

# other design parameters
months <- c(0, 6, 12, 18)
n <- 200 # per group
attrition_rate <- 0.05

## ----generate_data2, echo=TRUE-------------------------------------------
set.seed(20170701)

subjects <- data.frame(
  id = 1:(2*n),
  active = sample(c(rep(0,n), rep(1,n)), 2*n),
  female = sample(0:1, 2*n, replace=TRUE),
  age = rnorm(2*n, 75, 7.8),
  censor = rexp(2*n,rate=attrition_rate),
  ran.intercept = rnorm(2*n, sd=sigma_random_intercept),
  ran.slope     = rnorm(2*n, sd=sigma_random_slope))

trial <- right_join(subjects, 
  expand.grid(id = 1:(2*n), month=months)) %>%
  mutate(
    residual = rnorm(2*n*length(months), sd=sigma_residual),
    group = factor(active, 0:1, c('placebo', 'active')),
    missing = ifelse(month>censor, 1, 0)) %>%
  arrange(id, month)

## ----generate_data3, echo=TRUE-------------------------------------------
trial$ADAS13 <- round(
  model.matrix(~ female+age+month+month:active, data = trial)[, names(Beta)] %*% 
  Beta +
  with(trial, ran.intercept + ran.slope*month + residual), 
  digits = 0
)[,1]

## ------------------------------------------------------------------------
head(trial)

## ----generate_data4, echo=TRUE-------------------------------------------
# filter out the missing observations
trial_obs <- filter(trial, !missing)

# transfrom data from long to wide
trial_wide <- trial_obs %>%
  select(id, month, female, age, active, group, ADAS13) %>% 
  mutate(month = paste0('ADAS13.m', month)) %>%
  spread(month, ADAS13)

# data for MMRM
trial_mmrm <- right_join(
  select(trial_wide, id, ADAS13.m0), 
  filter(trial_obs, month>0))

## ------------------------------------------------------------------------
head(trial_obs)

## ----results='asis'------------------------------------------------------
latex(summary(group ~ female + age + ADAS13, data = trial_obs, subset = month == 0, 
  method='reverse'), prmsd=TRUE, file='')

## ----spaghetti_plot------------------------------------------------------
ggplot(trial_obs, aes(x=month, y=ADAS13, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2)

## ------------------------------------------------------------------------
summaryTable <- trial_obs %>% 
  group_by(group, month) %>%
  summarise(
    n=length(ADAS13),
    mean=mean(ADAS13),
    sd=sd(ADAS13),
    lower95 = smean.cl.normal(ADAS13)[['Lower']],
    upper95 = smean.cl.normal(ADAS13)[['Upper']],
    min=min(ADAS13),
    max=max(ADAS13))
print(as.data.frame(summaryTable), row.names=FALSE)

## ----meanplot------------------------------------------------------------
ggplot(summaryTable, aes(x=month, y=mean, color=group)) +
  geom_line() +
  geom_errorbar(aes(min=lower95, max=upper95), position=position_dodge(0.2), width=0)

## ------------------------------------------------------------------------
m1 <- subset(summaryTable, group=='active' & month==18)[['mean']]
m2 <- subset(summaryTable, group=='placebo' & month==18)[['mean']]
n1 <- subset(summaryTable, group=='active' & month==18)[['n']]
n2 <- subset(summaryTable, group=='placebo' & month==18)[['n']]
sd1 <- subset(summaryTable, group=='active' & month==18)[['sd']]
sd2 <- subset(summaryTable, group=='placebo' & month==18)[['sd']]
s <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))
tt <- (m2-m1)/(s*sqrt(1/n2 + 1/n1)) 
DF <- n1+n2-2

## ----echo=FALSE----------------------------------------------------------
print(t.test(ADAS13 ~ group, data = trial_obs, subset = month==18, 
  var.equal=TRUE), digits = 6)

## ------------------------------------------------------------------------
x <- seq(-5, 5, by = 0.01)
dens <- data.frame(
	x        =  x,
	density  = dt(x, df = DF)
)
shadel <- rbind(
  c(x=-5, y=0),
  filter(dens, x <= -abs(tt)),
  c(x=-abs(tt),0))
shadeu <- rbind(
  c(x=abs(tt),0),
  filter(dens, x >= abs(tt)),
  c(x=5, y=0))
ggplot(dens, aes(x=x, y=density)) + geom_line() +
  geom_polygon(data = shadel, aes(x=x, y=density)) +
  geom_polygon(data = shadeu, aes(x=x, y=density))

## ----fig.show='animate'--------------------------------------------------
lmfit <- lm(ADAS13.m18 ~ ADAS13.m0, data = trial_wide)
Fitted <- predict(lmfit)

Coef <- lmfit$coef + c(-3, 0.1)*10

for(i in 1:10){
  Coef <- Coef + c(3, -0.1)
  Fitted <- as.matrix(cbind(1, filter(trial_wide, !is.na(ADAS13.m18))['ADAS13.m0'])) %*% Coef
  p <- ggplot(filter(trial_wide, !is.na(ADAS13.m18)), aes(x=ADAS13.m0, y=ADAS13.m18)) + 
    geom_point() + geom_abline(intercept=Coef[1], slope=Coef[2]) +
    xlab('ADAS13 at baseline') +
    ylab('ADAS13 at 18 months') +
    geom_segment(aes(x=ADAS13.m0, y=ADAS13.m18, 
      xend=ADAS13.m0, yend=Fitted), alpha=0.25) +
    ylim(range(trial_wide$ADAS13.m18, na.rm=TRUE))
  print(p)
}

p <- ggplot(filter(trial_wide, !is.na(ADAS13.m18)), aes(x=ADAS13.m0, y=ADAS13.m18)) + 
  geom_point() + geom_smooth(method='lm') +
  xlab('ADAS13 at baseline') +
  ylab('ADAS13 at 18 months') +
  geom_segment(aes(x=ADAS13.m0, y=ADAS13.m18, 
    xend=ADAS13.m0, yend=Fitted), alpha=0.25) +
  ylim(range(trial_wide$ADAS13.m18, na.rm=TRUE))
print(p)

## ----ancovai, echo = FALSE, size = 'scriptsize'--------------------------
summary(lm(ADAS13.m18 ~ active + ADAS13.m0, data = trial_wide))

## ----echo = FALSE, size = 'scriptsize'-----------------------------------
center <- function(x) scale(x, scale = FALSE)

## ----ancovaii, echo = FALSE, size = 'scriptsize'-------------------------
summary(lm(ADAS13.m18 ~ active*center(ADAS13.m0), data = trial_wide))

## ----ancovaii2cov, echo = FALSE, size = 'scriptsize'---------------------
summary(lm(ADAS13.m18 ~ active*center(ADAS13.m0) + female + age, data = trial_wide))

## ----trial_stage1_plot, echo = FALSE-------------------------------------
ggplot(subset(trial, id %in% 1:4),
  aes(x=month, y=ADAS13, group = id, color = group)) +
  stat_smooth(method = 'lm') + geom_line() + geom_point() +
  facet_wrap(~id)

## ----trial_fit_stage1, eval = TRUE, echo = FALSE, size = 'scriptsize'----
trial_stage1 <- as.data.frame(do.call(rbind, lapply(unique(trial$id),
  function(i){
    fit <- lm(ADAS13 ~ month,
      data = trial_obs, subset = id == i)
    c(id = i, beta = fit$coef, sigma = summary(fit)$sigma)
})))
trial_stage1 <- right_join(trial_stage1,
  filter(trial, month == 0) %>% 
    select(id, active, group, age, female))
head(trial_stage1)

## ----size = 'scriptsize'-------------------------------------------------
summary(lm(ADAS13 ~ month, data = trial_obs, subset = id == 1))

## ----trial_plot_stage2---------------------------------------------------
ggplot(trial_stage1,
  aes(x=group, y=beta.month, group = group, color = group)) +
  geom_boxplot(alpha = 0.25) + 
  ylab('ADAS13 change per month')

## ----trial_stage2_fit----------------------------------------------------
summary(lm(beta.month ~ female + age + active, data = trial_stage1))

## ----trial_lme, size = 'tiny'--------------------------------------------
fit_lme <- lme(ADAS13 ~ month + month:active, data = trial_obs, random = ~month|id)
summary(fit_lme)

## ----trial_lme_age, size = 'tiny'----------------------------------------
fit_lme_cov <- lme(ADAS13 ~ age + female + month + month:active, data = trial_obs, random = ~month|id)
summary(fit_lme_cov)

## ----trial_lme_rcode, eval = FALSE, echo = TRUE--------------------------
## lme(ADAS13 ~ month + month:active, data = trial_obs,
##   random = ~month|id)
## 
## lme(ADAS13 ~ age + female + month + month:active, data = trial_obs,
##   random = ~month|id)

## ----trial_lme_profiles, echo = FALSE, size = 'scriptsize'---------------
trial_obs <- mutate(trial_obs, month.active = month*active)
fit_lme_cov_plot <- lme(ADAS13 ~ age + female + month + month.active, data = trial_obs, random = ~month|id)

# get profile for each apoe group
act_profile <- contrast(fit_lme_cov_plot,
  a = list(age = mean(trial_wide$age),
           female = 1,
           month = seq(0,18,6),
           month.active = seq(0,18,6)
))
pbo_profile <- contrast(fit_lme_cov_plot,
  a = list(age = mean(trial_wide$age),
           female = 1,
           month = seq(0,18,6),
           month.active = 0
))
# combine the profiles into one data.frame
pd <- rbind(
  filter(with(act_profile,
    data.frame(active=1, month.active, month, age, Estimate = Contrast, Lower, Upper)),
  month.active == month),
  with(pbo_profile,
    data.frame(active=0, month.active, month, age, Estimate = Contrast, Lower, Upper))
)

pd$group <- factor(pd$active, 0:1, c('placebo', 'active'))

## ----echo = FALSE, size = 'scriptsize'-----------------------------------
print(pd, row.names=FALSE)

## ----echo = FALSE, size = 'scriptsize'-----------------------------------
ggplot(pd, aes(x = month, y = Estimate, group = group))+
  geom_line(aes(color=group)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill=group), alpha=0.25) +
  ylim(c(25,45))

## ----echo = FALSE, size = 'scriptsize'-----------------------------------
ggplot(pd, aes(x = month, y = Estimate, group = group, color=group))+
  geom_line() +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0, position=position_dodge(0.2)) +
  ylim(c(25,45))

## ----varPar, echo=FALSE--------------------------------------------------
# simulate data with different variance parameters
varPar <- expand.grid(
  sigma_random_intercept = c(2, 10),
  sigma_random_slope = c(0.45, 0.75),
  sigma_residual = c(2, 8)
)

set.seed(20170701)

varPlot <- do.call(rbind, lapply(1:nrow(varPar), function(i){
  subjects <- data.frame(
    id = 1:(2*n),
    active = sample(c(rep(0,n), rep(1,n)), 2*n),
    female = sample(0:1, 2*n, replace=TRUE),
    age = rnorm(2*n, 75, 7.8),
    censor = rexp(2*n,rate=attrition_rate),
    sigma_random_intercept = varPar[i, 'sigma_random_intercept'],
    sigma_random_slope = varPar[i, 'sigma_random_slope'],
    sigma_residual = varPar[i, 'sigma_residual'],
    ran.intercept = rnorm(2*n, sd=varPar[i, 'sigma_random_intercept']),
    ran.slope     = rnorm(2*n, sd=varPar[i, 'sigma_random_slope']))

  trial <- right_join(subjects, 
    expand.grid(id = 1:(2*n), month=months)) %>%
    mutate(
      residual = rnorm(2*n*length(months), sd=varPar[i, 'sigma_residual']),
      group = factor(active, 0:1, c('placebo', 'active')),
      missing = ifelse(month>censor, 1, 0)) %>%
    arrange(id, month) %>%
    filter(!missing)
  trial$ADAS13 <- round(
    model.matrix(~ female+age+month+month:active, data = trial)[, names(Beta)] %*% 
    Beta +
    with(trial, ran.intercept + ran.slope*month + residual), 
    digits = 0
  )[,1]
  trial
}))

## ------------------------------------------------------------------------
ggplot(filter(varPlot, sigma_random_intercept==2 & sigma_random_slope==0.45), 
  aes(x=month, y=ADAS13, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~sigma_residual) +
  ylim(0,70)

## ------------------------------------------------------------------------
ggplot(filter(varPlot, sigma_residual==2 & sigma_random_slope==0.45), 
  aes(x=month, y=ADAS13, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~sigma_random_intercept) +
  ylim(0,70)

## ------------------------------------------------------------------------
ggplot(filter(varPlot, sigma_residual==2 & sigma_random_intercept==2), 
  aes(x=month, y=ADAS13, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~sigma_random_slope) +
  ylim(0,70)

## ----trial_lme_apoe_int, size = 'tiny'-----------------------------------
fit_lme_int <- update(fit_lme, random = ~1|id)
summary(fit_lme_int)

## ----trial_lme_apoe_int_vs_slope, size = 'footnotesize', echo = TRUE-----
anova(fit_lme_int, fit_lme)

