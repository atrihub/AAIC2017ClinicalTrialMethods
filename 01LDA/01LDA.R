options(digits = 3, stringsAsFactors = FALSE, width = 150)

## ----eval = TRUE, echo = FALSE, message=FALSE, warning=FALSE-------------
# required R packages
library(Hmisc)
library(multcomp)
library(tidyverse)
library(ggplot2)
library(nlme)
library(grid)
library(gridExtra)

# ggplot theme and palette settings
theme_set(theme_bw(base_size = 12))
cbbPalette <- c('#0072B2', '#D55E00', '#CC79A7', '#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442')
scale_colour_discrete <- function(...) scale_colour_manual(..., values=cbbPalette)
scale_fill_discrete <- function(...) scale_fill_manual(..., values=cbbPalette)

theme_table <- function(..., levs=2){
  theme_minimal(...) + 
    theme(
      panel.grid = element_blank(), 
      axis.text.x = element_blank(),
      axis.text.y = element_text(face='bold', color=cbbPalette[1:levs]),
      axis.title = element_blank())
} 

## ----eval = FALSE, echo = FALSE------------------------------------------
## # If you do not have them installed already, you will need to download them from CRAN via:
## install.packages(c('tidyverse', 'Hmisc', 'ggplot2', 'nlme', 'contrast', 'grid', 'gridExtra'))

## ----echo=FALSE, eval=FALSE----------------------------------------------
## ####################################################
## ### pilot estimates are from a model fit to ADNI ###
## ####################################################
## # library(ADNIMERGE) # available for loni.usc.edu
## fit_adni <- lme(ADAS11 ~ PTGENDER + scale(AGE,scale=FALSE) + M, data=adnimerge,
##   random=~M|RID, subset = M<=24 & DX.bl=='AD', na.action=na.omit)
## summary(fit_adni)

## ----generate_data1, echo=TRUE-------------------------------------------
# fixed effects parameters estimated from ADNI
Beta <- c(
   '(Intercept)'=19.52, # mean ADAS at baseline
        'female'=-0.29, # better scores for females
         'age_c'= 0.06, # worse change for older at baseline (age mean centered)
         'month'= 0.42, # worse per month post baseline
  'month:active'=-0.05) # improvement per month with treatment

# standard deviations for random effects
sigma_random_intercept <- 6.1
sigma_random_slope <- 0.38
sigma_residual <- 3.3

# other design parameters
months <- c(0, 6, 12, 18)
n <- 200 # per group
attrition_rate <- 0.40/18 # approx per month

## ------------------------------------------------------------------------
# set seed so that simulation is reproducible
set.seed(20170701)

# simulate subject specific data
subjects <- data.frame(
  id = 1:(2*n),
  active = sample(c(rep(0,n), rep(1,n)), 2*n),
  female = sample(0:1, 2*n, replace=TRUE),
  age = rnorm(2*n, 75, 7.8), 
  censor = rexp(2*n,rate=attrition_rate),
  ran.intercept = rnorm(2*n, sd=sigma_random_intercept),
  ran.slope     = rnorm(2*n, sd=sigma_random_slope)) %>%
  mutate(age_c = age-mean(age))

# simulate data over time
trial <- right_join(subjects, 
  expand.grid(id = 1:(2*n), month=months)) %>%
  mutate(
    residual = rnorm(2*n*length(months), sd=sigma_residual),
    group = factor(active, 0:1, c('placebo', 'active')),
    missing = ifelse(month>censor, 1, 0)) %>%
  arrange(id, month)

# calculate the ADAS scores with random effects and residuals and 
# round to nearest digit in 0-70
trial <- mutate(trial,
  ADAS11 = (model.matrix(~ female+age_c+month+month:active, trial)[, names(Beta)] %*% Beta)[,1],
  ADAS11 = round(ADAS11 + ran.intercept + ran.slope*month + residual, digits = 0),
  ADAS11 = replace(ADAS11, ADAS11<0, 0),
  ADAS11 = replace(ADAS11, ADAS11>70, 70))

# filter out the missing observations
trial_obs <- filter(trial, !missing)

# transfrom data from long to wide
trial_wide <- trial_obs %>%
  select(id, month, female, age, age_c, active, group, ADAS11) %>% 
  mutate(month = paste0('ADAS11.m', month)) %>%
  spread(month, ADAS11)

# data for MMRM
trial_mmrm <- right_join(
  select(trial_wide, id, ADAS11.m0), 
  filter(trial_obs, month>0)) %>%
  mutate(ADAS11.ch = ADAS11 - ADAS11.m0)

## ------------------------------------------------------------------------
head(select(trial_obs, -group, -missing))

## ----results='asis'------------------------------------------------------
latex(summary(group ~ female + age + ADAS11, data = trial_obs, subset = month == 0, 
  method='reverse'), prmsd=TRUE, file='')

## ----spaghetti_plot------------------------------------------------------
ggplot(trial_obs, aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.1, 0.85))

## ------------------------------------------------------------------------
summaryTable <- trial_obs %>% 
  group_by(group, month) %>%
  summarise(
    n=length(ADAS11),
    mean=mean(ADAS11),
    sd=sd(ADAS11),
    lower95 = smean.cl.normal(ADAS11)[['Lower']],
    upper95 = smean.cl.normal(ADAS11)[['Upper']],
    min=min(ADAS11),
    max=max(ADAS11))
print(as.data.frame(summaryTable), row.names=FALSE)

## ----meanplot------------------------------------------------------------
p <- ggplot(summaryTable, aes(x=month, y=mean, color=group)) +
  geom_line() +
  geom_errorbar(aes(min=lower95, max=upper95), position=position_dodge(0.2), width=0) +
  scale_x_continuous(breaks=months) +
  ylab('Mean ADAS (95% CI)') +
  theme(legend.position=c(0.1, 0.75))
countTab <- ggplot(summaryTable, aes(x=month, y=group, label=n)) + geom_text() + theme_table()
grid.draw(arrangeGrob(p,countTab,heights=c(3,1)))

## ------------------------------------------------------------------------
chsummaryTable <- trial_mmrm %>% 
  group_by(group, month) %>%
  summarise(
    n=length(ADAS11.ch),
    mean=mean(ADAS11.ch),
    sd=sd(ADAS11.ch),
    lower95 = smean.cl.normal(ADAS11.ch)[['Lower']],
    upper95 = smean.cl.normal(ADAS11.ch)[['Upper']],
    min=min(ADAS11.ch),
    max=max(ADAS11.ch))
chsummaryTable <- bind_rows(
  filter(summaryTable, month == 0) %>% 
    mutate(mean=0, sd=0, lower95=0, upper95=0, min=0, max=0),
  chsummaryTable) %>%
  arrange(group, month)

print(as.data.frame(chsummaryTable), row.names=FALSE)

## ----meanchplot----------------------------------------------------------
p <- ggplot(chsummaryTable, aes(x=month, y=mean, color=group)) +
  geom_line() +
  geom_errorbar(aes(min=lower95, max=upper95), position=position_dodge(0.2), width=0) +
  scale_x_continuous(breaks=months) +
  ylab('Mean ADAS change') +
  theme(legend.position=c(0.1, 0.75))
grid.draw(arrangeGrob(p,countTab,heights=c(3,1)))

## ------------------------------------------------------------------------
m1 <- subset(chsummaryTable, group=='active' & month==18)[['mean']]
m2 <- subset(chsummaryTable, group=='placebo' & month==18)[['mean']]
n1 <- subset(chsummaryTable, group=='active' & month==18)[['n']]
n2 <- subset(chsummaryTable, group=='placebo' & month==18)[['n']]
sd1 <- subset(chsummaryTable, group=='active' & month==18)[['sd']]
sd2 <- subset(chsummaryTable, group=='placebo' & month==18)[['sd']]
s <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))
tt <- (m2-m1)/(s*sqrt(1/n2 + 1/n1)) 
DF <- n1+n2-2

## ----echo=FALSE----------------------------------------------------------
print(t.test(ADAS11.ch ~ group, data = trial_mmrm, subset = month==18, 
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
  geom_polygon(data = shadel, aes(x=x, y=density), fill=cbbPalette[1]) +
  geom_polygon(data = shadeu, aes(x=x, y=density), fill=cbbPalette[1]) +
  geom_vline(xintercept=qt(c(0.025, 0.975), df=DF), color='grey', linetype='dashed') +
  annotate('text', x=-1.96, y=0.4, label='x=-1.96') +
  annotate('text', x=1.96, y=0.4, label='x=1.96')

## ----fig.show='animate'--------------------------------------------------
lmfit <- lm(ADAS11.m18 ~ ADAS11.m0, data = trial_wide)
Fitted <- predict(lmfit)

Coef <- lmfit$coef + c(-3, 0.1)*10

for(i in 1:10){
  Coef <- Coef + c(3, -0.1)
  Fitted <- as.matrix(cbind(1, filter(trial_wide, !is.na(ADAS11.m18))['ADAS11.m0'])) %*% Coef
  p <- ggplot(filter(trial_wide, !is.na(ADAS11.m18)), aes(x=ADAS11.m0, y=ADAS11.m18)) + 
    geom_point() + geom_abline(intercept=Coef[1], slope=Coef[2]) +
    xlab('ADAS11 at baseline') +
    ylab('ADAS11 at 18 months') +
    geom_segment(aes(x=ADAS11.m0, y=ADAS11.m18, 
      xend=ADAS11.m0, yend=Fitted), alpha=0.25) +
    coord_cartesian(ylim=c(range(trial_wide$ADAS11.m18, na.rm=TRUE)))
  print(p)
}

p <- ggplot(filter(trial_wide, !is.na(ADAS11.m18)), aes(x=ADAS11.m0, y=ADAS11.m18)) + 
  geom_point() + geom_smooth(method='lm') +
  xlab('ADAS11 at baseline') +
  ylab('ADAS11 at 18 months') +
  geom_segment(aes(x=ADAS11.m0, y=ADAS11.m18, 
    xend=ADAS11.m0, yend=Fitted), alpha=0.25) +
  ylim(range(trial_wide$ADAS11.m18, na.rm=TRUE))
print(p)

## ----ancovai, echo = FALSE, size = 'scriptsize'--------------------------
summary(lm(ADAS11.ch ~ active + ADAS11.m0, data = trial_mmrm))

## ----echo = FALSE, size = 'scriptsize'-----------------------------------
center <- function(x) scale(x, scale = FALSE)

## ----ancovaii, echo = FALSE, size = 'scriptsize'-------------------------
summary(lm(ADAS11.ch ~ active*center(ADAS11.m0), data = trial_mmrm))

## ----ancovaii2cov, echo = FALSE, size = 'scriptsize'---------------------
summary(lm(ADAS11.ch ~ active*center(ADAS11.m0) + female + age_c, data = trial_mmrm))

## ----trial_stage1_plot, echo = FALSE-------------------------------------
ggplot(subset(trial, id %in% 1:4),
  aes(x=month, y=ADAS11, group = id, color = group)) +
  stat_smooth(method = 'lm') + geom_line() + geom_point() +
  facet_wrap(~id)

## ----trial_fit_stage1, eval = TRUE, echo = FALSE, size = 'scriptsize'----
trial_stage1 <- as.data.frame(do.call(rbind, lapply(unique(trial$id),
  function(i){
    fit <- lm(ADAS11 ~ month,
      data = trial_obs, subset = id == i)
    c(id = i, beta = fit$coef, sigma = summary(fit)$sigma)
})))
trial_stage1 <- right_join(trial_stage1,
  filter(trial, month == 0) %>% 
    select(id, active, group, age_c, female))
head(trial_stage1)

## ----size = 'scriptsize'-------------------------------------------------
summary(lm(ADAS11 ~ month, data = trial_obs, subset = id == 1))

## ----trial_plot_stage2---------------------------------------------------
ggplot(trial_stage1,
  aes(x=group, y=beta.month, group = group, color = group)) +
  geom_boxplot(alpha = 0.25) + 
  ylab('ADAS11 change per month')

## ----trial_stage2_fit----------------------------------------------------
summary(lm(beta.month ~ female + age_c + active, data = trial_stage1))

## ----trial_lme, size = 'tiny'--------------------------------------------
fit_lme <- lme(ADAS11 ~ month + month:active, data = trial_obs, random = ~month|id)
summary(fit_lme)

## ----trial_lme_age, size = 'tiny'----------------------------------------
fit_lme_cov <- lme(ADAS11 ~ age_c + female + month + month:active, data = trial_obs, random = ~month|id)
summary(fit_lme_cov)

## ----trial_lme_rcode, eval = FALSE, echo = TRUE--------------------------
## lme(ADAS11 ~ month + month:active,
##   data = trial_obs, random = ~month|id)
## 
## lme(ADAS11 ~ age_c + female + month + month:active,
##   data = trial_obs, random = ~month|id)

## ----trial_lme_profiles, echo = FALSE, size = 'scriptsize'---------------
plotData0 <- filter(trial_obs, !duplicated(paste(month, active))) %>%
  arrange(active, month) %>%
  mutate(female = 1, age_c=0) %>%
  select(-age, -censor, -residual, -ADAS11)

plotMatrix <- model.matrix(~ age_c + female + month + month:active, data = plotData0)

plotData <- bind_cols(plotData0, 
  as.data.frame(confint(glht(fit_lme_cov, linfct = plotMatrix))$confint))


## ----echo = FALSE, size = 'scriptsize'-----------------------------------
print(select(plotData, group, month, Estimate, lwr, upr), row.names=FALSE)

## ----echo = FALSE, size = 'scriptsize'-----------------------------------
p <- ggplot(plotData, aes(x = month, y = Estimate, group = group))+
  geom_line(aes(color=group)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill=group), alpha=0.25) +
  ylim(c(15,30)) +
  scale_x_continuous(breaks=months) +
  ylab('Mean ADAS (95% CI)') +
  theme(legend.position=c(0.1, 0.75))
grid.draw(arrangeGrob(p,countTab,heights=c(3,1)))

## ----echo = FALSE, size = 'scriptsize'-----------------------------------
p <- ggplot(plotData, aes(x = month, y = Estimate, group = group, color=group))+
  geom_line() +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width=0, position=position_dodge(0.2)) +
  ylab('Mean ADAS (95% CI)') +
  ylim(c(15,30)) +
  scale_x_continuous(breaks=months) +
  ylab('Mean ADAS (95% CI)') +
  theme(legend.position=c(0.1, 0.75))
grid.draw(arrangeGrob(p,countTab,heights=c(3,1)))

## ----varPar, echo=FALSE--------------------------------------------------
# simulate data with different variance parameters
varPar <- expand.grid(
  sigma_random_intercept = c(2, 10),
  sigma_random_slope = c(0.2, 0.75),
  sigma_residual = c(2, 8)
)

varPlot <- do.call(rbind, lapply(1:nrow(varPar), function(i){
  set.seed(20170714)
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
    ran.slope     = rnorm(2*n, sd=varPar[i, 'sigma_random_slope'])) %>%
    mutate(age_c = age - mean(age))

  trial <- right_join(subjects, 
    expand.grid(id = 1:(2*n), month=months)) %>%
    mutate(
      residual = rnorm(2*n*length(months), sd=varPar[i, 'sigma_residual']),
      group = factor(active, 0:1, c('placebo', 'active')),
      missing = ifelse(month>censor, 1, 0)) %>%
    arrange(id, month) %>%
    filter(!missing)
  trial$ADAS11 <- round(
    model.matrix(~ female+age_c+month+month:active, data = trial)[, names(Beta)] %*% 
    Beta +
    with(trial, ran.intercept + ran.slope*month + residual), 
    digits = 0
  )[,1]
  trial
}))

## ------------------------------------------------------------------------
ggplot(filter(varPlot, sigma_random_intercept==2 & sigma_random_slope==0.2), 
  aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~paste('Residual SD =', sigma_residual)) +
  ylim(0,70) +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.1,0.8))

## ------------------------------------------------------------------------
ggplot(filter(varPlot, sigma_residual==2 & sigma_random_slope==0.2), 
  aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~paste('Random intercept SD =', sigma_random_intercept)) +
  ylim(0,70) +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.6,0.8))

## ------------------------------------------------------------------------
ggplot(filter(varPlot, sigma_residual==2 & sigma_random_intercept==2), 
  aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~paste('Random slope SD =', sigma_random_slope)) +
  ylim(0,70) +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.1,0.8))

## ----trial_lme_apoe_int, size = 'tiny'-----------------------------------
fit_lme_int <- update(fit_lme, random = ~1|id)
summary(fit_lme_int)

## ----trial_lme_apoe_int_vs_slope, size = 'footnotesize', echo = TRUE-----
anova(fit_lme_int, fit_lme)

