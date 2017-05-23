## ----packages, results = 'hide', echo = FALSE, warning = FALSE, message = FALSE, cache = FALSE, results = 'hide'----
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
library(mvtnorm)
library(mice)

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

center <- function(x){scale(x, scale=FALSE)} 

## ------------------------------------------------------------------------
# fixed effects parameters estimated from ADNI
Beta <- c(
   '(Intercept)'=19.60, # mean ADAS at baseline
        'female'=-0.78, # better scores for females
           'age'= 0.01, # worse per year of age at baseline
         'month'= 0.40, # worse per month post baseline
  'month:active'=-0.05) # improvement per month with treatment

# standard deviations for random effects
sigma_random_intercept <- 6.0
sigma_random_slope <- 0.42
sigma_residual <- 3.1

# other design parameters
months <- c(0, 6, 12, 18)
n <- 200 # per group
attrition_rate <- 0.40/18 # approx per month

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
  ran.slope     = rnorm(2*n, sd=sigma_random_slope))

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
  ADAS11 = (model.matrix(~ female+age+month+month:active, trial)[, names(Beta)] %*% Beta)[,1],
  ADAS11 = round(ADAS11 + ran.intercept + ran.slope*month + residual, digits = 0),
  ADAS11 = replace(ADAS11, ADAS11<0, 0),
  ADAS11 = replace(ADAS11, ADAS11>70, 70))

# filter out the missing observations
trial_mmrm <- filter(trial, !missing)

# transfrom data from long to wide
trial_wide <- trial_mmrm %>%
  select(id, month, female, age, active, group, ADAS11) %>% 
  mutate(month = paste0('ADAS11.m', month)) %>%
  spread(month, ADAS11) %>%
  select(id:group, ADAS11.m0, ADAS11.m6, ADAS11.m12, ADAS11.m18)

# data for MMRM
trial_mmrm <- right_join(
  select(trial_wide, id, ADAS11.m0), 
  filter(trial_mmrm, month>0)) %>%
  mutate(ADAS11.ch = ADAS11 - ADAS11.m0)

## ------------------------------------------------------------------------
means <- expand.grid(female=1, age=70, month=0:18, active=0, id=1:3) %>%
  filter(id %in% c(1,2) | month %in% months) 
means <- mutate(means,
    ADAS11 = (model.matrix(~ female+age+month+month:active, means)[, names(Beta)] %*% Beta)[,1],
    ADAS11 = replace(ADAS11, ADAS11<0, 0),
    ADAS11 = replace(ADAS11, ADAS11>70, 70),
    ADAS11 = replace(ADAS11, id == 2, ADAS11 + month*month*0.1),
    ADAS11 = replace(ADAS11, id == 3 & month == 0, ADAS11 + 0),
    ADAS11 = replace(ADAS11, id == 3 & month == 6, ADAS11 + 10),
    ADAS11 = replace(ADAS11, id == 3 & month == 12, ADAS11 + 2),
    ADAS11 = replace(ADAS11, id == 3 & month == 18, ADAS11 + 25)) %>%
  filter(id %in% 1:3) %>%
  mutate(Mean = factor(id, labels=c('linear', 'quadratic', 'categorical')))

ggplot(means, aes(x=month, y=ADAS11, group=Mean, color=Mean)) + 
  geom_line() +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.15,0.8))

## ----varPar, echo=FALSE--------------------------------------------------
# simulate data with different variance parameters
varPar <- expand.grid(
  variance = c('homogeneous', 'heterogeneous'),
  correlation = c('uncorrelated', 'correlated')
)

SD <- list(homogeneous = sqrt(rep(4, 4)), heterogeneous = (1:4)*2)
Cor <- list(uncorrelated = 0, correlated = 0.8)

varPlot <- do.call(rbind, lapply(1:nrow(varPar), function(i){
  set.seed(20170714)
  subjects <- data.frame(
    id = 1:(2*n),
    active = sample(c(rep(0,n), rep(1,n)), 2*n),
    female = sample(0:1, 2*n, replace=TRUE),
    age = rnorm(2*n, 75, 7.8),
    censor = rexp(2*n,rate=attrition_rate))
    
  vv <- diag(SD[[varPar[i,'variance']]])
  cc <- matrix(Cor[[varPar[i,'correlation']]], nrow=4, ncol=4)
  diag(cc) <- 1
  resids <- as.numeric(t(rmvnorm(nrow(subjects), mean=rep(0,4), sigma=vv%*%cc%*%vv)))

  trial <- right_join(subjects, 
    expand.grid(id = 1:(2*n), month=months)) %>%
    arrange(id, month) %>%
    mutate(residual = resids,
      group = factor(active, 0:1, c('placebo', 'active')),
      missing = ifelse(month>censor, 1, 0),
      variance = varPar[i,'variance'],
      correlation = varPar[i,'correlation']) %>%
    arrange(id, month) %>%
    filter(!missing)
  trial$ADAS11 <- round(
    model.matrix(~ female+age+month+month:active, data = trial)[, names(Beta)] %*% 
    Beta + trial$residual, 
    digits = 0
  )[,1]
  trial
}))

## ------------------------------------------------------------------------
ggplot(filter(varPlot, correlation=='correlated'), 
  aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~variance) +
  ylim(0,50) +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.1,0.8))

## ------------------------------------------------------------------------
ggplot(filter(varPlot, variance=='heterogeneous'), 
  aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~correlation) +
  ylim(0,50) +
  scale_x_continuous(breaks=months) +
  theme(legend.position='none')

## ----echo=FALSE, eval=FALSE----------------------------------------------
## ####################################################
## ### pilot estimates are from a model fit to ADNI ###
## ####################################################
## 
## # library(ADNIMERGE) # available for loni.usc.edu
## # adni_ad <- filter(adnimerge, M<=24 & M!=18 & !is.na(ADAS11) & DX.bl=='AD') %>%
## #   mutate(m = as.factor(M),
## #      visNo = as.numeric(m))
## #
## # with(adni_ad, table(m, visNo))
## # fit_adni <- gls(ADAS11 ~ PTGENDER + scale(AGE, scale=FALSE) + m, data=adni_ad,
## #   correlation = corSymm(form = ~ visNo | RID),
## #   weights = varIdent(form = ~ 1 | m) )
## # summary(fit_adni)

## ----echo=TRUE-----------------------------------------------------------
Beta <- c(
   '(Intercept)'= 19.8, # mean ADAS at baseline
        'female'=-0.51, # female perform better
         'age_c'= 0.04, # worse change for older at baseline (age mean centered)
            'm6'= 2.23, # worsening at month 6 in pbo
           'm12'= 4.46, # worsening at month 12 in pbo
           'm18'= 7.31, # worsening at month 18 in pbo
     'm6:active'=-0.20, # relative improvement at month 6 with treatment
    'm12:active'=-0.70, # relative improvement at month 12 with treatment
    'm18:active'=-1.75) # relative improvement at month 18 with treatment

# other design parameters
months <- c(0, 6, 12, 18)
n <- 200 # per group
attrition_rate <- 0.40/18 # approx per month

# var-cov parameters
SD <- 6.77                          # standard deviation scale parameter
vv <- diag(c(1, 1.2, 1.5, 1.8))          # hetergeneous variance weight matrix
cc <- matrix(0.75, nrow=4, ncol=4)   # correlation matrix
diag(cc) <- 1

## ------------------------------------------------------------------------
# set seed so that simulation is reproducible
set.seed(20170714)

# simulate subject specific data
subjects <- data.frame(
  id = 1:(2*n),
  active = sample(c(rep(0,n), rep(1,n)), 2*n),
  female = sample(0:1, 2*n, replace=TRUE),
  age = rnorm(2*n, 75, 7.8), 
  censor = rexp(2*n,rate=attrition_rate)) %>%
  mutate(age_c = age-mean(age))
  
# simulate vector of correlated residuals
resids <- as.numeric(t(rmvnorm(nrow(subjects), mean=rep(0,nrow(vv)), sigma=SD^2*vv%*%cc%*%vv)))

# simulate data over time
trial <- right_join(subjects,
  expand.grid(id = 1:(2*n), month=months)) %>%
  arrange(id, month) %>%    ## WARNING: data must be properly sorted by subject and time 
  mutate(residual = resids, ## prior to appending residuals
    group = factor(active, 0:1, c('placebo', 'active')),
    missing = ifelse(month>censor, 1, 0),
    m = as.factor(month),
    visNo = as.numeric(m)) %>%
  arrange(id, month)

# create visit indicators
trial <- cbind(trial, model.matrix(~ -1+m, data = trial))

# calculate the ADAS scores with random effects and residuals and 
# round to nearest digit in 0-70
trial <- mutate(trial,
  ADAS11 = (model.matrix(~ female+age_c+m6+m12+m18+(m6+m12+m18):active, data = trial)[, names(Beta)] %*% Beta)[,1],
  ADAS11 = round(ADAS11 + residual, digits = 0),
  ADAS11 = replace(ADAS11, ADAS11<0, 0),
  ADAS11 = replace(ADAS11, ADAS11>70, 70))

# filter out the missing observations
trial_obs <- filter(trial, !missing)

# transfrom data from long to wide
trial_wide <- trial_obs %>%
  select(id, month, female, age, age_c, active, group, ADAS11) %>% 
  mutate(month = paste0('ADAS11.m', month)) %>%
  spread(month, ADAS11) %>%
  select(id:group, ADAS11.m0, ADAS11.m6, ADAS11.m12, ADAS11.m18)

# data for MMRM
trial_mmrm <- right_join(
  select(trial_wide, id, ADAS11.m0), 
  filter(trial_obs, month>0)) %>%
  mutate(ADAS11.ch = ADAS11 - ADAS11.m0,
    m = as.factor(month),
    visNo = as.numeric(m))

## ------------------------------------------------------------------------
head(select(trial_mmrm, -group, -missing, -m))

## ----spaghetti_plot------------------------------------------------------
ggplot(trial, aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'loess', size = 2) +
  scale_x_continuous(breaks=months, lim=c(0,18)) +
  theme(legend.position=c(0.1, 0.85))

## ----echo = TRUE---------------------------------------------------------
# Symmetric correlation, hetergeneous variance
MMRMsymHet <- gls(ADAS11.ch ~ 
  -1+ADAS11.m0+female+age_c+(m6+m12+m18)+(m6+m12+m18):active,
  data=trial_mmrm, correlation = corSymm(form = ~ visNo | id),
  weights = varIdent(form = ~ 1 | m) )

# Compound Symmetric correlation, hetergeneous variance
MMRMcompSymHet <- gls(ADAS11.ch ~ 
  -1+ADAS11.m0+female+age_c+(m6+m12+m18)+(m6+m12+m18):active,
  data=trial_mmrm, correlation = corCompSymm(form = ~ visNo | id),
  weights = varIdent(form = ~ 1 | m) )

# see ?corClasses and ?varClasses for more options

## ----echo=TRUE, eval=FALSE-----------------------------------------------
## summary(MMRMsymHet)

## ------------------------------------------------------------------------
x <- summary(MMRMsymHet)
print(summary(x$modelStruct), sigma = x$sigma,
reEstimates = x$coef$random, verbose = verbose)

## ------------------------------------------------------------------------
cat('Coefficients:\n')
xtTab <- as.data.frame(x$tTable)
printCoefmat(x$tTable, eps=0.001, digits=3)
# cat("\nStandardized residuals:\n")
# print(x$residuals)
cat("\n")
cat("Residual standard error:", format(x$sigma),"\n")
# cat("Degrees of freedom:", dd[["N"]],"total;",dd[["N"]] - dd[["p"]],
#     "residual\n")

## ----echo = TRUE---------------------------------------------------------
# Symmetric correlation, hetergeneous variance
cLDAsymHet <- gls(ADAS11 ~ 
  -1+female+age_c+m0+(m6+m12+m18)+(m6+m12+m18):active,
  data=trial_obs, correlation = corSymm(form = ~ visNo | id),
  weights = varIdent(form = ~ 1 | m) )

# Compound Symmetric correlation, hetergeneous variance
cLDAcompSymHet <- gls(ADAS11 ~ 
  -1+female+age_c+m0+(m6+m12+m18)+(m6+m12+m18):active,
  data=trial_obs, correlation = corCompSymm(form = ~ visNo | id),
  weights = varIdent(form = ~ 1 | m) )

# see ?corClasses and ?varClasses for more options

## ------------------------------------------------------------------------
x <- summary(cLDAsymHet)
print(summary(x$modelStruct), sigma = x$sigma,
reEstimates = x$coef$random, verbose = verbose)

## ------------------------------------------------------------------------
cat('Coefficients:\n')
xtTab <- as.data.frame(x$tTable)
printCoefmat(x$tTable, eps=0.001, digits=3)
# cat("\nStandardized residuals:\n")
# print(x$residuals)
cat("\n")
cat("Residual standard error:", format(x$sigma),"\n")
# cat("Degrees of freedom:", dd[["N"]],"total;",dd[["N"]] - dd[["p"]],
#     "residual\n")

## ----echo = TRUE---------------------------------------------------------
# Linear time, random intercept
cLDAlin <- lme(ADAS11 ~ 
  female + age_c + month + month:active,
  data=trial_obs, random = ~1|id)

# Quadratic time, random intercept
cLDAquad <- lme(ADAS11 ~ 
  female + age_c + (month + I(month^2)) + (month + I(month^2)):active,
  data=trial_obs, random = ~1|id)

## ----echo = FALSE, size = 'scriptsize'-----------------------------------
plotData0 <- filter(trial_obs, !duplicated(paste(month, active))) %>%
  arrange(active, month) %>%
  mutate(female = 1, age_c=0) %>%
  select(-age, -censor, -residual, -ADAS11)

plotMatrix_quad <- model.matrix(~ female + age_c + (month + I(month^2)) + (month + I(month^2)):active, 
  data = plotData0)

plotMatrix_cat <- model.matrix(~ -1 + female + age_c + m0 + (m6 + m12 + m18) + (m6 + m12 + m18):active, 
  data = plotData0)

plotData <- bind_rows(
  bind_cols(plotData0, as.data.frame(confint(glht(cLDAquad, linfct = plotMatrix_quad))$confint)) %>%
    mutate(model = 'quadratic time'),
  bind_cols(plotData0, as.data.frame(confint(glht(cLDAsymHet, linfct = plotMatrix_cat))$confint)) %>%
    mutate(model = 'categorical time'))

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
countTab <- ggplot(summaryTable, aes(x=month, y=group, label=n)) + geom_text() + theme_table()

p <- ggplot(plotData, aes(x = month, y = Estimate))+
  geom_line(aes(color=group, linetype=model)) +
  geom_point(aes(color=group)) +
  ylim(c(15,30)) +
  scale_x_continuous(breaks=months) +
  ylab('Mean ADAS (95% CI)') +
  theme(legend.position='top')
grid.draw(arrangeGrob(p,countTab,heights=c(3,1)))

## ------------------------------------------------------------------------
contrastData0 <- filter(trial_obs, !duplicated(paste(month, active))) %>%
  arrange(active, month) %>%
  filter(active==1 & month>0) %>%
  mutate(female = 1, age_c=0) %>%
  select(-age, -censor, -residual, -ADAS11)

contrastMatrix_quad <- as.data.frame(model.matrix(~ female + age_c + (month + I(month^2)) + (month + I(month^2)):active, 
  data = contrastData0)) %>%
  mutate('(Intercept)' = 0, female=0, month=0, 'I(month^2)' = 0)

contrastMatrix_cat <- as.data.frame(model.matrix(~ -1 + female + age_c + m0 + (m6 + m12 + m18) + (m6 + m12 + m18):active, 
  data = contrastData0)) %>%
  mutate(female=0, m0=0, m6=0, m12=0, m18=0)

rownames(contrastMatrix_quad) <- rownames(contrastMatrix_cat) <- paste0('m', months[-1])

## ------------------------------------------------------------------------
summary(glht(cLDAquad, linfct = as.matrix(contrastMatrix_quad)))

## ------------------------------------------------------------------------
summary(glht(cLDAsymHet, linfct = as.matrix(contrastMatrix_cat)))

## ----results = 'hide', echo = FALSE--------------------------------------
# get default predictor matrix
ini_mi <- mice(trial_wide, maxit = 0, print = FALSE)
predictorMatrix <- ini_mi$predictorMatrix
# Don't want ADAS11.m12 predict by ADAS11.m18, etc.:
predictorMatrix['ADAS11.m12', 'ADAS11.m18'] <- 0
predictorMatrix['ADAS11.m6', 'ADAS11.m12'] <- 0
predictorMatrix['ADAS11.m6', 'ADAS11.m18'] <- 0

## ----results = 'hide', echo = TRUE---------------------------------------
trial_imp <- mice(trial_wide, predictorMatrix=predictorMatrix, seed = 20170714, maxit=100)

## ----echo = TRUE, size = 'scriptsize'------------------------------------
head(trial_wide)    # raw data with missing values:
head(complete(trial_imp)) # first complete version:

## ----echo = FALSE, size = 'scriptsize'-----------------------------------
fits_mi <- with(data=trial_imp, lm(ADAS11.m18~active*center(ADAS11.m0)))
summary(fits_mi)

## ----echo = FALSE--------------------------------------------------------
printCoefmat(summary(pool(fits_mi))[,1:5],eps=0.001,digits=3)

## ------------------------------------------------------------------------
post <- trial_imp$post
k_tipping <- seq(1, 2.5, 0.25)
est_tipping <- vector("list", length(k_tipping))
for (k in 1:length(k_tipping)){
  # increase imputed ADAS11.m18 in the active group by 
  # factor of k x MAR treatment estimate (3.1)
  # (nullify the imputed treatment effect to varying degrees)
  post["ADAS11.m18"] <- paste("imp[[j]][,i] <- imp[[j]][,i] + ", k_tipping[k], "* 4.0 * p$data$active[!r[,j]]")
  imp_k <- mice(trial_wide, post=post, predictorMatrix=predictorMatrix, seed = 20170714, maxit=100, print = FALSE)
  fit_k <- with(imp_k, lm(ADAS11.m18~active*center(ADAS11.m0)))
  est_tipping[[k]] <- summary(pool(fit_k))['active', 1:5]
}

## ----echo=FALSE----------------------------------------------------------
results_tipping <- do.call(rbind, est_tipping)
results_tipping <- cbind(tipping_factor = k_tipping, results_tipping)
printCoefmat(results_tipping, eps=0.001, digits=3)

## ------------------------------------------------------------------------
imputed_ids <- subset(trial_wide, is.na(ADAS11.m18))[1:5,]$id
# filter(trial_wide, id %in% imputed_ids)

## ------------------------------------------------------------------------
# Imputation assuming under MAR 
filter(complete(trial_imp), id %in% imputed_ids) %>%
  select(id, female, age, group, ADAS11.m0, ADAS11.m6, ADAS11.m12, ADAS11.m18)

## ------------------------------------------------------------------------
# Imputation under MNAR (k=0.5)
filter(complete(imp_k), id %in% imputed_ids) %>%
  select(id, female, age, group, ADAS11.m0, ADAS11.m6, ADAS11.m12, ADAS11.m18)

