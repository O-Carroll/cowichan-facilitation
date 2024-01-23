## O. Carroll. 2023. Interplay of plant competition and facilitation. ##
########################################################################
# LOAD PACKAGES: ####

# install.packages("")
library(tidyverse) 
library(readxl)
library(lme4)
library(lmerTest)
library(broom)
library(MuMIn)
library(emmeans)
library(multcomp)
library(cowplot)
library(GGally)
library(MASS)
library(bbmle)

theme_set(theme_bw())
options(scipen = 99)
select <- dplyr::select

# LOAD DATA: ####

# Main data frame:
data <- read_excel("2022_Final_Data/2022_Carroll_Cowichan.xlsx", na = "NA")

# Subset with only seed-producing phytometers (and removal of 1 influential outlier):
phytoms <- data |> filter(foc.seed.stndrd > 0, 
                          foc.seed.stndrd < 1000)

########################################################################
## ANALYSIS: ## ####
# 1: SUMMARISE SPECIES FECUNDITY, SUCCESS, FAILURE ----

# How many phytometers died, survived, and reproduced out of 480 total?
table(data$foc.anpp.stndrd > 0) 
# 97 completely died/gone, 383 had some biomass present by end of experiment
table(data$seed.logistic) 
# 56 survivors produced seed, 327 survivors failed to produce seed

# What was distribution of seed-producing individuals across each species?
table(phytoms$foc.id) 
# no F produced seed across whole experiment, B was best, C okay, L poor

# What was distribution of seed-producing individuals across each competition type?
table(phytoms$comp.type)

# What was distribution of seed-producing individuals across each sown.neighbor density?
table(phytoms$nei.density) 
boxplot(phytoms$foc.seed.stndrd ~ phytoms$nei.density)
table(data$nei.density) # i.e. of 96 no nei. phytoms., only 1 indiv produced seed
phytoms |> filter(nei.density == 0) |> select(foc.id) # Indiv. was Clarkia

# What was distribution of seed-producing individuals across each nut trt?
table(phytoms$nutrients)
boxplot(phytoms$foc.seed.stndrd ~ phytoms$nei.density)

# What was distribution of seed-producing individuals under each neighbor identity?
table(phytoms$nei.id) # C and B offered best protection
boxplot(phytoms$foc.seed.stndrd ~ phytoms$nei.id)
table(data$nei.id)

# 2: MODEL BERNOULLI REPRODUCTION/FAILURE PROCESS ----

#||||||||||||##
## OVERVIEW: ##

# Here, we use binomial GLMMs to test processes affecting the probability of seed production/failure. We assess 3 separate models corresponding to our 3 hypotheses (competition, facilitation, and environment). We include wholeplot.id as a random effect to account for the error structure introduced by our splitplot design. We compare models with independent and interaction effects of explanatory variables and evaluate which combination of parameters best explains variability in the probability of success/failure

## ||||||||||||||||||||||||||||||||||||||||##
# Standardize predictor variables (z scores):

data.mod <- data |>
  mutate(z.herbiv = 
           (nei.herbivore.destruction - mean(nei.herbivore.destruction, na.rm = T))/
           sd(nei.herbivore.destruction, na.rm = T),
         z.par = 
           (pct.par.interception - mean(pct.par.interception, na.rm = T))/
           sd(pct.par.interception, na.rm = T),
         z.water = 
           (water.use.index - mean(water.use.index, na.rm = T))/
           sd(water.use.index, na.rm = T),
         sown.nei.density = as.factor(sown.nei.density))

##||||||||||||||||||||||||||||||||||||||||||||##
## Build models specific to 3 core hypotheses ##

## 1 - Nutrients*Neighbor density:

# Build model:
mod.nei.dens <- glmer(seed.logistic ~ nutrients * as.numeric(nei.density) + (1|wholeplot.id), 
                      data    = (data.mod |> 
                                   mutate(nutrients = factor(nutrients, 
                                                             levels = c("low","high")))), 
                      family  = "binomial", 
                      control = glmerControl(optimizer = "bobyqa", 
                                             optCtrl   = list(maxfun = 100000)),
                      nAGQ    = 10)
            
summary(mod.nei.dens) 
r.squaredGLMM(mod.nei.dens)

mod.nei.dens.2 <- update(mod.nei.dens, 
                         seed.logistic ~ as.numeric(nei.density) + (1|wholeplot.id))

summary(mod.nei.dens.2)
r.squaredGLMM(mod.nei.dens.2)

mod.nei.dens.3 <- update(mod.nei.dens, 
                         seed.logistic ~ nutrients + as.numeric(nei.density) + (1|wholeplot.id))

summary(mod.nei.dens.3)
r.squaredGLMM(mod.nei.dens.3)

# Compare:
AIC(mod.nei.dens, mod.nei.dens.2, mod.nei.dens.3)
anova(mod.nei.dens, mod.nei.dens.2, mod.nei.dens.3)

# Inspect mod.nei.dens.3:
anova(mod.nei.dens.3)
car::Anova(mod.nei.dens.3)

# Extract predictions from mod.nei.dens.3:
preds.mod.nei.dens <- data.frame(wholeplot.id = rep(levels(as.factor(data.mod$wholeplot.id)), 
                                           each = 20),
                                 nutrients   = rep(levels(as.factor(data.mod$nutrients)), 
                                                   each = 10),
                                 nei.density = rep(seq(from = 0, to = 2, length.out = 10)))

preds.mod.nei.dens.2 <- data.frame(preds = predict(mod.nei.dens.3, 
                                            newdata = preds.mod.nei.dens, 
                                            type    = "response")) |>
  cbind(preds.mod.nei.dens) |>
  group_by(nutrients, nei.density) |>
  summarise(preds = mean(preds))

## 2 - Nutrients * Herbivore damage:

# Build models:
mod.herbiv <- glmer(seed.logistic ~ nutrients * nei.herbivore.destruction + (1|wholeplot.id), 
                    data    = (data.mod |> 
                                 mutate(nutrients = factor(nutrients, levels = c("low","high")))), 
                    family  = "binomial", 
                    control = glmerControl(optimizer = "bobyqa", 
                                           optCtrl   = list(maxfun = 100000)),
                    nAGQ    = 10)

summary(mod.herbiv) 
r.squaredGLMM(mod.herbiv)

mod.herbiv.2 <- update(mod.herbiv,
                       seed.logistic ~ nei.herbivore.destruction + (1|wholeplot.id))

summary(mod.herbiv.2)
r.squaredGLMM(mod.herbiv.2)

mod.herbiv.3 <- update(mod.herbiv,
                       seed.logistic ~ nutrients + nei.herbivore.destruction + (1|wholeplot.id))

summary(mod.herbiv.3)
r.squaredGLMM(mod.herbiv.3)

# Compare:
AIC(mod.herbiv, mod.herbiv.2, mod.herbiv.3)
anova(mod.herbiv, mod.herbiv.2, mod.herbiv.3)

# Inspect mod.herbiv.3:
anova(mod.herbiv.3)
car::Anova(mod.herbiv.3, type = "III")

# Extract preds:
preds.mod.herbiv <- data.frame(wholeplot.id = rep(levels(as.factor(data.mod$wholeplot.id)), 
                                                  each = 20),
                               nutrients    = rep(levels(as.factor(data.mod$nutrients)), 
                                                  each = 10),
                               nei.herbivore.destruction = rep(seq(from = 1, 
                                                                   to = 100, 
                                                                   length.out = 10)))

preds.mod.herbiv.2 <- data.frame(preds = predict(mod.herbiv.3, 
                                                 newdata = preds.mod.herbiv, 
                                                 type    = "response")) |>
  cbind(preds.mod.herbiv) |>
  group_by(nutrients, nei.herbivore.destruction) |>
  summarise(preds = mean(preds))

## 3 - Nutrients * Resources (water and light):

# Build models:

mod.res <- glmer(seed.logistic ~ (nutrients + z.par + z.water)^2 + (1|wholeplot.id), 
                 data    = (data.mod |> 
                              mutate(nutrients = factor(nutrients, levels = c("low","high")))), 
                 family  = "binomial", 
                 control = glmerControl(optimizer = "bobyqa", 
                                        optCtrl   = list(maxfun = 100000)),
                 nAGQ    = 10)

summary(mod.res)
r.squaredGLMM(mod.res)

mod.res.2 <- update(mod.res,
                    seed.logistic ~ nutrients + z.par + z.water + (1|wholeplot.id))

mod.res.3 <- update(mod.res,
                    seed.logistic ~ nutrients + z.par + (1|wholeplot.id))

mod.res.4 <- update(mod.res,
                    seed.logistic ~ nutrients + z.water + (1|wholeplot.id))

mod.res.5 <- update(mod.res,
                    seed.logistic ~ z.par + z.water + (1|wholeplot.id))

mod.res.6 <- update(mod.res,
                    seed.logistic ~ nutrients + (z.par + z.water)^2 + (1|wholeplot.id))

# Compare
AIC(  mod.res, mod.res.2, mod.res.3, mod.res.4, mod.res.5, mod.res.6)
anova(mod.res, mod.res.2, mod.res.3, mod.res.4, mod.res.5, mod.res.6)

# Inspect mod.res.6:
anova(mod.res.6)
car::Anova(mod.res.6, type = "III")

# Address singularity in raneff structure - 
# Test inclusion of wholeplot as fixeff, and compare with model that removes wholeplot altogether
mod.fix <- glm(seed.logistic ~ nutrients + wholeplot.id + (z.par + z.water)^2,
               data    = (data.mod |> 
                            mutate(nutrients = factor(nutrients, levels = c("low","high")))),
               family  = "binomial",
               control = glm.control(maxit = 100))

mod.fix.2 <- glm(seed.logistic ~ nutrients + (z.par + z.water)^2, # Dropping wholeplot
                 data    = (data.mod |> 
                              mutate(nutrients = factor(nutrients, levels = c("low","high")))), 
                 family  = "binomial",
                 control = glm.control(maxit = 100))

AIC(mod.fix, mod.fix.2) # Removing it is improvement
anova(mod.fix, mod.fix.2)

# Coeffs from fixed effect model and random effect model (with singularity) are identical:
summary(mod.fix)
summary(mod.fix.2)
summary(mod.res.6)

# Will proceed with fixed effect only model
summary(mod.fix.2)
r.squaredGLMM(mod.fix.2)
anova(mod.fix.2)
round(coef(summary(mod.fix.2)), 4)

# Extract predictions:
# PAR at different constant values of water (mean, and plus/minus 1 sd):
preds.res <- data.frame(nutrients = rep(levels(as.factor(data.mod$nutrients)), 
                                                   each = 100),
                                   z.par     = rep(seq(from       = min(data.mod$z.par),
                                                       to         = max(data.mod$z.par), 
                                                       length.out = 100), 
                                                   times = 2))

preds.res.2 <- data.frame(water.at.mean =
                                        predict(mod.fix.2, 
                                                newdata = (preds.res |> 
                                                             mutate(z.water = 0)),
                                                type    = "response"),
                                      water.minus.1.sd = 
                                        predict(mod.fix.2, 
                                                newdata = (preds.res |> 
                                                             mutate(z.water = -1)), 
                                                type    = "response"),
                                      water.plus.1.sd = 
                                        predict(mod.fix.2, 
                                                newdata = (preds.res |> 
                                                             mutate(z.water = +1)), 
                                                type = "response")) |>
  cbind(preds.res) |>
  mutate(backtransf.par = 
           sd(data.mod$pct.par.interception) * z.par + mean(data.mod$pct.par.interception)) |>
  gather(water.level, pred, -c(nutrients, z.par, backtransf.par))

## Clear unused models
rm(mod.nei.dens, mod.nei.dens.2,
   mod.herbiv, mod.herbiv.2,
   mod.res, mod.res.2, mod.res.3, mod.res.4, mod.res.5, mod.res.6,
   mod.fix)

## Extract predicted probabilities over each predictor for plotting:
# PAR at different constant values of water (mean, and plus/minus 1 sd):
water.par.interac.df <- data.frame(nutrients = rep(levels(as.factor(data.mod$nutrients)), 
                                                   each = 100),
                                   z.par     = rep(seq(from       = min(data.mod$z.par),
                                                       to         = max(data.mod$z.par), 
                                                       length.out = 100), times = 2))

water.par.interac.preds <- data.frame(water.at.mean = 
                                        predict(mod.fix.2, 
                                                newdata = (water.par.interac.df |> 
                                                             mutate(z.water = 0)),
                                                type    = "response"),
                                      water.minus.1.sd = 
                                        predict(mod.fix.2, 
                                                newdata = (water.par.interac.df |> 
                                                             mutate(z.water = -1)), 
                                                type    = "response"),
                                      water.plus.1.sd = 
                                        predict(mod.fix.2, 
                                                newdata = (water.par.interac.df |> 
                                                             mutate(z.water = +1)), 
                                                type = "response")) |>
  cbind(water.par.interac.df) |>
  mutate(backtransf.par = 
           sd(data.mod$pct.par.interception) * z.par + mean(data.mod$pct.par.interception)) |>
  gather(water.level, pred, -c(nutrients, z.par, backtransf.par))

# 3: MODEL NEG. BINOM. SEED PRODUCTION PROCESS ----

# Assess number of obs per wholeplot group:
table(phytoms$wholeplot.id)
# Can't include wholeplot as random effect as 3 groups only have 1 successful obs:
# Including as fixed eff instead.

# Dataframe with standardized continuous predictors:
phytom.mod.df <- phytoms |>
  mutate(z.par = (pct.par.interception - mean(pct.par.interception, na.rm = T))/
           sd(pct.par.interception, na.rm = T),
         z.water = (water.use.index - mean(water.use.index, na.rm = T))/
           sd(water.use.index, na.rm = T),
         sown.nei.density = as.numeric(sown.nei.density)) |>
  filter(id != 95) |>  # Identified as influential outlier in earlier mod diagnostics
  filter(foc.id != "l")

## Build model for neighbor resource use: ##
# Max model:
max.nb <- glm.nb(foc.seed.stndrd ~ nutrients + foc.id + (z.par + z.water)^2 + wholeplot.id, 
                 data = (phytom.mod.df |> 
                           mutate(nutrients = factor(nutrients, levels = c("low","high")))),
                 link = log,
                 control = glm.control(maxit = 100)) 

summary(max.nb)

# Assess influence of modifications:
max.nb.1 <- update(max.nb, foc.seed.stndrd ~ nutrients + foc.id + z.par + z.water + wholeplot.id)
max.nb.2 <- update(max.nb, foc.seed.stndrd ~ nutrients + foc.id + z.par + z.water)
max.nb.3 <- update(max.nb, foc.seed.stndrd ~ nutrients + foc.id + z.par + z.water)
max.nb.4 <- update(max.nb, foc.seed.stndrd ~ foc.id + z.par + z.water)
max.nb.5 <- update(max.nb, foc.seed.stndrd ~ foc.id)
max.nb.6 <- update(max.nb, foc.seed.stndrd ~ nutrients + foc.id + z.par)
max.nb.7 <- update(max.nb, foc.seed.stndrd ~ nutrients + foc.id)
max.nb.8 <- update(max.nb, foc.seed.stndrd ~ nutrients + foc.id  + z.par)
max.nb.9 <- update(max.nb, foc.seed.stndrd ~ foc.id + z.par)
max.nb.10 <- update(max.nb, foc.seed.stndrd ~ nutrients + z.par)
max.nb.11 <- update(max.nb, foc.seed.stndrd ~ (nutrients + foc.id + z.par)^2)
max.nb.12 <- update(max.nb, foc.seed.stndrd ~ 1)
max.nb.13 <- update(max.nb, foc.seed.stndrd ~ nutrients + foc.id + (z.par)^2)
max.nb.14 <- update(max.nb, foc.seed.stndrd ~ foc.id + z.par + wholeplot.id)
max.nb.15 <- update(max.nb, foc.seed.stndrd ~ 
                      nutrients + foc.id + z.par + nutrients:foc.id + foc.id:z.par)
max.nb.16 <- update(max.nb, foc.seed.stndrd ~ 
                      nutrients + (foc.id + z.par)^2)
max.nb.17 <- update(max.nb, foc.seed.stndrd ~ 
                      nutrients + foc.id + z.par + nutrients:foc.id + foc.id:z.par)

# Compare and select:
AIC(max.nb, max.nb.1, max.nb.2, max.nb.3, max.nb.4, max.nb.5, 
    max.nb.6, max.nb.7, max.nb.8, max.nb.9, max.nb.10,
    max.nb.11, max.nb.12, max.nb.13, max.nb.14, max.nb.15, 
    max.nb.16, max.nb.17)

# Evaluate:
anova(max.nb.15) # Chi squared test based on deviance
summary(max.nb.15) # Wald test
car::Anova(max.nb.15, type = "III")
r.squaredGLMM(max.nb.15)
confint(max.nb.15)
coef(max.nb.15)

# Extract predictions:
preds.nb <- data.frame(nutrients = rep(levels(as.factor(phytom.mod.df$nutrients)), 
                                       each = 100),
                       foc.id = "b",
                       z.par = rep(seq(from = phytom.mod.df |> 
                                         filter(foc.id == 'b') |> 
                                         select(z.par) |> 
                                         min(), 
                                       to = phytom.mod.df |> 
                                         filter(foc.id == 'b') |> 
                                         select(z.par) |> 
                                         max(), 
                                       length.out = 100), 
                                   times = 2)) |>
  rbind(data.frame(nutrients = rep(levels(as.factor(phytom.mod.df$nutrients)), 
                                   each = 100),
                   foc.id = "c",
                   z.par = rep(seq(from = phytom.mod.df |> 
                                     filter(foc.id == 'c') |> 
                                     select(z.par) |> 
                                     min(), 
                                   to = phytom.mod.df |> 
                                     filter(foc.id == 'c') |> 
                                     select(z.par) |> 
                                     max(), 
                                   length.out = 100), 
                               times = 2))) |>
  mutate(backtransf.par = 
           sd(phytom.mod.df$pct.par.interception) * z.par + 
           mean(phytom.mod.df$pct.par.interception))

preds.nb.2 <- preds.nb |>
  cbind(preds = exp(predict(max.nb.15, newdata = preds.nb)))

## Build models for response to planted neighbor density: ##
mod3a.1 <- glm.nb(foc.seed.stndrd ~ nutrients + foc.id + sown.nei.density,
                  data    = (phytom.mod.df |> 
                               mutate(nutrients = factor(nutrients, levels = c("low","high")))),
                  link    = log,
                  control = glm.control(maxit = 100)) 

mod3a.2 <- update(mod3a.1, foc.seed.stndrd ~ nutrients + sown.nei.density)
mod3a.3 <- update(mod3a.1, foc.seed.stndrd ~ nutrients)
mod3a.4 <- update(mod3a.1, foc.seed.stndrd ~ sown.nei.density)

# Select:
anova(mod3a.1,mod3a.2,mod3a.3,mod3a.4)
AIC(mod3a.1,mod3a.2,mod3a.3,mod3a.4)

# Evaluate
summary(mod3a.1)
car::Anova(mod3a.1, type = "III")
anova(mod3a.1)
r.squaredGLMM(mod3a.1)
# Not shown in figure because not significant, but shown in supp table 5a

## Build models for response to herbivore destruction of plot: ##
mod3b.1 <- glm.nb(foc.seed.stndrd ~ (nutrients + foc.id + nei.herbivore.destruction)^2,
                  data    = (phytom.mod.df |> 
                               mutate(nutrients = factor(nutrients, levels = c("low","high")))),
                  link    = log,
                  control = glm.control(maxit = 100)) 

mod3b.2 <- update(mod3b.1, foc.seed.stndrd ~ nutrients + foc.id + nei.herbivore.destruction)
mod3b.3 <- update(mod3b.1, foc.seed.stndrd ~ nutrients)
mod3b.4 <- update(mod3b.1, foc.seed.stndrd ~ nei.herbivore.destruction)
mod3b.5 <- update(mod3b.1, foc.seed.stndrd ~ foc.id)
mod3b.6 <- update(mod3b.1, foc.seed.stndrd ~ 
                    nutrients + foc.id + nei.herbivore.destruction + nutrients:foc.id)

# Select:
anova(mod3b.1,mod3b.2,mod3b.3,mod3b.4,mod3b.5,mod3b.6)
AIC(mod3b.1,mod3b.2,mod3b.3,mod3b.4,mod3b.5,mod3b.6)
AIC(mod3b.1,mod3b.2,mod3b.6)

# Evaluate:
summary(mod3b.6)
car::Anova(mod3b.6, type = "III")
anova(mod3b.6)
r.squaredGLMM(mod3b.6)

# Extract predictions for figures:
preds.fig3b <- data.frame(nutrients = rep(levels(as.factor(phytom.mod.df$nutrients)), each = 100),
                          foc.id    = "b",
                          nei.herbivore.destruction = rep(seq(from = 0, 
                                                              to   = 100, 
                                                              length.out = 100), 
                                                          times = 2)) |>
  rbind(data.frame(nutrients = rep(levels(as.factor(phytom.mod.df$nutrients)), each = 100),
                   foc.id    = "c",
                   nei.herbivore.destruction = rep(seq(from = 0, 
                                                       to   = 100, 
                                                       length.out = 100), 
                                                   times = 2)))

preds.fig3b.2 <- preds.fig3b |>
  cbind(preds = exp(predict(mod3b.6, newdata = preds.fig3b)))

# 4: COMPETITION MODELS (BY NEIGHBOR ID) ----

lambda.data <- data |>
  select(id, wholeplot.id, nutrients, nei.id, nei.density, foc.id, foc.seed.stndrd) |>
  filter(nei.id == 'n') |>
  group_by(wholeplot.id, nutrients, foc.id) |>
  summarise(foc.seed.stndrd = round(mean(foc.seed.stndrd), 0)) |> 
  mutate(nei.density = as.factor(0))

# Each focal-neighbor pair gets a no-neighbor measure of fecundity (repeated so same for focal across all neighbors)
lambda.longer <- rbind(lambda.data, lambda.data, lambda.data, lambda.data) |>
  cbind(nei.id = rep(c('b','c','f','l'), each = 48)) |>
  full_join(data |> 
              select(wholeplot.id, nutrients, foc.id, foc.seed.stndrd, nei.density, nei.id) |>
              filter(nei.id != 'n') ) |>
  mutate(sown.nei.density = ifelse(nei.density == 0, 0, ifelse(nei.density == 1, 2, 10)))  |>
  group_by(foc.id, nei.id, nutrients) |> 
  mutate(seed.logistic = ifelse(foc.seed.stndrd == 0, 0, 1))

comp.nested <- lambda.longer |>
  group_by(nutrients, foc.id, nei.id) |>
  nest()

# Fit negbinom model with starts as mean:

nbinom.nei <- comp.nested |>
  mutate(coeffs = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA,
                          purrr::map(data, ~ mle2(foc.seed.stndrd ~ dnbinom(mu   = mu, 
                                                                            size = size),   
                  # 'map' applies mle2 param estimation function to each element of the nested list
                                                  start = list(mu   = 
                                                                 mean(data[[1]]$foc.seed.stndrd), 
                                                               size = 1), 
                  # gives mean within levels of list (i.e. for each focal sp.)
                                                  data = .x, 
                                                  control = list(maxit = 300000), 
                                                  method = "Nelder-Mead")) ),
         convergence = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA,
                               coeffs[[1]]@details[["convergence"]] ), 
         aic         = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                               AIC(coeffs[[1]]) ),
         aicc        = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                               AICc(coeffs[[1]]) ),
         mu          = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                               coeffs[[1]]@coef[["mu"]] ),
         size        = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                               coeffs[[1]]@coef[["size"]] ) )

# Fit comp model using preds as starts:

comp.nei <- nbinom.nei |>
  select(nutrients:coeffs) |>
  mutate(comp.coeffs = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA,
                               purrr::map(data, ~ mle2(foc.seed.stndrd ~ 
                                                         dnbinom(mu = 
                                                                   lambda/(1 + a12*sown.nei.density),
                                                                 size = size),
                                                       start = list(lambda =
                                                                      coeffs[[1]]@coef[["mu"]], 
                                                                    a12  = 0,
                                                                    size =
                                                                      coeffs[[1]]@coef[["size"]]), 
                                                       data = .x, 
                                                       control = list(maxit = 300000),
                                                       method = "Nelder-Mead")) ),
         convergence = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                               comp.coeffs[[1]]@details[["convergence"]] ), 
         aic = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                       AIC(comp.coeffs[[1]]) ),
         aicc = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                        AICc(comp.coeffs[[1]]) ),
         lambda = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                          comp.coeffs[[1]]@coef[["lambda"]] ),
         alpha = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                         comp.coeffs[[1]]@coef[["a12"]] ),
         se = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                      stdEr(comp.coeffs[[1]]) ) )

comp.coeffs.nei <- comp.nei |>
  select(nutrients, foc.id, nei.id, lambda, alpha, se) |>
  mutate(lambda = round(lambda, 0),
         alpha  = round(alpha, 2),
         se     = round(se, 3)) |> 
# Manually setting coeffs. highlighted by reviewer to NA, Jan 2024:
  mutate(across(.cols = c("lambda", "alpha", "se"),
                .fns  = ~ case_when(
                  nutrients == "low" & foc.id == "b" & nei.id == "f" ~ NA,
                  nutrients == "low" & foc.id == "l" & nei.id == "c" ~ NA,
                  .default = .x)))

# 5: COMPETITION MODELS (IGNORING NEIGHBOR ID) ----

# Nest data:
test.comp.nonei <- data |>
  group_by(nutrients, foc.id) |>
  mutate(seed.logistic = ifelse(foc.seed.stndrd == 0, 0, 1)) |> 
  nest()

# Fit negbinom model with starts as mean:
nbinom.nonei <- test.comp.nonei |>
  mutate(coeffs = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA,
                          purrr::map(data, ~ mle2(foc.seed.stndrd ~ dnbinom(mu = mu, size = size),
                                                  start = list(mu = mean(data[[1]]$foc.seed.stndrd), 
                                                               size = 1),
                                                  data = .x, 
                                                  control = list(maxit = 300000),
                                                  method = "Nelder-Mead")) ),
         convergence = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA,
                               coeffs[[1]]@details[["convergence"]] ), 
         aic = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, AIC(coeffs[[1]]) ),
         aicc = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, AICc(coeffs[[1]]) ),
         mu = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, coeffs[[1]]@coef[["mu"]] ),
         size = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, coeffs[[1]]@coef[["size"]] ) )

# Fit comp model using preds as starts:
comp.nonei <- nbinom.nonei |>
  select(nutrients:coeffs) |>
  mutate(comp.coeffs = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA,
                               purrr::map(data, ~ mle2(foc.seed.stndrd ~ 
                                                         dnbinom(mu = lambda/(1 + a12*sown.nei.density),
                                                                 size = size),
                                                       start = list(lambda = coeffs[[1]]@coef[["mu"]], 
                                                                    a12 = 0,
                                                                    size = coeffs[[1]]@coef[["size"]]), 
                                                       data = .x, 
                                                       control = list(maxit = 300000),
                                                       method = "Nelder-Mead")) ),
         convergence = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                               comp.coeffs[[1]]@details[["convergence"]] ), 
         aic = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                       AIC(comp.coeffs[[1]]) ),
         aicc = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA,
                        AICc(comp.coeffs[[1]]) ),
         lambda = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                          comp.coeffs[[1]]@coef[["lambda"]] ),
         alpha = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                         comp.coeffs[[1]]@coef[["a12"]] ),
         se = ifelse( sum(data[[1]][["seed.logistic"]]) < 3, NA, 
                      stdEr(comp.coeffs[[1]]) ) )

comp.coeffs.nonei <- comp.nonei |>
  select(nutrients, foc.id, lambda, alpha) |>
  mutate(lambda = round(lambda, 0),
         alpha = round(as.numeric(alpha), 2),
         nut.foc = paste(foc.id, nutrients, sep = ".")) |> 
# Manually setting coeffs. highlighted by reviewer to NA, Jan 2024:
  mutate(across(.cols = c("lambda", "alpha"),
                .fns  = ~ case_when(
                  nutrients == "low" & foc.id == "l" ~ NA,
                  .default = .x)))

# 6: COMPARE NEIGHBOR TRAITS ----

## CANOPY HEIGHT ##
# Build model:
nei.canopy <- lmer(log(av.nei.canopy.height + 1) ~ nei.id*nutrients + (1|wholeplot.id), 
                   data = (data |> filter(nei.id!= 'n')) )

# Evaluate:
summary(nei.canopy)
plot(nei.canopy)
anova(nei.canopy)  # output
r.squaredGLMM(nei.canopy)
plot(emmeans(nei.canopy, "nei.id"), comparisons = TRUE)
emmeans(nei.canopy, list(pairwise ~ nei.id:nutrients), adjust = "tukey")
pwpm(emmeans(nei.canopy, ~ nei.id:nutrients))
cld(emmeans(nei.canopy, ~ nei.id:nutrients), Letters = letters)
round(coef(summary(nei.canopy)),3)

# Summarise data for output:
nei.canopy.summary <- data |>
  filter(nei.id != 'n') |>
  group_by(nutrients, nei.id) |>
  summarise(mean = mean(av.nei.canopy.height, na.rm = T),
            sd = sd(av.nei.canopy.height, na.rm = T),
            n = n(),
            se = sd/sqrt(n)) |> 
  merge(as.data.frame(cld(emmeans(nei.canopy, ~ nei.id:nutrients), Letters = letters)) |> 
          select(nutrients, nei.id, .group),
        by = c("nei.id", "nutrients")) |> 
  mutate(.group = str_remove(.group, " "),
         .group = str_remove(.group, " "),
         .group = str_remove(.group, " "))

## Water use ##
# Build model:
nei.water <- lmer(log(water.use.index + 1) ~ nei.id*nutrients + (1|wholeplot.id), 
                  data = (data |> filter(nei.id!= 'n')) )

# Evaluate:
summary(nei.water)
plot(nei.water)
anova(nei.water) # output
r.squaredGLMM(nei.water)
plot(emmeans(nei.water, "nei.id"), comparisons = TRUE)
pairs(emmeans(nei.water, "nei.id"))
emmeans(nei.water, list(pairwise ~ nei.id:nutrients), adjust = "tukey")
pwpm(emmeans(nei.water, ~ nei.id:nutrients))
cld(emmeans(nei.water, ~ nei.id:nutrients), Letters = letters)
round(coef(summary(nei.water)), 3)

# Summarise data for output
nei.water.summary <- data |>
  filter(nei.id != 'n') |>
  group_by(nutrients, nei.id) |>
  summarise(mean = mean(water.use.index, na.rm = T),
            sd = sd(water.use.index, na.rm = T),
            n = n(),
            se = sd/sqrt(n)) |> 
  merge(as.data.frame(cld(emmeans(nei.water, ~ nei.id:nutrients), Letters = letters)) |> 
          select(nutrients, nei.id, .group),
        by = c("nei.id", "nutrients")) |> 
  mutate(.group = str_remove(.group, "  "),
         .group = str_remove(.group, " "),
         .group = str_remove(.group, " "),
         .group = str_remove(.group, " "))

## PAR interception ##
# Build model:
nei.light <- lmer(log(pct.par.interception + 1) ~ nei.id*nutrients + (1|wholeplot.id), 
                  data = (data |> filter(nei.id!= 'n')) )

# Evaluate:
summary(nei.light)
plot(nei.light)
anova(nei.light) # output
r.squaredGLMM(nei.light)
plot(emmeans(nei.light, "nei.id"), comparisons = TRUE)
pairs(emmeans(nei.light, "nei.id"))
emmeans(nei.light, list(pairwise ~ nei.id:nutrients), adjust = "tukey")
pwpm(emmeans(nei.light, ~ nei.id:nutrients))
cld(emmeans(nei.light, ~ nei.id:nutrients), Letters = letters)
round(coef(summary(nei.light)),3)

# Summarise data for output:
nei.light.summary <- data |>
  filter(nei.id != 'n') |>
  group_by(nutrients, nei.id) |>
  summarise(mean = mean(pct.par.interception, na.rm = T),
            sd = sd(pct.par.interception, na.rm = T),
            n = n(),
            se = sd/sqrt(n)) |> 
  merge(as.data.frame(cld(emmeans(nei.light, ~ nei.id:nutrients), Letters = letters)) |> 
          select(nutrients, nei.id, .group),
        by = c("nei.id", "nutrients")) |> 
  mutate(.group = str_remove(.group, "  "),
         .group = str_remove(.group, " "),
         .group = str_remove(.group, " "))

## Herbivore destruction ##
# Build model:
nei.herbiv <- lmer(nei.herbivore.destruction ~ nei.id*nutrients + (1|wholeplot.id), 
                   data = (data |> filter(nei.id != 'n')) )

# Evaluate:
summary(nei.herbiv)
confint(nei.herbiv)
anova(nei.herbiv) # output
r.squaredGLMM(nei.herbiv)
plot(nei.herbiv)
pairs(emmeans(nei.herbiv, c("nei.id","nutrients")))
pwpm(emmeans(nei.herbiv, ~ nei.id:nutrients))
cld(emmeans(nei.herbiv, ~ nei.id:nutrients), Letters = letters)
coef(summary(nei.herbiv))

# Summarise data for figures/tables
nei.herbiv.summary <- data |>
  filter(nei.id!= 'n') |>
  group_by(nutrients, nei.id) |>
  summarise(mean = mean(nei.herbivore.destruction, na.rm = T),
            sd = sd(nei.herbivore.destruction, na.rm = T),
            n = n(),
            se = sd/sqrt(n)) |> 
  merge(as.data.frame(cld(emmeans(nei.herbiv, ~ nei.id:nutrients), Letters = letters)) |> 
          select(nutrients, nei.id, .group),
        by = c("nei.id", "nutrients")) |> 
  mutate(.group = str_remove(.group, " "),
         .group = str_remove(.group, " "))

########################################################################
## FIGURES  ## ----
# FIGURE 1 - Summary of outcomes  ----

# Barplot of survival by species
data.figs <- data |>
  mutate(planted = 1) |>
  select(id, foc.id, planted, anpp.logistic, seed.logistic) |>
  gather(stage, binary, -c(id, foc.id)) |>
  filter(binary != 0) 

fig.1 <- ggplot(data = data.figs, 
                 aes(x = stage)) + 
  geom_bar(aes(fill = foc.id), width = .8) +
  scale_x_discrete(limits = c("planted", "anpp.logistic", "seed.logistic"),
                   labels = c("Planted","Survived","Produced\nseed")) +
  ylab("Count") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,500), xlim = c(0.5,3.5), expand = FALSE) +
  scale_fill_manual(name = "", 
                    #values = c("b" = "forestgreen", "c" = "hotpink1", 
                    #          "f" = "dodgerblue3", "l" = "black"),
                    values = c("b" = "lightgray", "c" = "lightblue4", 
                               "f" = "lightsteelblue", "l" = "black"),
                    labels = c("b" = expression(italic("B. hordeaceus")), 
                               "c" = expression(italic("C. amoena")), 
                               "f" = expression(italic("F. myuros")),
                               "l" = expression(italic("L. multiflorum")))) +
  theme(legend.text.align = 0) + 
  theme(plot.margin = margin(.4,.4,.4,.4, "cm"))

legend.1 <- plot_grid( 
  get_legend( fig.1 + theme(legend.text = element_text(size  = 23, 
                                                       color = "black"),
                            legend.box.background = element_rect(colour = "black", 
                                                                 size   = 1))))
                                     

fig.1.l <- plot_grid(
  fig.1 + 
    theme(legend.position ="none") +
    draw_plot(legend.1, x = 2.45, y = 430, width = 1, height = 1))

save_plot("Figures/Figure_1.png", fig.1.l, base_width = 9, base_height = 9, dpi = 1000)

# FIGURE 2 - Probability of seed production ----

fig.2a <- ggplot(data = preds.mod.nei.dens.2,
                 aes(x = nei.density, y = preds, 
                     color = nutrients)) +
  geom_line(size = 3) +
  theme_bw() +
  xlab("Neighbor density level") + 
  ylab("Probability of seed production") +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,1), xlim = c(0, 2.1), expand = FALSE) +
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(labels = c(0, 1, 2), breaks = c(0, 1, 2)) +
  scale_color_manual(name = "Nutrient treatment", 
                     values = c("high" = "forestgreen", 
                                "low" = "orchid3"),
                     labels = c("high" = "NPK+", 
                                "low" = "Control")) +
  theme(plot.margin = margin(.4,.4,.4,.4, "cm"))

fig.2b <- ggplot(data  = preds.mod.herbiv.2,
                 aes(x = nei.herbivore.destruction, y = preds, 
                     color = nutrients)) +
  geom_line(size = 3) +
  theme_bw() +
  xlab("% Herbivore damage") + 
  ylab("") +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,1), xlim = c(0, 105), expand = FALSE) +
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(name = "Nutrient treatment", 
                     values = c("high" = "forestgreen", 
                                "low" = "orchid3"),
                     labels = c("high" = "NPK+", 
                                "low" = "Control")) +
  theme(plot.margin = margin(.4,.4,.4,.4, "cm"))

fig.2c <- ggplot(data = (preds.res.2 |>
                           group_by(water.level, backtransf.par) |>
                           summarise(pred = mean(pred)) |>
                           mutate(water.level = factor(water.level, levels = 
                                                         c("water.minus.1.sd",
                                                           "water.plus.1.sd",
                                                           "water.at.mean")))),
                 aes(x = backtransf.par, y = pred, 
                     color = water.level, linetype = water.level, alpha = water.level)) +
  geom_line(size = 3) +
  theme_bw() +
  xlab("% PAR interception") + 
  ylab("") +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,1), xlim = c(-5, 105), expand = FALSE) +
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(name = "Water use", 
                     values = c("water.at.mean" = "black",
                                "water.minus.1.sd" = "red",
                                "water.plus.1.sd" = "darkblue"),
                     labels = c("water.at.mean" = "Mean",
                                "water.minus.1.sd" = "Mean - 1 SD",
                                "water.plus.1.sd" = "Mean + 1 SD")) +
  scale_linetype_manual(name = "Water use", 
                        values = c("water.at.mean" = 1, 
                                   "water.minus.1.sd" = 5,
                                   "water.plus.1.sd" = 5),
                        labels = c("water.at.mean" = "Mean",
                                   "water.minus.1.sd" = "Mean - 1 SD",
                                   "water.plus.1.sd" = "Mean + 1 SD")) +
  scale_alpha_manual(name = "Water use", 
                     values = c("water.at.mean" = 1,
                                "water.minus.1.sd" = .5,
                                "water.plus.1.sd" = .5),
                     labels = c("water.at.mean" = "Mean",
                                "water.minus.1.sd" = "Mean - 1 SD",
                                "water.plus.1.sd" = "Mean + 1 SD")) +
  theme(plot.margin = margin(.4,.4,.4,.4, "cm"))

fig2a.legend <- plot_grid( get_legend( fig.2a +
                                         theme(legend.text=element_text(size=18, color = "black"),
                                               legend.title=element_text(size=18, color = "black"),
                                               legend.box.background = element_rect(colour = "black")))) 

fig2a.l <- plot_grid(fig.2a + theme(legend.position="none") + draw_plot(fig2a.legend, x = 0.2, y = 0.35),
                     labels = c("a"), label_size = 22,label_x = 0.08, label_y = 1)

fig2b.l <- plot_grid(fig.2b + theme(legend.position="none") + draw_plot(fig2a.legend, x = 33, y = 0.36), 
                     labels = c("b"), label_size = 22, label_x = 0.1, label_y = 1)

fig2c.legend <- plot_grid( get_legend( fig.2c +
                                         theme(legend.text=element_text(size=16, color = "black"),
                                               legend.title=element_text(size=16, color = "black"),
                                               legend.box.background = element_rect(colour = "black"))))
fig.2.2c <- plot_grid(
  fig.2c + theme(legend.position="none") + draw_plot(fig2c.legend, x = 32, y = 0.33), 
  labels = c("c"), label_size = 22, label_x = 0.1, label_y = 1)


fig2all <- plot_grid(fig2a.l, fig2b.l, fig.2.2c, ncol = 3)
fig2all

save_plot("Figures/Figure_2abc.png", fig2all, base_width = 16, base_height = 16/3, dpi = 1000)

# FIGURE 3 - % PAR interception (Neighbor species crowding) ----

fig3a <- ggplot(data = (preds.nb.2 |> filter(foc.id == 'b')),
                aes(x = backtransf.par, y = preds, color = nutrients)) +
  geom_line(size = 3) +
  theme_bw() +
  xlab("% PAR interception") + 
  ylab("Predicted seed production") +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,500), xlim = c(-5, 105), expand = FALSE) +
  scale_color_manual(name = "", 
                     values = c("high" = "forestgreen", 
                                "low" = "orchid3"),
                     labels = c("high" = "NPK+", 
                                "low" = "Control")) + 
  theme(plot.margin = margin(.4,.4,.4,.4, "cm"))

fig3b <- ggplot(data = (preds.nb.2 |> filter(foc.id == 'c')),
                aes(x = backtransf.par, y = preds, color = nutrients)) +
  geom_line(size = 3) +
  theme_bw() +
  xlab("% PAR interception") + 
  ylab("") +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,500), xlim = c(-5, 105), expand = FALSE) +
  scale_color_manual(name = "", 
                     values = c("high" = "forestgreen", 
                                "low" = "orchid3"),
                     labels = c("high" = "NPK+", 
                                "low" = "Control")) + 
  theme(plot.margin = margin(.4,.4,.4,.4, "cm"))

fig3.legend <- plot_grid( get_legend( fig3a +
                                        theme(legend.position = "bottom",
                                              legend.text=element_text(size=23, color = "black"),
                                              legend.box.background = element_rect(colour = "black")) +
                                        guides(colour = guide_legend(nrow = 2))) )

fig3a.l <- plot_grid(fig3a + theme(legend.position="none") + draw_plot(fig3.legend, x = 15, y = 450), 
                     labels = c("a"), label_size = 22)
fig3b.l <- plot_grid(fig3b + theme(legend.position="none"), labels = c("b"), label_size = 22)

fig3all <- plot_grid(fig3a.l, fig3b.l, ncol = 2)
fig3all

save_plot("Figures/Figure_3.png", fig3all, base_width = 11, base_height = 11/2, dpi = 1000)

# FIGURE 4 - % Herbivore damage ####

fig4a <- ggplot(data = (preds.fig3b.2  |> filter(foc.id == 'b')),
                   aes(x = nei.herbivore.destruction, y = preds, color = nutrients)) +
  geom_line(size = 3) +
  theme_bw() +
  xlab("% Herbivore damage") + 
  ylab("Predicted seed production") +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,500), xlim = c(-5, 105), expand = FALSE) +
  scale_color_manual(name = "", 
                     values = c("high" = "forestgreen", 
                                "low" = "orchid3"),
                     labels = c("high" = "NPK+", 
                                "low" = "Control")) + 
  theme(plot.margin = margin(.4,.4,.4,.4, "cm"))

fig4b <- ggplot(data = (preds.fig3b.2 |> filter(foc.id == 'c')),
                   aes(x = nei.herbivore.destruction, y = preds, color = nutrients)) +
  geom_line(size = 3) +
  theme_bw() +
  xlab("% Herbivore damage") + 
  ylab("") +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,500), xlim = c(-5, 105), expand = FALSE) +
  scale_color_manual(name = "", 
                     values = c("high" = "forestgreen", 
                                "low" = "orchid3"),
                     labels = c("high" = "NPK+", 
                                "low" = "Control")) + 
  theme(plot.margin = margin(.4,.4,.4,.4, "cm"))

fig4.legend <- plot_grid( get_legend( fig4a +
                                          theme(legend.position = "bottom",
                                                legend.text=element_text(size=23, color = "black"),
                                                legend.box.background = element_rect(colour = "black")) +
                                          guides(colour = guide_legend(nrow = 2))) )

fig4a.l <- plot_grid(fig4a + theme(legend.position="none") + draw_plot(fig4.legend, x = 15, y = 450), 
                        labels = c("a"), label_size = 22)
fig4b.l <- plot_grid(fig4b + theme(legend.position="none"), labels = c("b"), label_size = 22)

fig4.all <- plot_grid(fig4a.l, fig4b.l, ncol = 2)
fig4.all

save_plot("Figures/Figure_4.png", fig4.all, base_width = 11, base_height = 11/2, dpi = 1000)

# FIGURE 5 - Neighbor Resource Use ----

nei.height.fig <- ggplot(data = nei.canopy.summary,
                         aes(x = nei.id, y = mean, 
                             fill = nutrients, color = nutrients)) + 
  geom_bar(stat = "identity", position = "dodge", width = .9, lty = 1, color = "black") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), 
                stat     = "identity", 
                position = position_dodge(0.9), 
                width = .3, size = .5, color = "black") + # error bars show +/- 1 se
  theme_bw() +
  xlab("") + 
  ylab("Canopy height (mm)") +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,250), xlim = c(0.3,4.7), expand = FALSE) +
  scale_fill_manual(name = "", 
                    values = c("high" = "lightblue4", "low" = "gray"),
                    labels = c("high" = "NPK+", 
                               "low" = "Control")) + 
  theme(plot.margin = margin(.6,.4,.4,.4, "cm")) +
  scale_x_discrete(labels= c("B", "C", "F", "L")) +
  geom_text(aes(label = .group), 
            position  = position_dodge(0.9), 
            size   = 4, 
            vjust  = -4, 
            hjust  = .5, 
            colour = "black")

nei.water.fig <- ggplot(data = nei.water.summary,
                        aes(x = nei.id, y = mean, fill = nutrients, color = nutrients)) + 
  geom_bar(stat = "identity", position = "dodge", width = .9, lty = 1, color = "black") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), 
                stat     = "identity", 
                position = position_dodge(0.9), 
                width = .3, size = .5, color = "black") + # error bars show +/- 1 se
  theme_bw() +
  xlab("") + 
  ylab("Water use index") +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,1.5), xlim = c(0.3,4.7), expand = FALSE) +
  scale_fill_manual(name = "", 
                    values = c("high" = "lightblue4", "low" = "gray"),
                    labels = c("high" = "NPK+", 
                               "low" = "Control")) + 
  theme(plot.margin = margin(.6,.4,.4,.4, "cm")) +
  scale_x_discrete(labels= c("B", "C", "F", "L")) +
  geom_text(aes(label = .group), 
            position  = position_dodge(0.9), 
            size   = 4, 
            vjust  = -3.5, 
            hjust  = .5, 
            colour = "black")

nei.light.fig <- ggplot(data = nei.light.summary,
                        aes(x = nei.id, y = mean, 
                            fill = nutrients, color = nutrients)) + 
  geom_bar(stat = "identity", position = "dodge", width = .9, lty = 1, color = "black") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), 
                stat     = "identity", 
                position = position_dodge(0.9), 
                width = .3, size = .5, color = "black") + # error bars show +/- 1 se
  theme_bw() +
  xlab("") + 
  ylab("% PAR interception") +
  theme(axis.text.x=element_text(size=20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,60), xlim = c(0.3,4.7), expand = FALSE) +
  scale_fill_manual(name = "", 
                    values = c("high" = "lightblue4", "low" = "gray"),
                    labels = c("high" = "NPK+", 
                               "low" = "Control")) + 
  theme(plot.margin = margin(.6,.4,.4,.4, "cm")) +
  scale_x_discrete(labels= c("B", "C", "F", "L")) +
  geom_text(aes(label = .group), 
            position  = position_dodge(0.9), 
            size   = 4, 
            vjust  = -3, 
            hjust  = .5, 
            colour = "black")

nei.herbiv.fig <- ggplot(data = nei.herbiv.summary,
                         aes(x = nei.id, y = mean, 
                             fill = nutrients, color = nutrients)) + 
  geom_bar(stat = "identity", position = "dodge", width = .9, lty = 1, color = "black") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), 
                stat     = "identity", 
                position = position_dodge(0.9), 
                width    = .3, size = .5, color = "black") + # error bars show +/- 1 se
  theme_bw() +
  xlab("") + 
  ylab("% Herbivory") +
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y=element_text(size=18, color = "black"),
        axis.title=element_text(size=25, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  coord_cartesian(ylim = c(0,100), xlim = c(0.3,4.7), expand = FALSE) +
  scale_fill_manual(name = "", 
                    values = c("high" = "lightblue4", "low" = "gray"),
                    labels = c("high" = "NPK+", 
                               "low" = "Control")) + 
  theme(plot.margin = margin(.6,.4,.4,.4, "cm")) +
  scale_x_discrete(labels= c("B", "C", "F", "L")) +
  geom_text(aes(label = .group), 
            position  = position_dodge(0.9), 
            size   = 4, 
            vjust  = -3, 
            hjust  = .5, 
            colour = "black")

fig5.legend <- plot_grid( get_legend( nei.herbiv.fig +
                                        theme(legend.position = "bottom",
                                              legend.text=element_text(size=22, color = "black"),
                                              legend.box.background = element_rect(colour = "black")) +
                                        guides(colour = guide_legend(nrow = 2))) ) 

fig5a <- plot_grid(nei.height.fig + theme(legend.position="none"), labels = c("a"), label_size = 22)
fig5b <- plot_grid(nei.light.fig + theme(legend.position="none"), labels = c("b"), label_size = 22)
fig5c <- plot_grid(nei.water.fig + theme(legend.position="none"), labels = c("c"), label_size = 22)
fig5d <- plot_grid(nei.herbiv.fig + theme(legend.position="none") +
                     draw_plot(fig5.legend, x = 1.5, y = 90, width = 2, height = 2),
                   labels = c("d"), label_size = 22)

fig.5 <- plot_grid(fig5a, fig5b, fig5c, fig5d, nrow = 2, ncol = 2)
fig.5

save_plot("Figures/2023_Figure_5.png", fig.5, base_width = 10, base_height = 10, dpi = 1000)

########################################################################
## SUPPLEMENTARY INFORMATION ## ####
# Table S1 ####

# Not made using R 

# Table S2 ####

total.neinut <- data |>
  select(foc.id, nei.id, nutrients, foc.seed.stndrd, seed.logistic) |>
  group_by(foc.id, nei.id, nutrients) |>
  summarize(total.seed = sum(foc.seed.stndrd, na.rm = T),
            total.success = sum(seed.logistic, na.rm = T))

#writexl::write_xlsx(total.neinut, "FIGURES/total.neinut.xlsx")

total.neinut.density <- data |>
  select(foc.id, nei.id, nutrients, foc.seed.stndrd, seed.logistic, nei.density) |>
  group_by(foc.id, nei.id, nei.density, nutrients) |>
  summarize(total.seed = sum(foc.seed.stndrd, na.rm = T)) |>
  spread(nei.density, total.seed)

#writexl::write_xlsx(total.neinut.density, "FIGURES/total.neinut.density.xlsx")

comp.coeffs.nei 

# Table S3a ####

summary(mod.nei.dens.3)

# Table S3b ####

summary(mod.herbiv.3)

# Table S3c ####
bernoulli.glm <- data.frame(round(coef(summary(mod.fix.2)), 4)) |> 
  rownames_to_column("Parameter") # output

writexl::write_xlsx(bernoulli.glm, "FIGURES/step.one.glm.xlsx")

bernoulli.anovatable <- data.frame(round(car::Anova(mod.fix.2, type = "III"), 3)) |> 
  rownames_to_column("Factor")

writexl::write_xlsx(bernoulli.anovatable, "FIGURES/step.one.anovatable.xlsx")

# Table S4 ####

view(comp.coeffs.nonei)

# Table S5a ####

summary(mod3a.1)
car::Anova(mod3a.1, type = "III")

# Table S5b ####

neg.binomial.glm <- data.frame(round(coef(summary(max.nb.15)),4)) |> 
  rownames_to_column("Parameter") # output

#writexl::write_xlsx(neg.binomial.glm, "FIGURES/step.two.glm.xlsx")

neg.bin.glm.anovatable <- data.frame(round(car::Anova(max.nb.15, type = "III"), 3)) |> 
  rownames_to_column("Factor")

#writexl::write_xlsx(neg.bin.glm.anovatable, "FIGURES/step.two.anovatable.xlsx")

# Table S5c ####

summary(mod3b.6)

# Table S6a-d ####
supp.nei.canopy <- anova(nei.canopy) |> rownames_to_column("Parameter")
supp.nei.herbiv <- anova(nei.herbiv) |> rownames_to_column("Parameter") 
supp.nei.light <- anova(nei.light) |> rownames_to_column("Parameter")
supp.nei.water <- anova(nei.water) |> rownames_to_column("Parameter")

writexl::write_xlsx(supp.nei.canopy, "FIGURES/anova.nei.canopy.xlsx")
writexl::write_xlsx(supp.nei.herbiv, "FIGURES/anova.nei.herbiv.xlsx")
writexl::write_xlsx(supp.nei.light, "FIGURES/anova.nei.light.xlsx")
writexl::write_xlsx(supp.nei.water, "FIGURES/anova.nei.water.xlsx")

########################################################################