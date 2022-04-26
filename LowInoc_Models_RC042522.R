##-----------------------------------------------------------
## Low inoculation paper data analysis
## Project: NYSG R/SHH-18
## Date: Apr 25, 2022
##-----------------------------------------------------------
## Load packages

# Install essential packages if haven't already.
required_packages <- c('plyr', 'tidyverse', 'janitor', 'lme4', 'broom.mixed', 'emmeans', 'NCmisc')
for (p in required_packages) {
  if(!require(p,character.only = TRUE)) {
    install.packages(p, dep = TRUE)
  }
}
##-----------------------------------------------------------
setwd("~/OneDrive - Cornell University/Cornell document/Food safety lab/Salmon_Project/low_inoc_experiment/GitHub")
##-----------------------------------------------------------
## Load the dataset
li_df <- read.csv("LowInoc_Data_RC042522.csv") # 216
li_df$temperature <- as.factor(li_df$temperature)
li_df$day <- as.factor(li_df$day)
li_df$agegroup <- as.factor(li_df$agegroup)
li_df$nc <- as.factor(li_df$nc)
li_df$max_tmt <- as.factor(li_df$max_tmt)
##-----------------------------------------------------------
## Get the 250 ppm data
li_250 <- li_df %>% filter(max_tmt=="yes") # 54
li_250 <- li_250 %>% mutate(elimination_label = case_when(elimination == 0 ~ "presence", elimination == 1 ~ "absence"))
##-----------------------------------------------------------
## Chi-Squared test of 250 ppm data
tbl <- table(li_df$max_tmt, li_df$elimination)
chisq <- chisq.test(tbl)
# X-squared = 65.424, df = 1, p-value = 6.039e-16
# The complete elimination of L. monocytogenes is significantly associated with the 250 ppm nisin treatment
##-----------------------------------------------------------
# The Serotype Elimination Model
li_melo_sero <- glmer(elimination ~ serotype+ (1|agegroup), data = li_250, family = binomial(link = "logit"))
sero_or <- tidy(li_melo_sero, conf.int=TRUE, exponentiate=TRUE, effects="fixed", confint.method = "profile")

# The Temperature Elimination Model
li_melo_temp <- glmer(elimination ~ temperature + (1|agegroup), data = li_250, family = binomial(link = "logit"))
temp_or <- tidy(li_melo_temp,conf.int=TRUE,exponentiate=TRUE,effects="fixed", confint.method = "profile")

# The Day Elimination Model
li_melo_day <- glmer(elimination ~ day + (1|agegroup), data = li_250, family = binomial(link = "logit"))
day_or <- tidy(li_melo_day,conf.int=TRUE,exponentiate=TRUE,effects="fixed", confint.method = "profile")
##-----------------------------------------------------------
## Get the data of samples with detectable L. monocytogenes levels
li_df_tmd <- subset(li_df, elimination == 0)
##-----------------------------------------------------------
## The Level Model
li_lme_17 <- lmer(MPN_CVT ~ nc+temperature+serotype+day + nc:serotype + nc:day + temperature:serotype + temperature:day +
                    (1|agegroup), data = li_df_tmd, REML = TRUE)
# Type III ANOVA
anova(li_lme_17)

# Post hoc pairwise comparison - nc:serotype
li_nc_ser.emm <- emmeans(li_lme_17, ~nc+serotype)
li_nc_ser.emm_df <- as.data.frame(li_nc_ser.emm)
li_nc_ser.ctr1 <- contrast(li_nc_ser.emm, interaction = "pairwise", simple = "serotype", by = "nc", combine = TRUE, adjust = "tukey")

# Post hoc pairwise comparison - nc:day
li_nc_day.emm <- emmeans(li_lme_17, ~nc+day)
li_nc_day.emm_df <- as.data.frame(li_nc_day.emm)
li_nc_day.ctr1 <- contrast(li_nc_day.emm, interaction = "pairwise", simple = "nc", by = "day", combine = TRUE, adjust = "tukey")

# Post hoc pairwise comparison - temperature:serotype
li_tem_ser.emm <- emmeans(li_lme_17, ~serotype+temperature)
li_tem_ser.emm_df <- as.data.frame(li_tem_ser.emm)
li_tem_ser.ctr1 <- contrast(li_tem_ser.emm, interaction = "pairwise", simple = "temperature", by = "serotype", combine = TRUE, adjust = "tukey")
li_tem_ser.ctr2 <- contrast(li_tem_ser.emm, interaction = "pairwise", simple = "serotype", by = "temperature", combine = TRUE, adjust = "tukey")

# Post hoc pairwise comparison - temperature:day
li_tem_day.emm <- emmeans(li_lme_17, ~temperature+day)
li_tem_day.emm_df <- as.data.frame(li_tem_day.emm)
li_tem_day.ctr1 <- contrast(li_tem_day.emm, interaction = "pairwise", simple = "temperature", by = "day", combine = TRUE, adjust = "tukey")
##-----------------------------------------------------------


