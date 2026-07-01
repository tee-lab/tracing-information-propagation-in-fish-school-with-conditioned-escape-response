## Code to reproduce fig S2

library(lme4)
library(emmeans)
library(ggplot2)
library(dplyr)
library(ggeffects)

## calculate mean speeds for a given cr for initial phase

df.mean.spd.cr.ini <- read.csv('mean_spd_cr_ini_video_data.csv', header = F)
col.names <- c("video", "f.id", "mean.spd.cr.ini")
colnames(df.mean.spd.cr.ini) <- col.names

df.mean.spd.cr.ini$video <- factor(df.mean.spd.cr.ini$video)
df.mean.spd.cr.ini$f.id <- factor(df.mean.spd.cr.ini$f.id, levels = c("1", "2", "3", "4", "5"))

lme.mean.spd.cr.ini <- lmer(mean.spd.cr.ini ~ f.id + (1|video), data = df.mean.spd.cr.ini)
lme.mean.spd.cr.ini.null <- lmer(mean.spd.cr.ini ~ 1 + (1|video), data = df.mean.spd.cr.ini)
anova(lme.mean.spd.cr.ini, lme.mean.spd.cr.ini.null, test = "F")

## calculate mean speeds for a given cr for a relaxation phase.

df.mean.spd.cr.relax <- read.csv('mean_spd_cr_relax_video_data.csv', header = F)
col.names <- c("video", "f.id", "mean.spd.cr.relax")
colnames(df.mean.spd.cr.relax) <- col.names

df.mean.spd.cr.relax$video <- factor(df.mean.spd.cr.relax$video)
df.mean.spd.cr.relax$f.id <- factor(df.mean.spd.cr.relax$f.id, levels = c("1", "2", "3", "4", "5"))

lme.mean.spd.cr.relax <- lmer(mean.spd.cr.relax ~ f.id + (1|video), data = df.mean.spd.cr.relax)
lme.mean.spd.cr.relax.null <- lmer(mean.spd.cr.relax ~ 1 + (1|video), data = df.mean.spd.cr.relax)
anova(lme.mean.spd.cr.relax, lme.mean.spd.cr.relax.null, test = "F")

emm.spd.cr.relax <- emmeans(lme.mean.spd.cr.relax, ~ f.id)
f.id.comp.spd.cr.relax <- pairs(emm.spd.cr.relax, adjust = "tukey")

df.emm.spd.cr.relax <- as.data.frame(emm.spd.cr.relax)
df.comp.spd.cr.relax <- as.data.frame(f.id.comp.spd.cr.relax)
