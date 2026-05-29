## Code to do statistical tests

start.time <- Sys.time()
library(lme4)
library(emmeans)
library(ggplot2)
library(dplyr)
library(ggeffects)

K <- 4
k <- 2
omega.ini <- 0.7
sight <- 2.36

fname_dcf <- paste0('lmer_dtcf_mual_0.1_muat_2.3_muesp_1.25_K_', K, '_k_', k, '_sight_', sight, '_gamma_0.25_omegaini_', omega.ini, '.csv')
df.dist.to.cf <- read.csv(fname_dcf, header = F)
col.names <- c("video", "f.id", "dcf")
colnames(df.dist.to.cf) <- col.names

df.dist.to.cf$video <- factor(df.dist.to.cf$video)
df.dist.to.cf$f.id <- factor(df.dist.to.cf$f.id, levels = c("2", "3", "4", "5"))

lme.dcf <- lmer(dcf ~ f.id + (1|video), data = df.dist.to.cf)
lme.dcf.null <- lmer(dcf ~ 1 + (1|video), data = df.dist.to.cf)
anova(lme.dcf, lme.dcf.null, test = 'F')

emm.dcf <- emmeans(lme.dcf, ~ f.id)
f.id.comp.dcf <- pairs(emm.dcf, adjust = "tukey")

df.emm.dfc <- as.data.frame(emm.dcf)
df.comp.dcf <- as.data.frame(f.id.comp.dcf)

plt.dfc <- ggplot(df.emm.dfc, aes(x = f.id, y = emmean, group = 1)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  # Labels and Theme
  labs(x = "Fish crossing order",
       y = "Distance to conditioned fish") +
  theme_minimal()
print(plt.dfc)

## distance to tank wall

fname_dtw <- paste0('lmer_dttw_mual_0.1_muat_2.3_muesp_1.25_K_', K, '_k_', k, '_sight_', sight, '_gamma_0.25_omegaini_', omega.ini, '.csv')
df.dist.to.tw <- read.csv(fname_dtw, header = F)
col.names <- c("video", "f.id", "dtw")
colnames(df.dist.to.tw) <- col.names

df.dist.to.tw$video <- factor(df.dist.to.tw$video)
df.dist.to.tw$f.id <- factor(df.dist.to.tw$f.id, levels = c("1" ,"2", "3", "4", "5"))

lme.dtw <- lmer(dtw ~ f.id + (1|video), data = df.dist.to.tw)
lme.dtw.null <- lmer(dtw ~ 1 + (1|video), data = df.dist.to.tw)
anova(lme.dtw, lme.dtw.null, test = 'F')

emm.dtw <- emmeans(lme.dtw, ~ f.id)
f.id.comp.dtw <- pairs(emm.dtw, adjust = "tukey")

df.emm.dtw <- as.data.frame(emm.dtw)
df.comp.dtw <- as.data.frame(f.id.comp.dtw)

plt.dtw <- ggplot(df.emm.dtw, aes(x = f.id, y = emmean, group = 1)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  # Labels and Theme
  labs(x = "Fish crossing order",
       y = "Distance to tank wall at t = 0 (cm)") +
  theme_minimal()
print(plt.dtw)

# ggplot(df.comp.dcf, aes(x = contrast, y = estimate)) +
#   # Add a vertical line at 0 (the 'null' hypothesis)
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
#   # Add the points and error bars
#   geom_errorbar(aes(ymin = -SE, ymax = SE), color = "black") +
#   # Labeling
#   labs(title = "Pairwise Comparisons for Distance to CF",
#        x = "Estimate of Difference (dcf)",
#        y = "Contrast (Fish ID Comparison)") +
#   theme_bw()

## analysis for viewing angle

fname_psi <- paste0('lmer_psi_mual_0.1_muat_2.3_muesp_1.25_K_', K, '_k_', k, '_sight_', sight, '_gamma_0.25_omegaini_', omega.ini, '.csv')
df.psi.ic <- read.csv(fname_psi, header = F)
col.names <- c("video", "f.id", "psi")
colnames(df.psi.ic) <- col.names
df.psi.ic$psi <- (180*df.psi.ic$psi)/pi

df.psi.ic$video <- factor(df.psi.ic$video)
df.psi.ic$f.id <- factor(df.psi.ic$f.id, levels = c("2", "3", "4", "5"))

lme.psi <- lmer(psi ~ f.id + (1|video), data = df.psi.ic)
lme.psi.null <- lmer(psi ~ 1 + (1|video), data = df.psi.ic)
anova(lme.psi, lme.psi.null, test = "F")

emm.psi <- emmeans(lme.psi, ~ f.id)
f.id.comp.psi <- pairs(emm.psi, adjust = "tukey")

df.emm.psi <- as.data.frame(emm.psi)

plt.psi <- ggplot(df.emm.psi, aes(x = f.id, y = emmean, group = 1)) + 
  geom_point(size = 3, color = "darkblue") + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) + 
  labs(x = "Fish crossing order",
       y = "Viewing angle (in degree))") +
  theme_minimal()
print(plt.psi)

## analysis for relative orientation

fname_phi <- paste0('lmer_phi_mual_0.1_muat_2.3_muesp_1.25_K_', K, '_k_', k, '_sight_', sight, '_gamma_0.25_omegaini_', omega.ini, '.csv')
df.phi.ic <- read.csv(fname_phi, header = F)
col.names <- c("video", "f.id", "phi")
colnames(df.phi.ic) <- col.names
df.phi.ic$phi <- (180*df.phi.ic$phi)/pi

df.phi.ic$video <- factor(df.phi.ic$video)
df.phi.ic$f.id <- factor(df.phi.ic$f.id, levels = c("2", "3", "4", "5"))

lme.phi <- lmer(phi ~ f.id + (1|video), data = df.phi.ic)
lme.phi.null <- lmer(phi ~ 1 + (1|video), data = df.phi.ic)
anova(lme.phi, lme.phi.null, test = "F")

emm.phi <- emmeans(lme.phi, ~ f.id)
f.id.comp.phi <- pairs(emm.phi, adjust = "tukey")

df.emm.phi <- as.data.frame(emm.phi)

plt.phi <- ggplot(df.emm.phi, aes(x = f.id, y = emmean, group = 1)) + 
  geom_point(size = 3, color = "darkblue") + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) + 
  labs(x = "Fish crossing order",
       y = "Relative orientation (in degree))") +
  theme_minimal()
print(plt.phi)

## analysis for time since last fish crossed

fname_tslfc <- paste0('lmer_tslfc_mual_0.1_muat_1.2_muesp_0.8_K_', K, '_k_', k, '_sight_2.36_gamma_0.25_omegaini_', omega.ini, '.csv')
df.time.last.fish.crossed <- read.csv(fname_tslfc, header = F)
col.names <- c("video", "f.id", "tslfc")
colnames(df.time.last.fish.crossed) <- col.names

df.time.last.fish.crossed$video <- factor(df.time.last.fish.crossed$video)
df.time.last.fish.crossed$f.id <- factor(df.time.last.fish.crossed$f.id, levels = c("1", "2", "3", "4", "5"))

lme.tslfc <- lmer(tslfc ~ f.id + (1|video), data = df.time.last.fish.crossed)
lme.tslfc.null <- lmer(tslfc ~ 1 + (1|video), data = df.time.last.fish.crossed)
anova(lme.tslfc, lme.tslfc.null, test = 'F')

emm.tslfc <- emmeans(lme.tslfc, ~ f.id)
f.id.comp.tslfc <- pairs(emm.tslfc, adjust = "tukey")

df.emm.tslfc <- as.data.frame(emm.tslfc)

plt.tslfc <- ggplot(df.emm.tslfc, aes(x = f.id, y = emmean, group = 1)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  # Labels and Theme
  labs(x = "Fish crossing order",
       y = "Time since previous fish crossed") +
  theme_minimal()
print(plt.tslfc)

# print(f.id.comp.dcf)
# print(f.id.comp.tslfc)
