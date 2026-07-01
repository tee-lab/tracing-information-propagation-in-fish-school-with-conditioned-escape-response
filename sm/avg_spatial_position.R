## Code to calculate average distance to conditioned fish

start.time <- Sys.time()
library(lme4)
library(emmeans)
library(ggplot2)
library(dplyr)
library(ggeffects)

df.dist.to.cf <- read.csv('dist_to_cf_data.csv', header = F)
col.names <- c("video", "f.id", "dcf")
colnames(df.dist.to.cf) <- col.names

df.dist.to.cf$video <- factor(df.dist.to.cf$video)
df.dist.to.cf$f.id <- factor(df.dist.to.cf$f.id, levels = c("2", "3", "4", "5"))

lme.dcf <- lmer(dcf ~ f.id + (1|video), data = df.dist.to.cf)
# lme.dcf <- glmer(dcf ~ f.id + (1|video), family = Gamma(link = "log"), data = df.dist.to.cf)
lme.dcf.null <- lmer(dcf ~ 1 + (1|video), data = df.dist.to.cf)
# lme.dcf.null <- glmer(dcf ~ 1 + (1|video), family = Gamma(link = "log"), data = df.dist.to.cf)
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
       y = "Distance to conditioned fish (cm)") +
  theme_minimal()
print(plt.dfc)

my_color_palette <- c(
  "#DC3B3B", "#DC863B", "#DCD03B", "#9EDC3B", "#54DC3B", "#3BDC6D",
  "#3BDCB7", "#3BB7DC", "#3B6DDC", "#543BDC", "#9E3BDC", "#DC3BD0",
  "#DC3B86", "#EBA092", "#EBC992", "#E5EB92", "#BBEB92", "#92EB92",
  "#92EBBB", "#92EBE5", "#92C9EB", "#92A0EB", "#AE92EB", "#D792EB",
  "#EB92D7", "#EB92AE", "#97411A", "#977A1A", "#7A971A", "#41971A",
  "#1A972D", "#1A9767", "#1A8E97", "#1A5497", "#1A1A97", "#541A97",
  "#8E1A97", "#971A67", "#971A2D"
)
unique_videos <- unique(df.dist.to.cf$video)
auto_video_colors <- setNames(my_color_palette[1:length(unique_videos)], unique_videos)

plt.dfc <- ggpredict(lme.dcf, terms = c("f.id", "video"), type = "random") %>%
  plot() +
  scale_color_manual(values = auto_video_colors) +
  labs(x = "Crossing rank order (f.id)", y = "Distance to conditioned fish",
       title = "", color = "Video") + 
  theme_bw() + 
  theme(
    text = element_text(size = 25),
    legend.text = element_text(size = 10)
  )
print(plt.dfc)

## analysis to distance to tank wall

df.dist.to.tw <- read.csv('dist_to_tw_data.csv', header = F)
col.names <- c("video", "f.id", "dtw")
colnames(df.dist.to.tw) <- col.names

df.dist.to.tw$video <- factor(df.dist.to.tw$video)
df.dist.to.tw$f.id <- factor(df.dist.to.tw$f.id, levels = c("1" ,"2", "3", "4", "5"))

lme.dtw <- lmer(dtw ~ f.id + (1|video), data = df.dist.to.tw)
# lme.dtw <- glmer(dtw ~ f.id + (1|video), family = Gamma(link = "log"), data = df.dist.to.tw)
lme.dtw.null <- lmer(dtw ~ 1 + (1|video), data = df.dist.to.tw)
# lme.dtw.null <- glmer(dtw ~ 1 + (1|video), family = Gamma(link = "log"), data = df.dist.to.tw)
anova(lme.dtw, lme.dtw.null, test = 'F')

emm.dtw <- emmeans(lme.dtw, ~ f.id)
f.id.comp.dtw <- pairs(emm.dtw, adjust = "tukey")

df.emm.dtw <- as.data.frame(emm.dtw)
df.comp.dtw <- as.data.frame(f.id.comp.dtw)

plt.dtw <- ggplot(df.emm.dtw, aes(x = f.id, y = emmean, group = 1)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  # Labels and Theme
  labs(x = "Fish crossing order",
       y = "Distance to tank wall (cm)") +
  theme_minimal()
print(plt.dtw)

## analysis for viewing angle

df.psi <- read.csv('psi_data.csv', header = F)
col.names <- c("video", "f.id", "psi")
colnames(df.psi) <- col.names

df.psi$video <- factor(df.psi$video)
df.psi$f.id <- factor(df.psi$f.id, levels = c("2", "3", "4", "5"))

lme.psi <- lmer(psi ~ f.id + (1|video), data = df.psi)
lme.psi.null <- lmer(psi ~ 1 + (1|video), data = df.psi)
anova(lme.psi, lme.psi.null, test = 'F')

emm.psi <- emmeans(lme.psi, ~ f.id)
f.id.comp.psi <- pairs(emm.psi, adjust = "tukey")

df.emm.psi <- as.data.frame(emm.psi)

plt.psi <- ggplot(df.emm.psi, aes(x = f.id, y = emmean, group = 1)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  # Labels and Theme
  labs(x = "Fish crossing order",
       y = "Viewing angle") +
  theme_minimal()
print(plt.psi)

## analysis for relative orientation

df.phi <- read.csv('phi_data.csv', header = F)
col.names <- c("video", "f.id", "phi")
colnames(df.phi) <- col.names

df.phi$video <- factor(df.phi$video)
df.phi$f.id <- factor(df.phi$f.id, levels = c("2", "3", "4", "5"))

lme.phi <- lmer(phi ~ f.id + (1|video), data = df.phi)
lme.phi.null <- lmer(phi ~ 1 + (1|video), data = df.phi)
anova(lme.phi, lme.phi.null, test = 'F')

emm.phi <- emmeans(lme.phi, ~ f.id)
f.id.comp.phi <- pairs(emm.phi, adjust = "tukey")

df.emm.phi <- as.data.frame(emm.phi)

plt.phi <- ggplot(df.emm.phi, aes(x = f.id, y = emmean, group = 1)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  # Labels and Theme
  labs(x = "Fish crossing order",
       y = "Relative orientation") +
  theme_minimal()
print(plt.phi)

## analysis for time since last fish crossed

df.time.last.fish.crossed <- read.csv('time_since_last_fish_crossed_data.csv', header = F)
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
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  # Labels and Theme
  labs(x = "Fish crossing order",
       y = "Time since previous fish crossed (s)") +
  theme_minimal()
print(plt.tslfc)

plt.tslfc <- ggpredict(lme.tslfc, terms = c("f.id", "video"), type = "random") %>%
  plot() +
  scale_color_manual(values = auto_video_colors) +
  labs(x = "Cross rank order (f.id)", y = "TSLFC", title = "", color = "Video") + 
  theme_bw() +
  theme(
    text = element_text(size = 25),
    legend.text = element_text(size = 10)
  )
print(plt.tslfc)

