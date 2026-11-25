## Code to calculate average distance to conditioned fish

start.time <- Sys.time()
library(lme4)
library(emmeans)
library(ggplot2)
library(dplyr)
library(ggeffects)

# DISTANCE TO CONDITIONED FISH

# data

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

df.emm.dcf <- as.data.frame(emm.dcf)
df.comp.dcf <- as.data.frame(f.id.comp.dcf)

df.dist.to.cf.x <- read.csv('dist_to_cf_x_data.csv', header = F)
col.names <- c("video", "f.id", "dcf.x")
colnames(df.dist.to.cf.x) <- col.names

df.dist.to.cf.x$video <- factor(df.dist.to.cf.x$video)
df.dist.to.cf.x$f.id <- factor(df.dist.to.cf.x$f.id, levels = c("2", "3", "4", "5"))

lme.dcf.x <- lmer(dcf.x ~ f.id + (1|video), data = df.dist.to.cf.x)
# lme.dcf.x <- glmer(dcf.x ~ f.id + (1|video), family = Gamma(link = "log"), data = df.dist.to.cf.x)
lme.dcf.x.null <- lmer(dcf.x ~ 1 + (1|video), data = df.dist.to.cf.x)
# lme.dcf.x.null <- glmer(dcf.x ~ 1 + (1|video), family = Gamma(link = "log"), data = df.dist.to.cf.x)
anova(lme.dcf.x, lme.dcf.x.null, test = 'F')

emm.dcf.x <- emmeans(lme.dcf.x, ~ f.id)
f.id.comp.dcf.x <- pairs(emm.dcf.x, adjust = "tukey")

df.emm.dcf.x <- as.data.frame(emm.dcf.x)
df.comp.dcf.x <- as.data.frame(f.id.comp.dcf.x)

# model

df.dist.to.cf.model <- read.csv('dist_to_cf_model.csv', header = F)
col.names <- c("video", "f.id", "dcf")
colnames(df.dist.to.cf.model) <- col.names

df.dist.to.cf.model$video <- factor(df.dist.to.cf.model$video)
df.dist.to.cf.model$f.id <- factor(df.dist.to.cf.model$f.id, levels = c("2", "3", "4", "5"))

lme.dcf.model <- lmer(dcf ~ f.id + (1|video), data = df.dist.to.cf.model)
# lme.dcf.model <- glmer(dcf ~ f.id + (1|video), family = Gamma(link = "log"), data = df.dist.to.cf.model)
lme.dcf.model.null <- lmer(dcf ~ 1 + (1|video), data = df.dist.to.cf.model)
# lme.dcf.model.null <- glmer(dcf ~ 1 + (1|video), family = Gamma(link = "log"), data = df.dist.to.cf.model)
anova(lme.dcf.model, lme.dcf.model.null, test = 'F')

emm.dcf.model <- emmeans(lme.dcf.model, ~ f.id)
f.id.comp.dcf.model <- pairs(emm.dcf.model, adjust = "tukey")

df.emm.dcf.model <- as.data.frame(emm.dcf.model)
df.comp.dcf.model <- as.data.frame(f.id.comp.dcf.model)

df.dist.to.cf.x.model <- read.csv('dist_to_cf_x_model.csv', header = F)
col.names <- c("video", "f.id", "dcf.x")
colnames(df.dist.to.cf.x.model) <- col.names

df.dist.to.cf.x.model$video <- factor(df.dist.to.cf.x.model$video)
df.dist.to.cf.x.model$f.id <- factor(df.dist.to.cf.x.model$f.id, levels = c("2", "3", "4", "5"))

lme.dcf.x.model <- lmer(dcf.x ~ f.id + (1|video), data = df.dist.to.cf.x.model)
# lme.dcf.x.model <- glmer(dcf.x ~ f.id + (1|video), family = Gamma(link = "log"), data = df.dist.to.cf.x.model)
lme.dcf.x.model.null <- lmer(dcf.x ~ 1 + (1|video), data = df.dist.to.cf.x.model)
# lme.dcf.x.model.null <- glmer(dcf.x ~ 1 + (1|video), family = Gamma(link = "log"), data = df.dist.to.cf.x.model)
anova(lme.dcf.x.model, lme.dcf.x.model.null, test = 'F')

emm.dcf.x.model <- emmeans(lme.dcf.x.model, ~ f.id)
f.id.comp.dcf.x.model <- pairs(emm.dcf.x.model, adjust = "tukey")

df.emm.dcf.x.model <- as.data.frame(emm.dcf.x.model)
df.comp.dcf.x.model <- as.data.frame(f.id.comp.dcf.x.model)

## Plotting distance to conditioned fish as a function of crossing rank

# 1. Combine all four dataframes with series labels and a shared dodge_group per pair

df.sum.dist.cf <- bind_rows(
  as.data.frame(df.emm.dcf) %>% mutate(series = "2D", dodge_group = "Data"),
  as.data.frame(df.emm.dcf.x) %>% mutate(series = "X", dodge_group = "Data"),
  as.data.frame(df.emm.dcf.model) %>% mutate(series = "2D Model", dodge_group = "Model"),
  as.data.frame(df.emm.dcf.x.model) %>% mutate(series = "X Model", dodge_group = "Model")
)

# 2. Position dodge separates Pair 1 (1st & 2nd) from Pair 2 (3rd & 4th)
pd <- position_dodge(width = 0.5)

# 3. Plot
plt.dfc <- ggplot(df.sum.dist.cf, aes(x = f.id, y = emmean, color = series, group = dodge_group)) +
  geom_point(aes(shape = series), size = 7, position = pd) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = pd) +
  
  scale_color_manual(
    name = NULL,
    values = c(
      "2D" = "#fa4d56", # 1st dataframe hex
      "X" = "#fa4d56",
      "2D Model" = "#33b1ff",
      "X Model" = "#33b1ff" # 2nd dataframe hex
    ),
    breaks = c("2D", "X"),
    labels = c(
      "2D" = expression(r[C]^i),
      "X" = expression(x[C]^i)
    )
  ) +
  
  # Same shape for 1st & 3rd (16 = circle), same shape for 2nd & 4th (17 = triangle)
  scale_shape_manual(
    name = NULL,
    values = c(
      "2D" = 16,
      "X" = 17,
      "2D Model" = 16,
      "X Model" = 17
    ),
    breaks = c("2D", "X"),
    labels = c(
      "2D" = expression(r[C]^i),
      "X" = expression(x[C]^i)
    )
  ) +
  
  coord_cartesian(ylim = c(0, 6)) +
  scale_x_discrete(limits = c("2", "3", "4", "5"),
                   expand = expansion(add = 0.3)) +
  labs(x = "Crossing rank",
       y = "Distance to conditioned fish (cm)") +
  theme_bw(base_family = "Helvetica") +
  theme(
    
    axis.title = element_text(size = 28, color = "black"),
    axis.text = element_text(size = 40, color = "black"),
    legend.text = element_text(size = 40),
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.12),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(colour = "black", linewidth = 2),
    axis.ticks = element_line(colour = "black", linewidth = 1.5),
    axis.ticks.length = unit(-0.2, "cm")
  )

print(plt.dfc)

ggsave(filename = "dist_cf.pdf", plot = plt.dfc, device = cairo_pdf, width = 8, 
       height = 7, units = "in", dpi = 600)

## analysis to distance to tank wall

# data

df.dist.to.tw <- read.csv('dist_to_tw_data.csv', header = F)
col.names <- c("video", "f.id", "dtw")
colnames(df.dist.to.tw) <- col.names

df.dist.to.tw$video <- factor(df.dist.to.tw$video)
df.dist.to.tw$f.id <- factor(df.dist.to.tw$f.id, levels = c("2", "3", "4", "5"))

lme.dtw <- lmer(dtw ~ f.id + (1|video), data = df.dist.to.tw)
# lme.dtw <- glmer(dtw ~ f.id + (1|video), family = Gamma(link = "log"), data = df.dist.to.tw)
lme.dtw.null <- lmer(dtw ~ 1 + (1|video), data = df.dist.to.tw)
# lme.dtw.null <- glmer(dtw ~ 1 + (1|video), family = Gamma(link = "log"), data = df.dist.to.tw)
anova(lme.dtw, lme.dtw.null, test = 'F')

emm.dtw <- emmeans(lme.dtw, ~ f.id)
f.id.comp.dtw <- pairs(emm.dtw, adjust = "tukey")

df.emm.dtw <- as.data.frame(emm.dtw)
df.comp.dtw <- as.data.frame(f.id.comp.dtw)

# model

df.dist.to.tw.model <- read.csv('dist_to_tw_model.csv', header = F)
col.names <- c("video", "f.id", "dtw")
colnames(df.dist.to.tw.model) <- col.names

df.dist.to.tw.model$video <- factor(df.dist.to.tw.model$video)
df.dist.to.tw.model$f.id <- factor(df.dist.to.tw.model$f.id, levels = c("2", "3", "4", "5"))

lme.dtw.model <- lmer(dtw ~ f.id + (1|video), data = df.dist.to.tw.model)
# lme.dtw.model <- glmer(dtw ~ f.id + (1|video), family = Gamma(link = "log"), data = df.dist.to.tw.model)
lme.dtw.model.null <- lmer(dtw ~ 1 + (1|video), data = df.dist.to.tw.model)
# lme.dtw.model.null <- glmer(dtw ~ 1 + (1|video), family = Gamma(link = "log"), data = df.dist.to.tw.model)
anova(lme.dtw.model, lme.dtw.model.null, test = 'F')

emm.dtw.model <- emmeans(lme.dtw.model, ~ f.id)
f.id.comp.dtw.model <- pairs(emm.dtw.model, adjust = "tukey")

df.emm.dtw.model <- as.data.frame(emm.dtw.model)
df.comp.dtw.model <- as.data.frame(f.id.comp.dtw.model)

# plotting distance to tankwall

df.sum.dist.tw <- bind_rows(
  as.data.frame(df.emm.dtw) %>% mutate(group_label = "Experiments"),
  as.data.frame(df.emm.dtw.model) %>% mutate(group_label = "Model")
)

plt.dtw <- ggplot(df.sum.dist.tw, aes(x = f.id,y = emmean, 
                                      color = group_label, group = group_label)) +
  geom_point(size = 7, position = pd) + 
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = pd) + 
  scale_color_manual(
    name = NULL,
    values = c(
      "Experiments" = "#fa4d56",
      "Model" = "#33b1ff"
    )
  ) +
  coord_cartesian(ylim = c(0,17)) +
  scale_x_discrete(limits = c("2", "3", "4", "5"),
                   expand = expansion(add = 0.3)) +
  labs(x = "Crossing rank",
       y = "Shortest distance to hurdle (cm)") +
  theme_bw(base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 28, color = "black"),
    axis.text = element_text(size = 40, color = "black"),
    # legend.text = element_text(size = 30),
    # legend.position = "inside",
    # legend.position.inside = c(0.5, 0.12),
    # legend.direction = "vertical",
    # legend.background = element_rect(fill = "transparent", color = NA),
    # legend.key = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(colour = "black", linewidth = 2),
    axis.ticks = element_line(colour = "black", linewidth = 1.5),
    axis.ticks.length = unit(-0.2, "cm"),
    legend.position = "none"
  )

print(plt.dtw)

ggsave(filename = "dist_tw.pdf", plot = plt.dtw, device = cairo_pdf, width = 8, 
       height = 7, units = "in", dpi = 600)

## analysis for viewing angle

# data

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
df.comp.psi <- as.data.frame(f.id.comp.psi)

# model

df.psi.model <- read.csv('psi_model.csv', header = F)
col.names <- c("video", "f.id", "psi")
colnames(df.psi.model) <- col.names

df.psi.model$video <- factor(df.psi.model$video)
df.psi.model$f.id <- factor(df.psi.model$f.id, levels = c("2", "3", "4", "5"))

lme.psi.model <- lmer(psi ~ f.id + (1|video), data = df.psi.model)
lme.psi.model.null <- lmer(psi ~ 1 + (1|video), data = df.psi.model)
anova(lme.psi.model, lme.psi.model.null, test = 'F')

emm.psi.model <- emmeans(lme.psi.model, ~ f.id)
f.id.comp.psi.model <- pairs(emm.psi.model, adjust = "tukey")

df.emm.psi.model <- as.data.frame(emm.psi.model)
df.comp.psi.model <- as.data.frame(f.id.comp.psi.model)

# plotting viewing angle to conditioned fish

df.sum.psi <- bind_rows(
  as.data.frame(df.emm.psi) %>% mutate(group_label = "Experiments"),
  as.data.frame(df.emm.psi.model) %>% mutate(group_label = "Model")
)

plt.psi <- ggplot(df.sum.psi, aes(x = f.id, y = emmean, 
                                  color = group_label, group = group_label)) +
  geom_point(size = 7, position = pd) + 
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = pd) + 
  scale_color_manual(
    name = NULL,
    values = c(
      "Experiments" = "#fa4d56",
      "Model" = "#33b1ff"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 120, by = 20)
  ) +
  coord_cartesian(ylim = c(0,125)) +
  scale_x_discrete(limits = c("2", "3", "4", "5"),
                   expand = expansion(add = 0.3)) +
  labs(x = "Crossing rank",
       y = "Viewing angle (in degree)") +
  theme_bw(base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 28, color = "black"),
    axis.text = element_text(size = 30, color = "black"),
    # legend.text = element_text(size = 30),
    # legend.position = "inside",
    # legend.position.inside = c(0.5, 0.12),
    # legend.direction = "vertical",
    # legend.background = element_rect(fill = "transparent", color = NA),
    # legend.key = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(colour = "black", linewidth = 2),
    axis.ticks = element_line(colour = "black", linewidth = 1.5),
    axis.ticks.length = unit(-0.2, "cm"),
    legend.position = "none"
  )

print(plt.psi)

ggsave(filename = "psi.pdf", plot = plt.psi, device = cairo_pdf, width = 8, 
       height = 7, units = "in", dpi = 600)

## analysis for relative orientation

# data

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
df.comp.phi <- as.data.frame(f.id.comp.phi)

# model

df.phi.model <- read.csv('phi_model.csv', header = F)
col.names <- c("video", "f.id", "phi")
colnames(df.phi.model) <- col.names

df.phi.model$video <- factor(df.phi.model$video)
df.phi.model$f.id <- factor(df.phi.model$f.id, levels = c("2", "3", "4", "5"))

lme.phi.model <- lmer(phi ~ f.id + (1|video), data = df.phi.model)
lme.phi.model.null <- lmer(phi ~ 1 + (1|video), data = df.phi.model)
anova(lme.phi.model, lme.phi.model.null, test = 'F')

emm.phi.model <- emmeans(lme.phi.model, ~ f.id)
f.id.comp.phi.model <- pairs(emm.phi.model, adjust = "tukey")

df.emm.phi.model <- as.data.frame(emm.phi.model)
df.comp.phi.model <- as.data.frame(f.id.comp.phi.model)

## plotting relative orientation to conditioned fish

df.sum.phi <- bind_rows(
  as.data.frame(df.emm.phi) %>% mutate(group_label = "Experiments"),
  as.data.frame(df.emm.phi.model) %>% mutate(group_label = "Model")
  )

plt.phi <- ggplot(df.sum.phi, aes(x = f.id, y = emmean, color = group_label, group = group_label)) +
  geom_point(size = 7, position = pd) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = pd) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Experiments" = "#fa4d56",
      "Model" = "#33b1ff"
    )
  ) + 
  scale_y_continuous(
    breaks = seq(0,120,by=20)
  ) +
  coord_cartesian(ylim = c(0,125)) +
  scale_x_discrete(limits = c("2", "3", "4", "5"),
                   expand = expansion(add = 0.3)) +
  labs(x = "Crossing rank",
       y = "Relative orientation (in degree)") +
  theme_bw(base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 28, color = "black"),
    axis.text = element_text(size = 30, color = "black"),
    # legend.text = element_text(size = 30),
    # legend.position = "inside",
    # legend.position.inside = c(0.5, 0.12),
    # legend.direction = "vertical",
    # legend.background = element_rect(fill = "transparent", color = NA),
    # legend.key = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(colour = "black", linewidth = 2),
    axis.ticks = element_line(colour = "black", linewidth = 1.5),
    axis.ticks.length = unit(-0.2, "cm"),
    legend.position = "none"
  )

print(plt.phi)

ggsave(filename = "phi.pdf", plot = plt.phi, device = cairo_pdf, width = 8,
       height = 7, units = "in", dpi = 600)

## relative orientation to tank wall

df.phi.tw <- read.csv('phi_tw_data.csv', header = F)
col.names <- c("video", "f.id", "phi.tw")
colnames(df.phi.tw) <- col.names

df.phi.tw$video <- factor(df.phi.tw$video)
df.phi.tw$f.id <- factor(df.phi.tw$f.id, levels = c("2", "3", "4", "5"))

lme.phi.tw <- lmer(phi.tw ~ f.id + (1|video), data = df.phi.tw)
lme.phi.tw.null <- lmer(phi.tw ~ 1 + (1|video), data = df.phi.tw)
anova(lme.phi.tw, lme.phi.tw.null, test = 'F')

emm.phi.tw <- emmeans(lme.phi.tw, ~ f.id)
f.id.comp.phi.tw <- pairs(emm.phi.tw, adjust = "tukey")

df.emm.phi.tw <- as.data.frame(emm.phi.tw)
df.comp.phi.tw <- as.data.frame(f.id.comp.phi.tw)

## plotting viewing angle to tankwall

df.sum.phi.tw <- bind_rows(
  as.data.frame(df.emm.phi.tw) %>% mutate(group_label = "Experiments")
)

plt.phi.tw <- ggplot(df.sum.phi.tw, aes(x = f.id, y = emmean, color = group_label, group = group_label)) +
  geom_point(size = 7, position = pd) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.2, position = pd) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Experiments" = "#fa4d56",
      "Model" = "#33b1ff"
    )
  ) + 
  scale_y_continuous(
    breaks = seq(0,120,by=20)
  ) +
  coord_cartesian(ylim = c(0,125)) +
  scale_x_discrete(limits = c("2", "3", "4", "5"),
                   expand = expansion(add = 0.3)) +
  labs(x = "Crossing rank",
       y = "Viewing angle to hurdle (in degree)") +
  theme_bw(base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 28, color = "black"),
    axis.text = element_text(size = 30, color = "black"),
    # legend.text = element_text(size = 30),
    # legend.position = "inside",
    # legend.position.inside = c(0.5, 0.12),
    # legend.direction = "vertical",
    # legend.background = element_rect(fill = "transparent", color = NA),
    # legend.key = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(colour = "black", linewidth = 2),
    axis.ticks = element_line(colour = "black", linewidth = 1.5),
    axis.ticks.length = unit(-0.2, "cm"),
    legend.position = "none"
  )

print(plt.phi.tw)

ggsave(filename = "phi_tw.pdf", plot = plt.phi.tw, device = cairo_pdf, width = 8,
       height = 7, units = "in", dpi = 600)