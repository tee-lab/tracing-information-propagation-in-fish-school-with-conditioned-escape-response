## Code to reproduce fig S2

library(lme4)
library(emmeans)
library(ggplot2)
library(dplyr)

## calculate mean speeds for a given cr for initial phase

df.mean.spd.cr.ini <- read.csv('mean_spd_cr_ini_video_data.csv', header = F)
col.names <- c("video", "f.id", "mean.spd.cr.ini")
colnames(df.mean.spd.cr.ini) <- col.names

df.mean.spd.cr.ini$video <- factor(df.mean.spd.cr.ini$video)
df.mean.spd.cr.ini$f.id <- factor(df.mean.spd.cr.ini$f.id, levels = c("1", "2", "3", "4", "5"))

lme.mean.spd.cr.ini <- lmer(mean.spd.cr.ini ~ f.id + (1|video), data = df.mean.spd.cr.ini)
lme.mean.spd.cr.ini.null <- lmer(mean.spd.cr.ini ~ 1 + (1|video), data = df.mean.spd.cr.ini)
anova(lme.mean.spd.cr.ini, lme.mean.spd.cr.ini.null, test = "F")

emm.spd.cr.ini <- emmeans(lme.mean.spd.cr.ini, ~f.id)
f.id.comp.spd.cr.ini <- pairs(emm.spd.cr.ini, adjust = "tukey")

df.emm.spd.cr.ini <- as.data.frame(emm.spd.cr.ini)
df.comp.spd.cr.ini <- as.data.frame(f.id.comp.spd.cr.ini)

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

## plotting them together

df.sum.spd.ini.relax <- bind_rows(
  as.data.frame(emm.spd.cr.ini) %>% mutate(group_label = "Initial phase"),
  as.data.frame(emm.spd.cr.relax) %>% mutate(group_label = "Relaxation phase")
)

pd = position_dodge(width = 0.5)

plt.spd.ini.relax <- ggplot(df.sum.spd.ini.relax, aes(x = f.id, y = emmean, colour = group_label,
                                                      group = group_label)) +
  geom_point(size = 7, position = pd) + 
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), position = pd, 
                width = 0.2, linewidth = 1.2) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Initial phase" = "#fa4d56",
      "Relaxation phase" = "#33b1ff"
    )
  ) + 
  scale_y_continuous(
    breaks = seq(0, 7, by = 1)
  ) +
  coord_cartesian(ylim = c(0,7)) +
  scale_x_discrete(limits = c("1","2", "3", "4", "5"),
                   expand = expansion(add = 0.3)) +
  labs(x = "Crossing rank",
       y = "Mean speed (cm/s)") +
  theme_bw(base_family = "Helvetica") +
  theme(
    axis.title = element_text(size = 35, color = "black"),
    axis.text = element_text(size = 35, color = "black"),
    legend.text = element_text(size = 30),
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.12),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(colour = "black", linewidth = 2),
    axis.ticks = element_line(colour = "black", linewidth = 1.5),
    axis.ticks.length = unit(-0.2, "cm")
  )

print(plt.spd.ini.relax)

ggsave(filename = "speed_ini_relax_data.pdf", plot = plt.spd.ini.relax, 
       device = cairo_pdf, width = 8, height = 7, units = "in", dpi = 600)



