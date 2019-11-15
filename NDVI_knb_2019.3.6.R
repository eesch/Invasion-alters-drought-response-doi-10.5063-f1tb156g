##########
# Set up
#########
library("tidyverse")
library("lubridate")
library("zoo")
library("car")
library("lme4")
library("piecewiseSEM")
library("cowplot")
library("scales")
library("gridExtra")
library("nlme")

setwd("/Users/ellen/Desktop/Ellen/2019_Esch et al_NDVI/knb")
raw.ndvi.data <- read_csv("./raw.ndvi.data_esch2019.csv", col_type = cols())
biomass.data <- read_csv("./biomass.data_esch2019.csv", col_type = cols())

#############
# Phenological Dates
# extract dates by fitting cubic splines to the measured NDVI means
#############
emptydataframe <- tibble(
  x = numeric(0),
  first.deriv = numeric(0),
  second.deriv = numeric(0),
  ndvi = numeric(0),
  Date = character(0),
  norm.second.deriv = numeric(0),
  norm.first.deriv = numeric(0),
  sum.derivs = numeric(0),
  dif.derivs = numeric(0),
  change.sum = numeric(0),
  Type = numeric(0),
  gs = numeric(0),
  Plot = numeric(0)
)

plot_summary <- emptydataframe
for (i in c(1:30)) {
  plot <- i
  plot.filter <- raw.ndvi.data %>%
    filter(PLOT == plot)
  spline <- smooth.spline(x = plot.filter$DATE, y = plot.filter$mean) # fit a cubic spline
  ndvi <- as_tibble(predict(spline, as.numeric(ymd("2012-11-16")):as.numeric(ymd("2016-07-28")), deriv = 0)) %>%
    rename(ndvi = y) # create a smoothed NDVI curve between our sampling start and end dates over the 4 years
  first.deriv <- as_tibble(predict(spline, as.numeric(ymd("2012-11-16")):as.numeric(ymd("2016-07-28")), deriv = 1)) %>% # first derivative of fitted spline function
    rename(first.deriv = y) %>%
    mutate(norm.first.deriv = ((1 - -1) / (max(first.deriv) - min(first.deriv)) * (first.deriv - max(first.deriv)) + 1)) # centered and scaled first derivative
  second.deriv <- as_tibble(predict(spline, as.numeric(ymd("2012-11-16")):as.numeric(ymd("2016-07-28")), deriv = 2)) %>% # second derivative of fitted spline function
    rename(second.deriv = y) %>%
    mutate(norm.second.deriv = ((1 - -1) / (max(second.deriv) - min(second.deriv)) * (second.deriv - max(second.deriv)) + 1)) # centered and scaled second derivative
  smoothed.data <- full_join(ndvi, first.deriv, by = "x") %>%
    full_join(second.deriv, by = "x") %>%
    mutate(sum.derivs = norm.second.deriv + norm.first.deriv) %>%
    mutate(dif.derivs = norm.first.deriv - norm.second.deriv) %>%
    mutate(Date = as.Date(x)) %>%
    mutate(change.sum = (sum.derivs - lag(sum.derivs, n = 1L)))

  gs_summary <- emptydataframe
  for (i in c(2013:2016)) {
    gs <- i
    lowndvi <- 0.049
    gs.filter <- smoothed.data %>%
      filter(Date >= (paste((gs - 1), "-9-1", sep = "")), Date <= (paste(gs, "-7-30", sep = "")))
    sos <- if (nrow(gs.filter %>% filter(norm.first.deriv > 0, ndvi > lowndvi)) > 0) { # scaled & centered 1st deriv must be positive #fitted NDVI must be > lower-limit threshold
          (gs.filter %>% filter(norm.first.deriv > 0, ndvi > lowndvi, change.sum < 0) %>% filter(row_number() == 1)) # 3)sum of the scaled and centered 1st & 2nd derivs must indicate a local maxima
    } else { emptydataframe }
    if (nrow(sos) > 0) { sos$Type <- "sos" } else { sos$Type <- character(0) }
    if (nrow(sos) > 0) { sos$gs <- gs } else { sos$gs <- numeric(0) }
    if (nrow(sos) > 0) { sos$Plot <- plot } else { sos$Plot <- character(0) }
    max <- gs.filter %>% filter(Date < (paste(gs, "-07-30", sep = ""))) %>% .[which.max(.$ndvi), ] # maximum NDVI
    if (nrow(max) > 0) { max$Type <- "max" } else { max$Type <- character(0) }
    if (nrow(max) > 0) { max$gs <- gs } else { max$gs <- character(0) }
    if (nrow(max) > 0) { max$Plot <- plot } else { max$Plot <- character(0) }
    eos <- if (nrow(sos) > 0) {
        (gs.filter %>% filter(Date > max$Date, ndvi > lowndvi) %>% .[which.min(.$dif.derivs), ]) # Senescence = min value of the difference between 1st & 2nd derivs, occur after max NDVI, and be > lower cutoff value
        } else { emptydataframe }
    if (nrow(eos) > 0) { eos$Type <- "eos" } else { eos$Type <- character(0) }
    if (nrow(eos) > 0) { eos$gs <- gs } else { eos$gs <- character(0) }
    if (nrow(eos) > 0) { eos$Plot <- plot } else { eos$Plot <- character(0) }
    pheno.dates <- rbind(sos, max, eos)
    gs_summary <- bind_rows(pheno.dates, gs_summary)
  }
  plot_summary <- bind_rows(gs_summary, plot_summary)
}
phenology <- full_join(plot_summary, biomass.data, by = c("Plot", "gs")) %>%
  mutate(YDAY = yday(Date), linDAY = ifelse(YDAY < 213, (YDAY + 365), YDAY)) # must make date a linear variable (linDAY) which does not start over Jan 1 (since growing season occurs over winter-spring); Here if date happened before Aug, we added a "year of dates" (+365)

##########
# Seasonally Integrated NDVI
##########
# date sequence between green-up to senescence
area.dateseq <- tibble(x = numeric(0), Plot = numeric(0), gs = numeric(0))
for (i in c(1:30)) {
  p <- i
  for (g in c(2013:2016)) {
    GS <- g
    area.bounds <- phenology %>%
      filter(gs == GS, Plot == p) %>%
      select(Plot, x, ndvi, Type, gs, TREAT, CC) %>%
      filter(Type == "sos" | Type == "eos")
    if (nrow(area.bounds) > 0) {
      (bounds <- tibble(value = seq(area.bounds$x[area.bounds$Type == "sos"], (area.bounds$x[area.bounds$Type == "eos"]))))
    } else {
      (bounds <- tibble(value = 1))
    }
    dateseq <- bounds %>%
      mutate(Plot = p, gs = GS) %>%
      rename(x = value)
    area.dateseq <- bind_rows(dateseq, area.dateseq)
  }
}

# date sequence between green-up to senescence which has NDVI above the lower-limit threshold of 0.049
smoothed.ndvi <- tibble(x = numeric(0), ndvi = numeric(0), Date = character(0), Plot = numeric(0))
for (i in c(1:30)) {
  p <- i
  plot.filter <- raw.ndvi.data %>%
    filter(PLOT == p)
  spline <- smooth.spline(x = plot.filter$DATE, y = plot.filter$mean)
  ndvi <- as.tibble(predict(spline, as.numeric(ymd("2012-11-16")):as.numeric(ymd("2016-07-28")), deriv = 0)) %>%
    rename(ndvi = y) %>%
    mutate(Date = as.Date(x), Plot = p)
  smoothed.ndvi <- bind_rows(ndvi, smoothed.ndvi)
}
green.dateseq <- full_join(smoothed.ndvi, area.dateseq, by = c("x", "Plot")) %>%
  filter(gs != "NA") %>%
  filter(ndvi > lowndvi)

# seasonally integrated ndvi is approximated by looking at rectanges with a time interval (width) of 1 day, and height of the mean ndvi between x and x+1.
area_under_curve <- tibble(area = numeric(0), gsl = numeric(0), gs = numeric(0), Plot = numeric(0))
for (i in c(1:30)) {
  p <- i
  for (g in c(2013:2016)) {
    GS <- g
    greendates <- green.dateseq %>%
      filter(Plot == p, gs == GS)
    data.frame(greendates)
    datediff <- tibble(datediff = diff(greendates$x)) %>%
      mutate(datediff = ifelse(datediff > 1, 0, 1)) # if there is a break in "green dates", need to exclude that rectange
    if (nrow(greendates) > 0) {
      (AUC <- tibble(area = sum(datediff * rollmean((greendates$ndvi - lowndvi), 2)), gsl = sum(datediff), gs = g, Plot = p))
    } else {
      (AUC <- tibble(area = 0, gsl = 0, gs = g, Plot = p))
    }
    area_under_curve <- bind_rows(AUC, area_under_curve)
  }
}

# #####
# Join phenological dates, max ndvi and seasonally intergrated NDVI for final data
# #####
lineardates <- phenology %>%
  select(Plot, gs, Type, linDAY) %>%
  filter(Type != "max") %>%
  spread(Type, linDAY)
sosDATE <- phenology %>%
  filter(Type == "sos") %>%
  select(Date, Plot, gs) %>%
  rename(SOSDATE = Date)
eosDATE <- phenology %>%
  filter(Type == "eos") %>%
  select(Date, Plot, gs) %>%
  rename(EOSDATE = Date)
finaldata <- phenology %>%
  filter(Type == "max") %>%
  select(Plot, gs, ndvi) %>%
  full_join(area_under_curve, by = c("Plot", "gs")) %>%
  full_join(biomass.data, by = c("Plot", "gs")) %>%
  full_join(lineardates, by = c("Plot", "gs")) %>%
  full_join(sosDATE, by = c("Plot", "gs")) %>%
  full_join(eosDATE, by = c("Plot", "gs"))

# #########
# Table 1
# #########
# MaxNDVI
Anova(lmer(ndvi ~ AnRF * CC + (1 | Plot), data = filter(finaldata)), type = 2)
aggregate(ndvi ~ CC, data = finaldata, FUN = mean)
0.4840218 / 0.2364530 # shrub 2x higher ndvi

Anova(lmer(ndvi ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "shrub")), type = 2)
summary(lmer(ndvi ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "shrub")))
0.0008489 * 63 + 0.3613733 # = ndvi MIN for natives
0.0008489 * 237.95 + 0.3613733 # = ndvi MAX for natives
0.414854 / 0.5633691 # maintain 74% of max

Anova(lmer(ndvi ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "herbaceous")), type = 2)
summary(lmer(ndvi ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "herbaceous")))
0.0023817 * 63 - 0.1076597 # = ndvi MIN for invaded
0.0023817 * 237.95 - 0.1076597 # = ndvi MAX for invaded
0.0423874 / 0.4590658 # maintain 9% of max


# IntegratedNDVI
Anova(lmer(area ~ AnRF * CC + (1 | Plot), data = (finaldata)), type = 2)
aggregate(area ~ CC, data = finaldata, FUN = mean)
59.72376 / 13.76652 # shrubs 4.3x greater

Anova(lmer(area ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "shrub")), type = 2)
summary(lmer(area ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "shrub")))

Anova(lmer(area ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "herbaceous")), type = 2)
summary(lmer(area ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "herbaceous")))


# Green-up (start of season = sos)
Anova(lmer(sos ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "shrub")), type = 2)
summary(lmer(sos ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "shrub")))

Anova(lmer(sos ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "herbaceous", !is.na(sos))), type = 2)
summary(lmer(sos ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "herbaceous", !is.na(sos))))


#Senescence (end of season = eos)
Anova(lmer(eos ~ AnRF + (1 | Plot), data = (finaldata %>% filter(CC == "shrub"))), type = 2)
summary(lmer(eos ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "shrub")))
as.Date(482)

Anova(lmer(eos ~ AnRF + (1 | Plot), data = (finaldata %>% filter(CC == "herbaceous", !is.na(eos)))), type = 2)
summary(lmer(eos ~ AnRF + (1 | Plot), data = filter(finaldata, CC == "herbaceous", !is.na(eos))))
as.Date(460)

aggregate(EOSDATE ~ AnRF + CC, data = filter(finaldata, AnRF == 63 | AnRF == 237.95), FUN = mean)

# ######
# SEM 
# ######
# shrub
shrub <- finaldata %>%
  filter(CC == "shrub") %>%
  select(Plot, gs, AnRF, sos, eos, ndvi, area) %>%
  rename(max = ndvi)

shrub_1 <- psem(
  lme(sos ~ AnRF, random = ~ 1 | Plot, data = shrub),
  lme(eos ~ AnRF, random = ~ 1 | Plot, data = shrub),
  lme(max ~ AnRF + sos, ~ 1 | Plot, data = shrub),
  lme(area ~ AnRF + max + sos + eos, random = ~ 1 | Plot, data = shrub)
)
summary(shrub_1) # add eos ~ sos

shrub_2 <- psem(
  lme(sos ~ AnRF, random = ~ 1 | Plot, data = shrub),
  lme(eos ~ AnRF + sos, random = ~ 1 | Plot, data = shrub),
  lme(max ~ AnRF + sos, ~ 1 | Plot, data = shrub),
  lme(area ~ AnRF + max + sos + eos, random = ~ 1 | Plot, data = shrub)
)
AIC(shrub_1, aicc = T) - AIC(shrub_2, aicc = T) # fit improved
summary(shrub_2) # remove sos~anrf

shrub_3 <- psem(
  lme(eos ~ AnRF + sos, random = ~ 1 | Plot, data = shrub),
  lme(max ~ AnRF + sos, ~ 1 | Plot, data = shrub),
  lme(area ~ AnRF + max + sos + eos, random = ~ 1 | Plot, data = shrub)
)
AIC(shrub_2, aicc = T) - AIC(shrub_3, aicc = T) # fit improved
summary(shrub_3) # nothing to remove

summary(shrub_3) # accept final shrub model

# herbaceous
grass <- finaldata %>%
  filter(CC == "herbaceous", !is.na(sos)) %>%
  select(Plot, gs, AnRF, sos, eos, ndvi, area) %>%
  rename(max = ndvi)

grass_1 <- psem(
  lme(sos ~ AnRF, random = ~ 1 | Plot, data = grass),
  lme(eos ~ AnRF, random = ~ 1 | Plot, data = grass),
  lme(max ~ AnRF + sos, ~ 1 | Plot, data = grass),
  lme(area ~ AnRF + max + sos + eos, random = ~ 1 | Plot, data = grass)
)
summary(grass_1) # remove max~sos

grass_2 <- psem(
  lme(sos ~ AnRF, random = ~ 1 | Plot, data = grass),
  lme(eos ~ AnRF, random = ~ 1 | Plot, data = grass),
  lme(max ~ AnRF, ~ 1 | Plot, data = grass),
  lme(area ~ AnRF + max + sos + eos, random = ~ 1 | Plot, data = grass)
)
AIC(grass_1, aicc = T) - AIC(grass_2, aicc = T) # fit NOT improved

summary(grass_1) # accept final grass model
	
#########
# Figure 1
#########
Fig1A <- (full_join(raw.ndvi.data, (biomass.data %>% 
                                      filter(gs == 2013) %>% 
                                      select(Plot, TREAT, CC) %>% 
                                      rename(PLOT = Plot)), by = "PLOT") %>%
  group_by(TREAT, CC, DATE) %>%
  summarise(Mean = mean(mean), Se = sd(mean) / sqrt(n())) %>%
  filter(CC == "herbaceous") %>%
  ggplot(aes(x = ymd(DATE), y = Mean, col = as.factor(TREAT), shape = as.factor(TREAT))) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - Se, ymax = Mean + Se), width = 0) +
  theme_cowplot() +
  labs(x = "Date", y = "Herbaceous NDVI") +
  scale_color_manual(
    values = c("brown3", "darkolivegreen3", "cornflowerblue"),
    breaks = c("150", "100", "50"),
    labels = c("Normal \nprecipitation\n", "Moderate \ndrought\n", "Severe \ndrought\n"),
    name = "Rainfall\nTreatment"
  ) +
  scale_x_date(labels = date_format("%d %b %Y"), 
               breaks = seq(as.Date("2013-1-1"), as.Date("2016-5-1"), by = "12 months")) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, .6), 
                     limits = c(-.09, .685)) +
  geom_line(aes(color = factor(TREAT))) +
  guides(col = F, shape = F) +
  scale_shape_manual(
    values = c(15, 16, 17),
    breaks = c("150", "100", "50"),
    labels = c("Normal \nprecipitation\n", "Moderate \ndrought\n", "Severe \ndrought\n"),
    name = "Rainfall\nTreatment"
  )) %>%
  ggdraw() +
  draw_label("(a)", .04, .975, fontface = "italic", size = 16)

Fig1B <- (full_join(raw.ndvi.data, (biomass.data %>% 
                                      filter(gs == 2013) %>% 
                                      select(Plot, TREAT, CC) %>% 
                                      rename(PLOT = Plot)), by = "PLOT") %>%
  group_by(TREAT, CC, DATE) %>%
  summarise(Mean = mean(mean), Se = sd(mean) / sqrt(n())) %>%
  filter(CC == "shrub") %>%
  ggplot(aes(x = ymd(DATE), y = Mean, col = as.factor(TREAT), shape = as.factor(TREAT))) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Mean - Se, ymax = Mean + Se)) +
  theme_cowplot() +
  labs(x = "Date", y = "Shrub NDVI") +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6), 
                     limits = c(-.09, .685)) +
  scale_color_manual(
    values = c("brown3", "darkolivegreen3", "cornflowerblue"),
    breaks = c("150", "100", "50"),
    labels = c("Normal \nprecipitation\n", "Moderate \ndrought\n", "Severe \ndrought\n"),
    name = "Rainfall\nTreatment"
  ) +
  scale_x_date(labels = date_format("%d %b %Y"), 
               breaks = seq(as.Date("2013-1-1"), as.Date("2016-5-1"), by = "12 months")) +
  geom_line(aes(color = factor(TREAT))) +
  scale_shape_manual(
    values = c(15, 16, 17),
    breaks = c("150", "100", "50"),
    labels = c("Normal \nprecipitation\n", "Moderate \ndrought\n", "Severe \ndrought\n"),
    name = "Rainfall\nTreatment"
  ) +
  guides(col = F, shape = F)) %>%
  ggdraw() +
  draw_label("(b)", .04, .975, fontface = "italic", size = 16)

Fig1Leg <- (full_join(raw.ndvi.data, (biomass.data %>% 
                                        filter(gs == 2013) %>% 
                                        select(Plot, TREAT, CC) %>% 
                                        rename(PLOT = Plot)), by = "PLOT") %>%
  group_by(TREAT, CC, DATE) %>%
  summarise(Mean = mean(mean), Se = sd(mean) / sqrt(n())) %>%
  filter(CC == "shrub") %>%
  ggplot(aes(x = ymd(DATE), y = Mean, col = as.factor(TREAT), shape = as.factor(TREAT))) +
  geom_point(size = 5) +
  scale_color_manual(
    values = c("brown3", "darkolivegreen3", "cornflowerblue"),
    breaks = c("150", "100", "50"),
    labels = c("Normal \nprecipitation\n", "Moderate \ndrought\n", "Severe \ndrought\n"),
    name = "Rainfall\nTreatment"
  ) +
  scale_shape_manual(
    values = c(15, 16, 17),
    breaks = c("150", "100", "50"),
    labels = c("Normal \nprecipitation\n", "Moderate \ndrought\n", "Severe \ndrought\n"),
    name = "Rainfall\nTreatment"
  )) %>%
  get_legend()

grid.arrange(Fig1A, Fig1Leg, Fig1B, ncol = 2, nrow = 2, layout_matrix = rbind(c(1, 2), c(3, 2)), widths = c(1, .15), heights = c(1, 1))

#######
# Figure 2
#######
Fig2a <- (ggplot(finaldata %>%
    group_by(gs, AnRF, CC) %>%
    summarise(MEAN = mean(ndvi), SE = (sd(ndvi) / sqrt(n()))),
  aes(x = AnRF, y = MEAN, col = CC) ) +
  geom_point(size = 3, aes(shape = factor(gs))) +
  theme_cowplot() +
  labs(x = "Rainfall (October - April, mm)", y = "Maximum NDVI") +
  scale_color_manual(values = c("orange2", "darkgreen"), breaks = c("herbaceous", "shrub")) +
  theme(strip.text.x = element_text(face = "bold")) +
  scale_shape_manual(values = c(16, 15, 18, 17), breaks = c(2013, 2014, 2015, 2016)) +
  geom_smooth(fill = NA, method = "lm") +
  geom_errorbar(aes(ymin = MEAN - SE, ymax = MEAN + SE), width = 2) +
  theme(legend.position = "none")) %>% 
  ggdraw() +
  draw_label("(a)", .04, .975, fontface = "italic", size = 16)

Fig2c <- (ggplot(finaldata %>%
    group_by(gs, AnRF, CC) %>%
    summarise(MEAN = mean(area), SE = (sd(area) / sqrt(n()))),
  aes(x = AnRF, y = MEAN, col = CC)) +
  geom_point(size = 3, aes(shape = factor(gs))) +
  theme_cowplot() +
  labs(x = "Rainfall (October - April, mm)", y = "Seasonally integrated  NDVI") +
  scale_color_manual(values = c("orange2", "darkgreen"), breaks = c("herbaceous", "shrub")) +
  theme(strip.text.x = element_text(face = "bold")) +
  scale_shape_manual(values = c(16, 15, 18, 17), breaks = c(2013, 2014, 2015, 2016)) +
  geom_smooth(fill = NA, method = "lm") +
  geom_errorbar(aes(ymin = MEAN - SE, ymax = MEAN + SE), width = 2) +
  theme(legend.position = "none")) %>%
  ggdraw() +
  draw_label("(c)", .04, .975, fontface = "italic", size = 16)

Fig2b <- (ggplot(finaldata %>%
  filter(!is.na(sos)) %>%
  group_by(TREAT, CC, gs, AnRF) %>%
  summarise(MEAN = mean(sos), SE = (sd(sos) / sqrt(n()))),
  aes(x = AnRF, y = MEAN, col = CC, fill = CC)) +
  geom_point(size = 3, aes(shape = as.factor(gs))) +
  geom_smooth(fill = NA, method = "lm", lty = 1, inherit.aes = F, aes(x = AnRF, y = sos),
    col = "orange2", data = (finaldata %>% filter(CC == "herbaceous", !is.na(sos))) ) +
  geom_errorbar(aes(ymin = (MEAN - SE), ymax = (MEAN + SE)), width = 2) +
  labs(y = "Green-up date", x = "Rainfall (October - April, mm)") +
  scale_color_manual(values = c("orange2", "darkgreen"), breaks = c("herbaceous", "shrub")) +
  scale_shape_manual(values = c(21, 22, 23, 24), breaks = c(2013, 2014, 2015, 2016)) +
  scale_fill_manual(values = c("orange2", "darkgreen"), breaks = c("herbaceous", "shrub")) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = c(275, 336, 397), labels = c("Oct.1", "Dec.1", "Feb.1"))) %>%
  ggdraw() +
  draw_label("(b)", .04, .975, fontface = "italic", size = 16)

Fig2d <- (ggplot(finaldata %>%
  filter(!is.na(eos)) %>%
  group_by(TREAT, CC, gs, AnRF) %>%
  summarise(MEAN = mean(eos), SE = (sd(eos) / sqrt(n()))),
  aes(x = AnRF, y = MEAN, col = CC, fill = CC)) +
  geom_point(size = 3, aes(shape = as.factor(gs))) +
  geom_smooth(fill = NA, method = "lm") +
  geom_errorbar(aes(ymin = (MEAN - SE), ymax = (MEAN + SE)), width = 2) +
  labs(y = "Senescence date", x = "Rainfall (October - April, mm)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("orange2", "darkgreen"), breaks = c("herbaceous", "shrub")) +
  scale_shape_manual(values = c(21, 22, 23, 24), breaks = c(2013, 2014, 2015, 2016)) +
  scale_fill_manual(values = c("orange2", "darkgreen"), breaks = c("herbaceous", "shrub")) +
  guides(fill = F) +
  scale_y_continuous(breaks = c(397, 457, 487, 518), labels = c("Feb.1", "Apr.1", "May1", "Jun.1"))) %>%
  ggdraw() +
  draw_label("(d)", .04, .975, fontface = "italic", size = 16)

legend_f2 <- get_legend(ggplot(finaldata %>%
  group_by(gs, AnRF, CC) %>%
  summarise(Mean = mean(area), SE = (sd(area) / sqrt(n()))), aes(x = AnRF, y = Mean, col = CC)) +
  geom_point(size = 3, aes(shape = factor(gs))) +
  scale_color_manual(
    name = "Community\nComposition",
    values = c("orange2", "darkgreen"),
    breaks = c("herbaceous", "shrub"),
    labels = c("Herbaceous", "Shrub")
  ) +
  theme(strip.text.x = element_text(face = "bold")) +
  scale_shape_manual(values = c(16, 15, 18, 17), breaks = c(2013, 2014, 2015, 2016), name = "Growing\nSeason"))

FIG2NEW<-grid.arrange(Fig2a, Fig2c,legend_f2,Fig2b,Fig2d,ncol=3,nrow=2, layout_matrix=rbind(c(1:3),c(4,5,3)),widths=c(2,2,.5),heights=c(3,3))


##############
## Supplemental  Figure 1
###############
rain <- read_csv("./rainfall_esch2019.csv", col_type = cols()) %>%
  mutate(date = mdy(date)) %>%
  mutate(m = month(ymd(date)), year = year(ymd(date)), gsm = paste(m, "1", year, sep = "-")) %>%
  mutate(gs = ifelse(m >= 1 & m <= 4, year, ifelse(m >= 10 & m <= 12, (year + 1), NA))) %>%
  mutate(numevents = ifelse(mm.100 > 0, 1, 0)) %>%
  mutate(numericmonth = NA)
rain$numericmonth[rain$m == "10"] <- 1 # oct
rain$numericmonth[rain$m == "11"] <- 2 # nov
rain$numericmonth[rain$m == "12"] <- 3 # dec
rain$numericmonth[rain$m == "1"] <- 4 # jan
rain$numericmonth[rain$m == "2"] <- 5 # feb
rain$numericmonth[rain$m == "3"] <- 6 # mar
rain$numericmonth[rain$m == "4"] <- 7 # apr

AmtMonthRF <- rain %>%
  filter(STATION == "Manipulation") %>%
  group_by(gs, numericmonth) %>%
  summarise(t100 = sum(mm.100)) %>%
  spread(numericmonth, t100) %>%
  replace(., is.na(.), 0) %>%
  gather(numericmonth, amount, -gs)

AmtMonthRF %>% group_by(gs) %>% summarise(sum(amount, na.rm = T)) %>% data.frame()

Sup1A <- (ggplot(AmtMonthRF,
  aes(x = as.numeric(numericmonth), y = amount, col = as.factor(gs), fill = as.factor(gs), shape = as.factor(gs))) +
  geom_point(size = 4) + geom_line() +
  scale_color_manual(values = c("#e41a1c", "#1f78b4", "#984ea3", "black"),
    labels = c("2013 = 156.80 mm", "2014 = 125.67 mm", "2015 = 139.38 mm", "2016 = 156.42 mm") ) +
  scale_fill_manual(values = c("#e41a1c", "#1f78b4", "#984ea3", "black"),
    labels = c("2013 = 156.80 mm", "2014 = 125.67 mm", "2015 = 139.38 mm", "2016 = 156.42 mm") ) +
  scale_shape_manual(values = c(21:24),
    labels = c("2013 = 156.80 mm", "2014 = 125.67 mm", "2015 = 139.38 mm", "2016 = 156.42 mm") ) +
  labs(x = "", y = "Precipitation (mm)", shape = "", fill = "", col = "") +
  theme(legend.position = c(.5, .9)) +
  scale_x_continuous(breaks = c(2, 3, 4, 5, 6, 7, 8),
    labels = c("Oct.", "Nov.", "Dec.", "Jan.", "Feb.", "Mar.", "Apr.") ) ) %>%
  ggdraw() + draw_label("(a)", .07, .975, fontface = "italic", size = 16)	

EventsMonthRF <- rain %>%
  filter(STATION == "Manipulation") %>%
  group_by(gs, numericmonth) %>%
  summarise(events = sum(numevents)) %>%
  spread(numericmonth, events) %>%
  replace(., is.na(.), 0) %>%
  gather(numericmonth, events, -gs)

EventsMonthRF %>% group_by(gs) %>% summarise(sum(events, na.rm = T)) %>% data.frame()

Sup1B <- (ggplot(EventsMonthRF, 
  aes(x = as.numeric(numericmonth), y = events, col = as.factor(gs), fill = as.factor(gs), shape = as.factor(gs))) +
  geom_point(size = 4) + geom_line() +
  scale_color_manual(values = c("#e41a1c", "#1f78b4", "#984ea3", "black"),
    labels = c("2013 = 12 events", "2014 = 13 events", "2015 = 12 events", "2016 = 12 events") ) +
  scale_fill_manual(
    values = c("#e41a1c", "#1f78b4", "#984ea3", "black"),
    labels = c("2013 = 12 events", "2014 = 13 events", "2015 = 12 events", "2016 = 12 events") ) +
  scale_shape_manual(values = c(21:24),
    labels = c("2013 = 12 events", "2014 = 13 events", "2015 = 12 events", "2016 = 12 events") ) +
  labs(x = "", y = "\nNumber of rain events", shape = "", fill = "", col = "") +
  theme(legend.position = c(.5, .9)) +
  scale_x_continuous(breaks = c(2, 3, 4, 5, 6, 7, 8),
    labels = c("Oct.", "Nov.", "Dec.", "Jan.", "Feb.", "Mar.", "Apr.") ) +
  scale_y_continuous(breaks = c(0:6)) ) %>%
  ggdraw() +
  draw_label("(b)", .07, .975, fontface = "italic", size = 16)

smerrain <- rain %>%
  filter(STATION == "Manipulation") %>%
  group_by(gs) %>%
  summarise(mm50 = sum(mm.50), mm100 = sum(mm.100), mm150 = sum(mm.150)) %>%
  gather(treatment, rain, -gs) %>%
  mutate(date = paste(gs, "-01-01", sep = ""))

yearhistoric <- rain %>%
  filter(STATION != "Manipulation") %>%
  group_by(gs, STATION) %>%
  summarise(rain = sum(mm.100), N = n()) %>%
  group_by(gs) %>%
  summarise(rain = mean(rain), N = sum(!is.na(STATION))) %>%
  mutate(treatment = "mm100") %>%
  mutate(date = paste(gs, "-01-01", sep = ""), gs = as.numeric(gs))

avline <- mean(yearhistoric$rain, na.rm = T) ; avline


Sup1C <- ggdraw(bind_rows(yearhistoric, smerrain) %>%
  filter(treatment == "mm100") %>%
  ggplot(aes(x = (ymd(date)), y = rain)) +
  geom_bar(stat = "identity", na.rm = T) +
  labs(x = "Growing Season year", y = "Biomass (g/m2)") +
  labs(x = "Growing season year", y = "Growing season rainfall\n(mm between Oct. - Apr.)") +
  geom_point(inherit.aes = FALSE, data = filter(smerrain), aes(x = ymd(date), y = rain, col = treatment, fill = treatment, size = 2, shape = treatment)) +
  scale_fill_manual(values = c("darkorange1", "gold1", "brown3"),
    breaks = c("mm150", "mm100", "mm50"),
    labels = c("Normal\nprecipitation", "Moderate\ndrought", "Severe\ndrought"),
    name = "Rainfall Treatment") +
  scale_color_manual(values = c("darkorange1", "gold1", "brown3"),
    breaks = c("mm150", "mm100", "mm50"),
    labels = c("Normal\nprecipitation", "Moderate\ndrought", "Severe\ndrought"),
    name = "Rainfall Treatment") +
  scale_shape_manual(values = c(21, 22, 23),
    breaks = c("mm150", "mm100", "mm50"),
    labels = c("Normal\nprecipitation", "Moderate\ndrought", "Severe\ndrought"),
    name = "Rainfall Treatment") +
  guides(size = F, colour = guide_legend(override.aes = list(size = 5))) +
  geom_hline(yintercept = avline)) +
  draw_label("(c)", .04, .975, fontface = "italic", size = 16)


grid.arrange(Sup1A, Sup1B, Sup1C, ncol = 2, nrow = 2, layout_matrix = rbind(c(1, 2), c(3)), widths = c(1, 1), heights = c(1, 1))

#############################
# Supplemental Fig 3, 
# NDVI ~ biomass
#############################
# smer herbaceous with NDVI
finaldata %>% filter(Biomass_g.m2 == 0) %>% summarise(mean(ndvi)) # lowndvi = 0.049 #cutoff threshold

# stats
# max ndvi ~ biomass
Anova(lm(ndvi ~ log(Biomass_g.m2) * CC, data = (finaldata %>% filter(Biomass_g.m2 > 0))), type = 2)
summary(lm(ndvi ~ log(Biomass_g.m2), data = (finaldata %>% filter(Biomass_g.m2 > 0))))


# area ~ biomass
Anova(lm(area ~ log(Biomass_g.m2) * CC, data = (finaldata %>% filter(Biomass_g.m2 > 0))), type = 2)

Anova(lm(area ~ log(Biomass_g.m2), data = (finaldata %>% filter(Biomass_g.m2 > 0, CC == "shrub"))), type = 2)
summary(lm(area ~ log(Biomass_g.m2), data = (finaldata %>% filter(Biomass_g.m2 > 0, CC == "shrub"))))

Anova(lm(area ~ log(Biomass_g.m2), data = (finaldata %>% filter(Biomass_g.m2 > 0, CC == "herbaceous"))), type = 2)
summary(lm(area ~ log(Biomass_g.m2), data = (finaldata %>% filter(Biomass_g.m2 > 0, CC == "herbaceous"))))


# ndvi~stem elongation
Anova(lm(ndvi ~ log(StemElongation_mm.month), data = finaldata))
LogElong <- finaldata %>% mutate(logElongation = log(StemElongation_mm.month))
cor.test(LogElong$ndvi, LogElong$logElongation)
cor.test(LogElong$area, LogElong$logElongation)


# fig
f1 <- ggdraw(ggplot(
  (finaldata %>% filter(Biomass_g.m2 > 0)),
  aes(y = ndvi, x = (Biomass_g.m2), col = CC) ) +
  geom_smooth(method = "lm", fill = NA, col = "black") +
  geom_point(aes(col = CC)) +
  labs(shape = "Species", y = "Maximum NDVI", x = expression(Biomass~production~(g/m^{2}))) +
  guides(col = F) +
  scale_color_manual(values = c("orange2", "darkgreen"), name = "Community\ncomposition", labels = c("Invaded", "Native")) +
  scale_x_continuous(trans = log_trans(), breaks = c(3, 5, 10, 50, 100, 500, 2000))) +
  draw_label("(a)", .04, .975, fontface = "italic", size = 16)

f2 <- ggdraw(ggplot((finaldata %>% filter(Biomass_g.m2 > 0)), aes(y = area, x = Biomass_g.m2, col = CC)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA) +
  labs(shape = "Species", y = "Seasonally integrated NDVI", x = expression(Biomass~production~(g/m^{2}))) +
  scale_color_manual(values = c("orange2", "darkgreen"), name = "Community\nComposition", labels = c("Exotic \nannual\n", "Perennial \nshrub")) +
  guides(col = F) +
  scale_x_continuous(trans = log_trans(), breaks = c(3, 5, 10, 50, 100, 500, 2000))) +
  draw_label("(c)", .04, .975, fontface = "italic", size = 16)

fsupLEG <- get_legend(ggplot((finaldata %>% filter(Biomass_g.m2 > (0))), aes(y = area, x = Biomass_g.m2, col = CC)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA) + labs(shape = "Species", y = "Seasonally integrated NDVI", x = expression(Biomass~production~(g/m^{2}))) +
  scale_color_manual(values = c("orange2", "darkgreen"), name = "Community\nComposition", labels = c("Herbaceous", "Shrub")) +
  scale_x_continuous(trans = log_trans(), breaks = c(3, 5, 10, 50, 100, 500, 2000)))

f3 <- ggdraw(ggplot(finaldata %>% filter(!is.na(StemElongation_mm.month)), aes(y = ndvi, x = StemElongation_mm.month)) +
  geom_point(col = "darkgreen") +
  geom_smooth(method = "lm", col = "darkgreen", fill = NA) +
  labs(x = "Monthly stem growth (mm)", y = "Maximum NDVI") +
  scale_x_continuous(trans = log_trans(), breaks = c(15, 25, 50, 100))) +
  draw_label("(b)", .04, .975, fontface = "italic", size = 16)

f4 <- ggdraw(ggplot(finaldata %>% filter(!is.na(StemElongation_mm.month)), aes(y = area, x = StemElongation_mm.month)) +
  geom_point(col = "darkgreen") +
  geom_smooth(method = "lm", col = "darkgreen", fill = NA) +
  labs(x = "Monthly stem growth (mm)", y = "Seasonally integrated NDVI") +
  scale_x_continuous(trans = log_trans(), breaks = c(15, 25, 50, 100))) +
  draw_label("(d)", .04, .975, fontface = "italic", size = 16)
	
grid.arrange(f1, f2, fsupLEG, f3, f4, ncol = 3, nrow = 2, layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6)), widths = c(1, 1, .3), heights = c(1, 1))
