#**************************************************************#
##     ------           UE terrain 2026             -------    #

# Master BE - Ecole Pratique des Hautes Etudes #

# Project : ecoacoustics in Vercors landscapes #
# analyse ecoacoustic indices from automatic recording units (ARU)
# indices computed in ephe-terrain-2026-indices-ecoacoustiques.R
# Author : Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
# Contact : Killian Rey (killian.rey@etu.ephe.psl.eu)
#**************************************************************#

## libraries -------------------------------------------------------------------

library(mapview)
library(sf)
library(leaflet)
library(ggspatial)

library(dplyr)
library(tidyverse)
library(patchwork)

library(hms)
library(suncalc)

library(viridis)
library(ade4)
library(factoextra)

library(mgcv)
library(ggeffects)

## data-------------------------------------------------------------------------

indices <- read.csv2("outputs/ephe-terrain-2026-acoustic-indices.csv")
meta <- read.csv2("data/sig/ephe-terrain-2026-coord-sm.csv", dec = ".")

all.acou <- merge(indices,
                  meta,
                  by.x = "SM",
                  by.y = "sm.code",
                  all = T)

# convert START to date and separate date and hour

all.acou$START_dt <- parse_date_time(all.acou$START,
                                     orders = c("dmy_HM", "dmy"),
                                     tz = "Europe/Paris")

all.acou <- all.acou %>%
  mutate(DATE = as.Date(START_dt), HOUR = as_hms(START_dt))

## check acoustic data ---------------------------------------------------------

# number and time range per ARU

synthesis.aru <- aggregate(all.acou$START_dt,
                           by = list(all.acou$SM),
                           FUN = "length")
time.aru.min <- aggregate(all.acou$START_dt,
                          by = list(all.acou$SM),
                          FUN = "min")
time.aru.max <- aggregate(all.acou$START_dt,
                          by = list(all.acou$SM),
                          FUN = "max")

## standardize time range ------------------------------------------------------

# analyze from 05/05/26 15:00 to 08/09/26 10:00

all.acou.sub <- subset(all.acou,
                       START_dt >= "2026-09-05 15:00:00" &
                         START_dt <= "2026-09-08 10:00:00")

## check time series -----------------------------------------------------------

aci.ts <- ggplot(all.acou.sub) +
  aes(x = START_dt, y = ACI) +
  geom_line() +
  geom_point() +
  facet_wrap( ~ placename) +
  theme_minimal()

ndsi.ts <- ggplot(all.acou.sub) +
  aes(x = START_dt, y = NDSI) +
  geom_line() +
  geom_point() +
  facet_wrap( ~ placename) +
  theme_minimal()

h.ts <- ggplot(all.acou.sub) +
  aes(x = START_dt, y = H) +
  geom_line() +
  geom_point() +
  facet_wrap( ~ placename) +
  theme_minimal()

## principal component analysis ------------------------------------------------

all.acou.sub$ID <- paste("SM", 1:nrow(all.acou.sub), sep = "")
rownames(all.acou.sub) <- all.acou.sub$ID

pc.data <- all.acou.sub[, c("BIOAC", "HT", "HF", "H", "ACI", "NDSI")]
pc.acou <- dudi.pca(pc.data, nf = 2, scannf = F)

s.label(pc.acou$li)
s.class(pc.acou$li, fac = factor(all.acou.sub$SM))
s.corcircle(pc.acou$co)

mean(all.acou.sub$ACI)
mean(subset(all.acou.sub, ID == "SM311")$ACI)

# fancy plots for introductory talk

cols <- rev(c(
  "#F2C14E",
  "#6BBF59",
  "#20A39E",
  "#21618C",
  "#3B1F5C"
)
)

fviz_pca_ind(pc.acou,
             label = "none", 
             geom.ind="point",
             pointshape=19,
             habillage = all.acou.sub$habitat, 
             palette = cols,
             addEllipses = TRUE,
             title = "",
             legend.position = "none"
)+
  scale_color_manual(name = "Habitats", values = cols, labels = c("Ouvert","Mosaïque","Clairière","Forêt claire","Forêt dense"))+
   guides(fill = "none")

ggsave(filename = "outputs/ephe-terrain-2026-pca-ind.png", width = 8, height = 6)

fviz_pca_var(pc.acou,
             title = "")

ggsave(filename = "outputs/ephe-terrain-2026-pca-var.png", width = 6, height = 6)

## compute day/night times -----------------------------------------------------

# compute mean and sd per ARU, separate night and day

# compute day / night times
sun_times <- getSunlightTimes(
  data = data.frame(
    date = all.acou.sub$DATE,
    lat = all.acou.sub$y,
    lon = all.acou.sub$x
  ),
  keep = c("sunrise", "sunset"),
  tz = "Europe/Paris"
)

# 2. Injecter les résultats et créer la colonne Jour/Nuit
all.acou.sub <- all.acou.sub %>%
  mutate(
    sunrise = sun_times$sunrise,
    sunset = sun_times$sunset,
    PERIOD = if_else(START_dt >= sunrise &
                       START_dt <= sunset, "day", "night")
  )

## maps for NDSI and ACI per day and night -------------------------------------

# summarize indices

avg.indices <- aggregate(
  all.acou.sub[, c("NDSI", "ACI")],
  by = list(all.acou.sub$SM, all.acou.sub$PERIOD),
  FUN = "mean"
)
colnames(avg.indices) <- c("SM", "PERIOD", "NDSI.avg", "ACI.avg")

sd.indices <- aggregate(
  all.acou.sub[, c("NDSI", "ACI")],
  by = list(all.acou.sub$SM, all.acou.sub$PERIOD),
  FUN = "sd"
)
colnames(sd.indices) <- c("SM", "PERIOD", "NDSI.sd", "ACI.sd")

indices.msd <- merge(avg.indices, sd.indices, by = c("SM", "PERIOD"))

indices.msd.meta <- merge(indices.msd, meta, by.x = "SM", by.y = "sm.code")

# make maps - all combined

indices_sf <- st_as_sf(indices.msd.meta, coords = c("x", "y"), crs = 4326)
indices_sf_day <- subset(indices_sf, PERIOD == "day")
indices_sf_night <- subset(indices_sf, PERIOD == "night")

m.aci.all <- ggplot(indices_sf) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = ACI.avg), size = 4) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "jour+nuit", color = "ACI moyen")+
  theme(legend.position = "none")

m.aci.day <- ggplot(indices_sf_day) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = ACI.avg), size = 4) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "jour", color = "ACI moyen")+
  theme(legend.position = "none")

m.aci.night <- ggplot(indices_sf_night) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = ACI.avg), size = 4) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "nuit", color = "ACI moyen")

m.ndsi.all <- ggplot(indices_sf) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = NDSI.avg), size = 4) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "jour+nuit", color = "NDSI moyen")+
  theme(legend.position = "none")

m.ndsi.day <- ggplot(indices_sf_day) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = NDSI.avg), size = 4) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "jour", color = "NDSI moyen")+
  theme(legend.position = "none")

m.ndsi.night <- ggplot(indices_sf_night) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = NDSI.avg), size = 4) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "nuit", color = "NDSI moyen")

plots.aci <- m.aci.all + m.aci.day + m.aci.night
plots.aci
ggsave(
  filename = "outputs/aci-all-periods-map.png",
  plot = plots.aci,
  width = 12,
  height = 6,
  dpi = 300
)

plots.ndsi <- m.ndsi.all + m.ndsi.day + m.ndsi.night
plots.ndsi
ggsave(
  filename = "outputs/ndsi-all-periods-map.png",
  plot = plots.ndsi,
  width = 12,
  height = 6,
  dpi = 300
)

## boxplots per indices --------------------------------------------------------

all.acou.sub$habitat <- factor(
  all.acou.sub$habitat,
  levels = c("open", "semi-open", "clearing", "clear forest", "dense forest")
)
all.acou.sub$PERIOD <- factor(all.acou.sub$PERIOD)

# translate to french for student's presentation

all.acou.sub <- all.acou.sub %>%
  mutate(habitat_fr = recode(habitat,
                      "open" = "ouvert",
                      "semi-open" = "semi-ouvert",
                      "clearing" = "clairière",
                      "clear forest" = "forêt claire", 
                      "dense forest" = "forêt dense"))


all.acou.sub <- all.acou.sub %>%
  mutate(period_fr = recode(PERIOD,
                             "day" = "jour",
                             "night" = "nuit"))

bx.aci <- ggplot(all.acou.sub) +
  aes(x = habitat_fr, y = ACI) +
  geom_boxplot(fill = "gray90") +
  facet_wrap( ~ period_fr) +
  theme_minimal()+ 
labs(x = "habitat",y="ACI")

bx.ndsi <- ggplot(all.acou.sub) +
  aes(x = habitat, y = NDSI) +
  geom_boxplot(fill = "gray90") +
  facet_wrap( ~ period_fr) +
  theme_minimal()

bx.aci
ggsave("outputs/ephe-terrain-2026-boxplot-aci.png",width = 10, height = 5)

bx.ndsi
ggsave("outputs/ephe-terrain-2026-boxplot-ndsi.png",width = 10, height = 5)

## some stats : do soundscapes differ btw habitats? ----------------------------

mycols <- c("goldenrod", "#440154")

# ACI
aci.lm <- lm(ACI ~ habitat * PERIOD, data = all.acou.sub)

par(mfrow = c(2, 2))
plot(aci.lm)

summary(aci.lm)
p.mod.aci <- ggpredict(aci.lm, terms = c("habitat", "PERIOD"))

plot(p.mod.aci, show_residuals = TRUE) +
  scale_color_manual(values = mycols) +
  scale_fill_manual(values = mycols)
ggsave("outputs/ephe-terrain-2026-lm-aci.png",width = 10, height = 5)


# NDSI

ndsi.lm <- lm(NDSI ~ habitat * PERIOD, data = all.acou.sub)

par(mfrow = c(2, 2))
plot(aci.lm)

summary(ndsi.lm)
p.mod.ndsi <- ggpredict(ndsi.lm, terms = c("habitat", "PERIOD"))

plot(p.mod.ndsi, show_residuals = TRUE) +
  scale_color_manual(values = mycols) +
  scale_fill_manual(values = mycols)
ggsave("outputs/ephe-terrain-2026-lm-ndsi.png",width = 10, height = 5)

## some stats : do time series differ among habitats? --------------------------

all.acou.sub$num.dat <- as.numeric(all.acou.sub$START_dt)

# ACI

ts.aci <- gam(ACI ~ habitat + s(num.dat, by = habitat), data = all.acou.sub)
summary(ts.aci)

ggplot(all.acou.sub, aes(x = START_dt, y = ACI, color = habitat)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(as.numeric(x)),
    se = FALSE,
    linewidth = 1.2
  ) +
  scale_color_viridis_d(option = "viridis", direction = -1) +
  theme_minimal() +
  labs(x = "Date", y = "ACI", color = "Habitat")
ggsave("outputs/ephe-terrain-2026-gam-aci.png",width = 10, height = 5)

# NDSI

ts.ndsi <- gam(NDSI ~ habitat + s(num.dat, by = habitat), data = all.acou.sub)
summary(ts.ndsi)

ggplot(all.acou.sub, aes(x = START_dt, y = NDSI, color = habitat)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(as.numeric(x)),
    se = FALSE,
    linewidth = 1.2
  ) +
  scale_color_viridis_d(option = "viridis", direction = -1) +
  theme_minimal() +
  labs(x = "Date", y = "NDSI", color = "Habitat")
ggsave("outputs/ephe-terrain-2026-gam-ndsi.png",width = 10, height = 5)
