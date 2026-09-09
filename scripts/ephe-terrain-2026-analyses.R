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
library(viridis)
library(ggspatial)
library(tidyverse)
library(hms)
library(ade4)
library(suncalc)
library(dplyr)

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

aci.ts <- ggplot(all.acou.sub)+
  aes(x = START_dt, y = ACI)+
  geom_line()+
  geom_point()+
  facet_wrap(~placename)+
  theme_minimal()

ndsi.ts <- ggplot(all.acou.sub)+
  aes(x = START_dt, y = NDSI)+
  geom_line()+
  geom_point()+
  facet_wrap(~placename)+
  theme_minimal()

h.ts <- ggplot(all.acou.sub)+
  aes(x = START_dt, y = H)+
  geom_line()+
  geom_point()+
  facet_wrap(~placename)+
  theme_minimal()

## principal component analysis ------------------------------------------------

all.acou.sub$ID <- paste("SM",1:nrow(all.acou.sub),sep="")
rownames(all.acou.sub) <- all.acou.sub$ID

pc.data <- all.acou.sub[,c("BIOAC","HT","HF","H","ACI","NDSI")]
pc.acou <- dudi.pca(pc.data, nf = 2, scannf = F)

s.label(pc.acou$li)
s.class(pc.acou$li, fac = factor(all.acou.sub$SM))
s.corcircle(pc.acou$co)

mean(all.acou.sub$ACI)
mean(subset(all.acou.sub, ID == "SM311")$ACI)

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
    PERIOD = if_else(START_dt >= sunrise & START_dt <= sunset, "day", "night")
  )

## maps for NDSI and ACI per day and night -------------------------------------

# summarize indices

avg.indices <- aggregate(
  all.acou.sub[, c("NDSI", "ACI")],
  by = list(all.acou.sub$SM, all.acou.sub$PERIOD),
  FUN = "mean"
)
colnames(avg.indices) <- c("SM","PERIOD","NDSI.avg","ACI.avg")

sd.indices <- aggregate(
  all.acou.sub[, c("NDSI", "ACI")],
  by = list(all.acou.sub$SM, all.acou.sub$PERIOD),
  FUN = "sd"
)
colnames(sd.indices) <- c("SM","PERIOD","NDSI.sd","ACI.sd")

indices.msd <- merge(avg.indices,sd.indices, by = c("SM","PERIOD"))

indices.msd.meta <- merge(indices.msd,meta,by.x = "SM",by.y = "sm.code")

# make maps - all combined

indices_sf <- st_as_sf(indices.msd.meta, coords = c("x", "y"), crs = 4326)

m.aci.all <- mapview(
  indices_sf, 
  zcol = "ACI.avg", 
  map.types = "OpenStreetMap",
  col.regions = viridis(10, option = "inferno")
)


s.aci.all <- mapview(
  indices_sf, 
  zcol = "ACI.sd", 
  map.types = "OpenStreetMap",
  col.regions = viridis(10, option = "inferno")
)

m.ndsi.all <- mapview(
  indices_sf, 
  zcol = "NDSI.avg", 
  map.types = "OpenStreetMap",
  col.regions = viridis(10, option = "inferno")
)

s.ndsi.all <- mapview(
  indices_sf, 
  zcol = "NDSI.sd", 
  map.types = "OpenStreetMap",
  col.regions = viridis(10, option = "inferno")
)

mapshot(m.aci.all, file = "outputs/map-aci-mean.png", delay = 5)
mapshot(s.aci.all, file = "outputs/map-aci-sd.png", delay = 5)

# make maps for all / day / night

indices_sf_day <- subset(indices_sf,PERIOD == "day")
indices_sf_night <- subset(indices_sf,PERIOD == "night")

m.aci.all <- ggplot(indices_sf) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = ACI.avg), size = 2) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "period= all", color = "mean ACI")

m.aci.day <- ggplot(indices_sf_day) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = ACI.avg), size = 2) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "period= day", color = "mean ACI")

m.aci.night <- ggplot(indices_sf_night) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = ACI.avg), size = 2) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "period= night", color = "mean ACI")

m.ndsi.all <- ggplot(indices_sf) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = NDSI.avg), size = 2) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "period= all", color = "mean NDSI")

m.ndsi.day <- ggplot(indices_sf_day) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = NDSI.avg), size = 2) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "period= day", color = "mean NDSI")

m.ndsi.night <- ggplot(indices_sf_night) +
  annotation_map_tile(type = "osm", zoom = 13) +
  geom_sf(aes(color = NDSI.avg), size = 2) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "period= night", color = "mean NDSI")

plots.aci <- m.aci.all + m.aci.day + m.aci.night

ggsave(
  filename = "outputs/aci-all-periods-map.png", 
plot = plots.aci,
    width = 12, 
  height = 6, 
  dpi = 300
)

plots.ndsi <- m.ndsi.all + m.ndsi.day + m.ndsi.night

ggsave(
  filename = "outputs/ndsi-all-periods-map.png", 
  plot = plots.ndsi,
  width = 12, 
  height = 6, 
  dpi = 300
)
