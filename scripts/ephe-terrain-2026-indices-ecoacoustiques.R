#**************************************************************#
##     ------           UE terrain 2026             -------    #

# Master BE - Ecole Pratique des Hautes Etudes #

# Project : ecoacoustics in Vercors landscapes #
# Compute ecoacoustic indices from automatic recording units (ARU)
# Author : Jean-Yves Barnagaud (jean-yves.barnagaud@ephe.psl.eu)
# Contact : Killian Rey (killian.rey@etu.ephe.psl.eu)
#**************************************************************#

## libraries -------------------------------------------------------------------

library(soundecology)
library(tuneR)
library(seewave)

library(lubridate) 
library(beepr) 

# the following R functions are retrieved from https://github.com/agasc/Soundscape-analysis-with-R
source("r-functions/Soundscape-analysis-with-R-master/AcouIndexAlpha.R")

## data ------------------------------------------------------------------------

path.folder <- "data/raw/"
content.folder <- dir(path = path.folder)
path.acou <- content.folder[grep("SMA", content.folder)]

## retrieve date and time ------------------------------------------------------

## compute indices -------------------------------------------------------------


ind.acou <- data.frame()

# loop over all ARU (disable to get just one file)
t1 <- Sys.time()

for(k in 1:length(path.acou)){

# run on one ARU 

#k <- 1

# all audio files on this ARU 

m <- path.acou[k]
path <- paste(path.folder, m, "Data", sep = "/")
audio.wav <- dir(path, pattern = ".wav")

# loop over all audio files on this ARU 

for (j in 1:length(audio.wav)){
  
  # open one file
  
  wave.name <- audio.wav[j]
  wave.path <- paste(path, "/", wave.name, sep = "")
  wave0 <- try(readWave(wave.path))
  
  # corrupt files are listed in a .txt in folder outputs, otherwise move on
  
  if(class(wave0)=="try-error"){
    write(wave.path,paste("outputs/failed_wave_Vercors-26","-",season,".txt",sep=""),append=T)
  } else {
    
    # start date retrieved from file name
    
    dt.start <- nchar(wave.name) - 18
    dt.file <-
      substring(wave.name, first = dt.start, last = nchar(wave.name))
    yr <- substring(dt.file, first = 1, last = 4)
    mth <- substring(dt.file, first = 5, last = 6)
    dy <- substring(dt.file, first = 7, last = 8)
    
    hr <- substring(dt.file, first = 10, last = 11)
    mn <- substring(dt.file, first = 12, last = 13)
    sc <- substring(dt.file, first = 14, last = 15)
    start <- ymd_hms(paste(yr, mth, dy, hr, mn, sc, sep = "-"), tz = "CET")
    
    # compute indices only if file duration = 5 min
    if(seewave::duration(wave0) == 300) { 
     
      s1 <- Sys.time() 
      Result <-
        AcouIndexAlpha(
          wave0,
          stereo = FALSE,
          min_freq = 500,
          max_freq = 12000,
          anthro_min = 500,
          anthro_max = 1500,
          bio_min = 1500,
          bio_max = 12000,
          wl = 512,
          j = 5,
          AcouOccupancy = F,
          Bioac = T,
          Hf = T,
          Ht = T,
          H = T,
          ACI = T,
          AEI_villa = F,
          M = F,
          NDSI = T,
          ADI = T,
          NP = T
        )
      s2 <- Sys.time()
      s2-s1
      
      cat(paste(wave.name,as.character(Sys.time()),"////", sep= " "),file=paste("outputs/outfile",k,".txt",sep=""),append=TRUE)
      
      tp <- Result$Mono_left
      tp$FILE <- wave.name
      tp$START <- start
      tp$SM <-
        substring(wave.name, first = 1, last = nchar(wave.name) - 20)
      
      ind.acou <- rbind(ind.acou, tp)
      
      
    } # if 5mn  
    rm("wave0")
  } # else
} #j

# stop loop over all sites

} # k sites
t2 <- Sys.time()
beep()

## compile data file -----------------------------------------------------------

colnames(ind.acou) <- c("BIOAC",
                        "HT",
                        "HF",
                        "H",
                        "ACI",
                        "NDSI",
                        "ADI",
                        "NP",
                        "FILE",
                        "START",
                        "SM"
)

# write results
write.csv2(ind.acou,file="outputs/ephe-terrain-2026-acoustic-indices.csv",row.names=F)



