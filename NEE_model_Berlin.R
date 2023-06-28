#=======================================================================================
#
# Title:       Net Ecosystem Exchange of CO2.
# Author:      Dr. Max Anjos (maxanjos@campus.ul.pt)
# Description: More details on the approach are available at:
#             https://github.com/ByMaxAnjos/CO2-vegetation-emissions
# Data: 09.11.2022
#
#' @param vegetation_fraction (required). Raster with fraction of vegetation
#' @param weather (required) Hourly meteorological. Global radiation in wm2 and air temperatura ºC
#' @return NEE, GPP, and Reco fluxes of CO2: csv. table, plots and sf multipolylines and raster map with 100 meters of resolution.
#' @examples
#=======================================================================================
#Install packages
# Load theses packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, sf, data.table, raster, tmap, osmdata, gstat, automap, openair, recipes, timetk)
library(tidyverse)
library(data.table)
library(tmap)
library(osmdata)
library(sf)
library(raster)
library(gstat)   # The most popular R-Package for Kriging (imho)
library(automap) # Automatize some (or all) parts of the gstat-workflow
library(openair)
library(recipes)
library(timetk)
library(terra)
library(stars)
library(tidymodels)

setwd("~/Documents/CO2CityMap/Berlin/Components/Vegetation/")

#Increasing the memory before the calculation (windows)
# mymemory <- memory.limit() * 3
# memory.limit(size = mymemory)

source("~/Documents/CO2CityMap/Berlin/Components/Transport/CO2-traffic-emissions/R/ZCCM_functions.R")

#================================================================
#Define the region of int/project system in UTM
#================================================================

#get ROI
city <- "Berlin"
mycrs <- "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"

study_area <- osmdata::getbb(city, format_out = "sf_polygon", limit = 1)$multipolygon %>%
  st_as_sf() %>% # convert to sf object
  st_transform(crs = mycrs)
study_area <- st_make_valid(study_area)
qtm(study_area) #plot map

st_write(study_area, "data/area_berlin.shp")

#================================================================
#Load LCZ map. Download the LCZ global map from https://zenodo.org/record/6364594/files/lcz_filter_v1.tif?download=1
#================================================================
lcz_map <- raster("~/Documents/CO2CityMap/Berlin/Components/Building/LCZ/lcz_berlin2.tif")
lcz_map <- raster::projectRaster(lcz_map, crs = mycrs)

#=================================================================
## Load and pre-processing rasters
#=================================================================

raster_veg <- raster("data/veg_frac_30m_ber.TIF")
raster_veg <- raster::projectRaster(raster_veg, crs = mycrs)
qtm(raster_veg, style = "cobalt")

#calculate vegetation fraction from raster with 1m resolution
# raster_veg[raster_veg>0] <- 1
# veg_fraction <- terra::aggregate(x = raster_veg, fact = 30, cores = 5,
#                                  fun = function(x, ...)
#                                    {(sum(x==1, ...)/900)*100})
# veg_ras = raster::crop(veg_fraction, raster::extent(study_area))
# veg_ras= raster::mask(veg_ras, study_area)
# raster::writeRaster(raster_veg, "data/veg_1m_ber.TIF", format="GTiff", overwrite = TRUE)
# mapview::mapview(veg_ras, alpha = 0.4)

#=================================================================
## Load air temperature data
#=================================================================

#UCON data
air_UCON <- fread("~/Documents/CO2CityMap/Berlin/Components/Building/inputs/data/airT_UCON_2015_2022.csv") %>%
  rename(longitude = Longitude,
         latitude = Latitude)

#=================================================================
## Global radiation
#=================================================================
Rg_data <- fread("data/GR_berlin_ERA5_2005_2020.csv")

importGRmax <- function(file = myfilename){
  gr_df <- janitor::clean_names(file)
  gr_df <- str_split_fixed(file$time, ":", 2)
  gr_df <-  tibble("date" = unlist(gr_df[,1]),
                   "hour" = unlist(gr_df[,2])) %>%
    mutate(date=as.Date(anytime::anydate(date), format="%Y-%m-%d")) %>%
    mutate(date= as.POSIXct(paste(date, hour),format="%Y-%m-%d %H"))
  my_gr <- bind_cols(file, gr_df) %>% janitor::clean_names() %>%
    rename(Rg = g_i, sun = h_sun, airT = t2m, ws = ws10m) %>%
    dplyr::select(date, Rg, sun, airT, ws)

  return(my_gr)
}
Rg_data <- importGRmax(file = Rg_data) %>%
  dplyr::select(-airT)


#=================================================================
# Calculate NEE
#=================================================================
#Define the period
#imonth <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
imonth <- "jul"
iyear <- c(2018)
myday <- 1
imput <- expand.grid(imonth, iyear)

NeeUHI <- function(imput, roi = study_area, veg_df = veg_ras, rg_df = Rg_data, air_df = air_UCON, spRes = 100,
                   NEE = TRUE, GPP = TRUE, Reco = TRUE) {

  imonth <- imput[1]
  iyear <- imput[2]

  #Select Rg data
  if(is.null(myday)) {
    rg_model <- rg_df %>%
      selectByDate(year = iyear, month = imonth) %>%
      mutate(Rg = ifelse(Rg < 5, 0, Rg))
  } else {
    rg_model <- rg_df %>%
      selectByDate(year = iyear, month = imonth, day = myday) %>%
      mutate(Rg = ifelse(Rg < 5, 0, Rg))
  }

  #Select air data
  if(is.null(myday)) {
    air_model <- air_df %>%
      openair::selectByDate(year = iyear, month = imonth)
  } else {
    air_model <- air_df %>%
      openair::selectByDate(year = iyear, month = imonth, day = myday)
  }

  meteo_model <- inner_join(air_model, rg_model, by="date")

  veg_model <- terra::as.data.frame(veg_df, xy=TRUE, na.rm = TRUE) %>% set_names(c("x", "y", "veg_fraction")) %>%
    #filter(veg_fraction >10) %>%
    mutate(FID=row_number())

    # resample iLCZ
    # Convert to spatial pixel
    grid_stations <- raster(raster::extent(roi), res = spRes)
    raster::values(grid_stations) <- 1:raster::ncell(grid_stations)
    raster::crs(grid_stations) <- mycrs

    iLCZ <- raster::resample(lcz_map, grid_stations)
      # raster::crop(lcz_map, raster::extent(roi)) %>%
      # raster::mask(roi)

    # rename iLCZ to "lcz"
    names(iLCZ) <- "lcz"

    # convert lcz_map to polygon
    lcz_shp <- terra::as.polygons(rast(iLCZ)) %>%
      st_as_sf() %>%
      st_transform(crs = mycrs)

    # convert to list of polygons
    poly_list <- lapply(st_geometry(lcz_shp), function(x)
      st_sfc(x, crs = mycrs))

    # calculate areas of each polygon
    lcz_areas <- lapply(1:length(poly_list), function(i)
      st_area(poly_list[[i]]))

    # sample the number of points according to area of the polygon and convert to data.frame
    lcz_poi <- lapply(1:length(poly_list), function(i)
      st_sample(poly_list[[i]], size=100, prob=lcz_areas[[i]], method = "random", exact = FALSE) %>%
        as.data.frame())

    # merge all dataframes and convert to sf points
    lcz_poi <- do.call(rbind.data.frame, lcz_poi) %>%
      st_as_sf() %>%
      st_transform(crs = mycrs)

    # intersect lcz poi with lcz shp
    lcz_poi_get <- sf::st_intersection(lcz_poi, lcz_shp)

  #Calculate by day
  iday <- meteo_model %>%
    mutate(day = lubridate::day(date)) %>%
    distinct(day, .keep_all = FALSE) %>%
    expand.grid()

    NEEday <- function(iday) {

      myday <- iday[1]

      mod_day <- meteo_model %>%
        mutate(day = lubridate::day(date)) %>%
        openair::selectByDate(day = myday)

      #Downscale to hour
      ihour <- mod_day %>%
        mutate(ihour = lubridate::hour(date)) %>%
        distinct(ihour, .keep_all = FALSE) %>%
        expand.grid()

      NEEhour <- function(ihour) {

        myhour <- ihour[1]

        rad_model <- mod_day %>%
          mutate(hour = lubridate::hour(date)) %>%
          openair::selectByDate(hour = myhour)

        date_cross <-  crossing(stuff = veg_model$FID, rad_model$date) %>%
          set_names("FID", "date")
        setDT(date_cross)
        date_cross[,FID:=as.numeric(FID)]
        setDT(veg_model)
        veg_model[,FID:=as.numeric(FID)]

        fcover_model <- inner_join(date_cross, veg_model, by="FID", all.x = TRUE, all.y = FALSE)
        fcover_model <- inner_join(fcover_model, rad_model %>% dplyr::select(-day), by="date", all.x = TRUE, all.y = FALSE)

        #Get station points
        shp_stations <- rad_model %>%
          distinct(latitude, longitude, .keep_all = TRUE) %>%
          dplyr::select(station, latitude, longitude, airT) %>%
          sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
          st_transform(crs = mycrs)

        #Intersect shp stations with lcz shp
        lcz_stations <- sf::st_intersection(shp_stations, lcz_shp)
        # #Create train and test sets
        # lcz_stations_split <- initial_split(lcz_stations)
        # lcz_stations_train <- training(lcz_stations_split)
        # lcz_stations_test <- testing(lcz_stations_split)

        #Merge table to create a sf to model
        lcz_poi_mod <- inner_join(lcz_poi_get,lcz_stations %>% as_tibble() %>%
                                    dplyr::select(lcz, station, airT), by=c("lcz")) %>%
          group_by(station) %>%
          mutate(FID_station=cur_group_id()) %>%
          dplyr::select(-station, -lcz)

        # Convert LCZ map to starts
        lcz_stars <- st_as_stars(iLCZ, dimensions = "XY")
        train_mod = st_join(lcz_poi_mod, st_as_sf(lcz_stars)) %>%
          mutate(lcz = as.integer(lcz))
        train_mod = train_mod[!is.na(train_mod$lcz),]
        st_crs(lcz_stars) <- st_crs(train_mod)

        #IDW
        idw_mod = gstat(formula = airT ~ 1, data = train_mod)
        idw_map = predict(idw_mod, lcz_stars)
        idw_map = idw_map["var1.pred",,]
        idw_map <- raster(rast(idw_map))
        mydate <- rad_model %>% distinct(date, .keep_all = FALSE) %>% as_tibble()
        mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
        names(idw_map) <- paste0("IDW_airT_",mydate)

        #return(idw_map)

        # Resample air map
        air_resample = raster::resample(idw_map, veg_df)
        #raster::writeRaster(air_resample, paste0("Building/outputs/maps/UHI/", mydate, "_UHI.TIF"), format="GTiff", overwrite = TRUE)
        #Merge building with air raster
        air_poi <- as_tibble(rasterToPoints(air_resample)) %>% set_names(c("x", "y", "airT"))
        #veg_poi <- as_tibble(rasterToPoints(veg_df)) %>% set_names(c("x", "y", "veg_fraction"))
        veg_cal <- inner_join(fcover_model %>% dplyr::select(x, y, Rg, veg_fraction), air_poi, by= c("x", "y")) %>%
          mutate(hour = paste0(myhour))
        #Calculate the NEE CO2 fluxes
        #NEE_data <- ECO2NEE_2(fcover_model)

        #Prepare the function
        ##Get Rre
          night_period <- meteo_model %>%
            filter(Rg <10)
          # Calculate seasonality of temperature with rolling window
          temp_seasonal <- zoo::rollapply(data = night_period$airT, width = 7, FUN = mean, by = 4)
          # Define Q10 function for respiration
          q10 <- function(T, Tref, Q10) {Q10 ^ ((T - Tref)/10)}
          # Set parameters for Q10 function
          Tref <- 15  # reference temperature in Celsius
          Q10 <- 2.0  # assumed change in respiration per 10 degrees Celsius increase in temperature
          # Calculate Rref for each time point using Q10 function
          Rref <- q10(temp_seasonal, Tref, Q10)
          Rref <- abs(mean(Rref))

        NEEboost <- function(x) {
          x <- x %>% #filter(veg_fraction > 15) %>%
            mutate(B = abs(-8.25 + 0.35*veg_fraction)
                   ,a = 0.0002*veg_fraction + 0.0052 # 0.005 + 0.016*veg_fraction
                   #,yc = if_else(Rg < 20, REddyProc::fLloydTaylor(10, 330, airT + 273.15)[1], 0)
                   ,Reco = ifelse(Rg <10, Rref * (exp(230 * (1/(15 - -46.02) - 1/(airT - -46.02))))*(veg_fraction/100), 0)
                   #,Reco = 1.3 * (exp(230 * (1/(15 - -46.02) - 1/(airT - -46.02))))*(veg_fraction/100)
                   #,yc = 1.3 * (exp(230 * (1/(15 - -46.02) - 1/(airT - -46.02))))*(veg_fraction/100)
                   ,o =0.96
                   ,PAR = 0.505*Rg
                   ,NEE = Reco - 1/2*o*(a*PAR+B-((a*PAR+B)^2 -4*a*B*o*PAR)^0.5)
                   ,GPP = -1/2*o*(a*PAR+B-((a*PAR+B)^2 -4*a*B*o*PAR)^0.5)
                   #,GPP = ifelse(GPP1>0, 0, GPP1)
                   #,NEE = GPP-Reco
                   )
          return(x)
        }

        NEE_data <- NEEboost(veg_cal)

        if(NEE == TRUE) {
          #Get NEE
          NEE_ras <- NEE_data %>% dplyr::select(x, y, NEE)
          NEE_ras = raster::rasterFromXYZ(xyz = NEE_ras,crs = mycrs)
          NEE_ras = raster::crop(NEE_ras, raster::extent(roi))
          NEE_ras= raster::mask(NEE_ras, roi)
          mydate <- rad_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
          mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
          names(NEE_ras) <- paste0("NEE_", mydate)
          #raster::writeRaster(NEE_ras, paste0("Vegetation/output/maps/", mydate, "_NEE.TIF"), format="GTiff", overwrite = TRUE)
          #return(NEE_ras)
        }
        if(GPP == TRUE) {
          #Get GPP
          GPP_ras <- NEE_data %>% dplyr::select(x, y, GPP)
          GPP_ras = raster::rasterFromXYZ(xyz = GPP_ras,crs = mycrs)
          GPP_ras = raster::crop(GPP_ras, raster::extent(roi))
          GPP_ras= raster::mask(GPP_ras, roi)
          mydate <- rad_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
          mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
          names(GPP_ras) <- paste0("GPP_", mydate)
          #raster::writeRaster(GPP_ras, paste0("Vegetation/output/maps/", mydate, "_NEE.TIF"), format="GTiff", overwrite = TRUE)
          #return(GPP_ras)
        }
        if(Reco ==TRUE) {

          #Get Reco
          Reco_ras <- NEE_data %>% dplyr::select(x, y, Reco)
          Reco_ras = raster::rasterFromXYZ(xyz = Reco_ras,crs = mycrs)
          Reco_ras = raster::crop(Reco_ras, raster::extent(roi))
          Reco_ras= raster::mask(Reco_ras, roi)
          mydate <- rad_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
          mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
          names(Reco_ras) <- paste0("Reco_", mydate)
          #raster::writeRaster(Reco_ras, paste0("Vegetation/output/maps/", mydate, "_NEE.TIF"), format="GTiff", overwrite = TRUE)
          #return(Reco_ras)
        }
        return(list(NEE=NEE_ras, GPP = GPP_ras, Reco = Reco_ras))

      }


      MapHour <-  apply(ihour, 1,NEEhour)
      Myhour <- unlist(MapHour)
      return(Myhour)

    }

    Mapday <-  pbapply::pbapply(iday, 1, NEEday)
    Myday <- unlist(Mapday)
    return(Myday)
}

myNEE <- apply(imput, 1, NeeUHI)
NEE_max <- unlist(myNEE)
NEE_max <- raster::stack(NEE_max)
tmap_mode("plot")
qtm(NEE_max)
tm_shape(NEE_max)+
  tm_raster(n=10, palette = "Spectral", auto.palette.mapping = TRUE,
            title = "CO2 µmol/m2/s", ) +
  tm_layout(scale = 1.5) +
  qtm(study_area, fill = NULL)

getLCZparameters <- function(lcz_map = lcz_map) {
  #lcz.id <- c(seq(1, 10, 1), seq(101, 107))
  lcz <- c(seq(1, 10, 1), seq(11, 17))
  lcz.code <- c(seq(1, 10, 1), "A", "B", "C", "D", "E", "F", "G")
  # lcz.name <- c('compact_high-rise', 'compact_midrise', 'compact_low-rise',
  #               'open_high-rise', 'open_midrise', 'open_low-rise',
  #               'lightweight_low-rise', 'large_low-rise', 'sparsely_built',
  #               'heavy_industry', 'dense_trees', 'scattered_trees', 'bush_scrub',
  #               'low_plants', 'bare_rock_paved', 'bare_soil_sand', 'water')

  lcz.name <- c("Compact highrise", "Compact midrise", "Compact lowrise", "Open highrise",
                "Open midrise", "Open lowrise", "Lightweight low-rise", "Large lowrise",
                "Sparsely built", "Heavy Industry", "Dense trees", "Scattered trees",
                "Bush, scrub", "Low plants", "Bare rock or paved", "Bare soil or sand", "Water")
  lcz.col <- c("#910613", "#D9081C", "#FF0A22", "#C54F1E", "#FF6628", "#FF985E",
               "#FDED3F", "#BBBBBB", "#FFCBAB", "#565656", "#006A18", "#00A926",
               "#628432", "#B5DA7F", "#000000", "#FCF7B1", "#656BFA")

  #LCZ parameters
  SVF.min <- c(0.2, 0.3, 0.2, 0.5, 0.5, 0.6, 0.2, 0.75, 0.85, 0.6, 0.35, 0.5, 0.7, rep(0.9, 4))
  SVF.max <- c(0.4, 0.6, 0.6, 0.7, 0.8, 0.9, 0.5, 0.75, 0.85, 0.9, 0.35, 0.8, 0.9, rep(0.9, 4))
  aspect.ratio.min <- c(3, 0.75, 0.75, 0.75, 0.3, 0.3, 1, 0.1, 0.1, 0.2, 1.5, 0.25, 0.25, rep(0.1, 4))
  aspect.ratio.max <- c(3, 2, 1.5, 1.25, 0.75, 0.75, 2, 0.3, 0.25, 0.5, 1.5, 0.75, 1.0, rep(0.1, 4))
  build.frac.min <- c(40, 40, 40, rep(20,3), 60, 30, 10, 20, rep(9, 7))
  build.frac.max <- c(60, 70, 70, rep(40,3), 90, 50, 20, 30, rep(9, 7))
  imp.frac.min <- c(40, 40, 40, rep(20, 3), 60, 30, 10, 20, rep(0, 7))
  imp.frac.max <- c(60, 70, 70, rep(40, 3), 90, 50, 20, 30, rep(10, 7))
  veg.frac.max <- c(10, 20, 30, 40, 40, 60, 30, 20, 80, 50, rep(100, 4), 10, 100, 100)
  veg.frac.min <- c(0, 0, 0, 30, 20, 30, 0, 0, 60, 40, 90, 90, 90, 90, 0, 90, 90)
  tree.frac.min <- c(rep(0, 10), 90, 90, rep(0, 5))
  tree.frac.max <- c(rep(0, 10), 100, 100, rep(0, 5))
  height.roug.min <- c(26, 10, 3, 26, 10, 3, 2, 3, 3, 5, 3, 3, 2.9, 0.9, 0.24, 0.23,  0)
  height.roug.max <- c(26, 25, 10, 26, 25, 10, 4, 10, 10, 15, 30, 15, 2.9, 0.9, 0.24, 0.23, 0)
  terra.roug.min <- c(8, 6, 6, 7, 5, 5, 4, 5, 5, 5, 8, 5, 4, 3, 1, 1, 1)
  terra.roug.max <- c(8, 7, 6, 8, 6, 6, 5, 5, 6, 6, 8, 6, 5, 4, 2, 2, 1)
  surf.admit.min <- c(1.500, 1.500, 1.200, 1.400, 1.400, 1.200, 800, 1.200, 1.000, 1.000, 0, 1.000, 700, 1.200, 1.200, 600, 1.500)
  surf.admit.max <- c(1.800, 2.000, 1.800, 1.800, 2.000, 1.800, 1.500, 1.800, 1.800, 2.5000, 0, 1.800, 1.500, 1.600, 2.500, 1.400, 1.500)
  surf.albedo.min <- c(rep(0.10, 3), rep(0.12, 3), rep(0.15, 2), rep(0.12, 2), 0.10, rep(0.15, 4), 0.20, 0.02)
  surf.albedo.max <- c(rep(0.20, 3), rep(0.25, 3), 0.35, 0.25, 0.25, 0.20, 0.20, 0.25, 0.30, 0.25, 0.30, 0.35, 0.10)
  antrop.heat.min <- c(50, 74, 74, 49, 24, 24, 34, 49, 9, 310, rep(0, 7))
  antrop.heat.max <- c(300, 74, 74, 49, 24, 24, 34, 49, 9, 310, rep(0, 7))

  # lcz.col <- c('#8c0000', '#d10000', '#ff0100', '#be4d01', '#ff6602', '#ff9955',
  #              '#faee05', '#bcbcbc', '#ffccaa', '#555555', '#006a01', '#01aa00',
  #              '#648526', '#b9db79', '#000000', '#fbf7ae', '#6a6aff')
  lcz.df <- data.frame(lcz, lcz.name, lcz.code, lcz.col, SVF.min, SVF.max, aspect.ratio.min, aspect.ratio.max, build.frac.min, build.frac.max,
                       imp.frac.min, imp.frac.max, veg.frac.min, veg.frac.max, tree.frac.min, tree.frac.max,
                       height.roug.min, height.roug.max, terra.roug.min, terra.roug.max, surf.admit.min, surf.admit.max, surf.albedo.min, surf.albedo.max,
                       antrop.heat.min, antrop.heat.max,
                       stringsAsFactors = F) %>%
    mutate(z0 = ifelse(lcz.code %in% c("G"), 0.0002, #Get z0
                       ifelse(lcz.code %in% c("E", "F"), 0.0005,
                              ifelse(lcz.code=="D", 0.03,
                                     ifelse(lcz.code %in% c(7, "C"), 0.10,
                                            ifelse(lcz.code %in% c(8, "B"), 0.25,
                                                   ifelse(lcz.code %in% c(2, 3, 5, 6, 9, 10), 0.5,
                                                          ifelse(lcz.code %in% c(2, 4), 1.0,
                                                                 ifelse(lcz.code %in% c(1, "A"), 2, ""))))))))) %>%
    mutate(SVF.mean = round((SVF.min + SVF.max)/2, digits = 2),
           aspect.ratio.mean = (aspect.ratio.min + aspect.ratio.max)/2,
           build.frac.mea = (build.frac.min + build.frac.max)/2,
           imp.frac.mean = (imp.frac.min + imp.frac.max)/2,
           veg.frac.mean = (veg.frac.min + veg.frac.max)/2,
           tree.frac.mean = (tree.frac.min +tree.frac.max)/2,
           height.roug.mean = (height.roug.min + height.roug.max)/2,
           terra.roug.mean = (terra.roug.min + terra.roug.max)/2,
           surf.admit.mean = (surf.admit.min + surf.admit.max)/2,
           surf.albedo.mean = (surf.albedo.min + surf.albedo.max)/2,
           antrop.heat.mean = (antrop.heat.min + antrop.heat.max)/2
           )
  #Preprocessing raster
  names(lcz_map) <- "lcz"
  lcz_shp <- terra::as.polygons(rast(lcz_map)) %>% st_as_sf()
  lcz_result <- inner_join(lcz_shp, lcz.df, by="lcz")

  return(lcz_result)
}

lcz_par <- getLCZparameters(lcz_map)

tm_shape(lcz_par) +
  tm_polygons("aspect.ratio.mean")



