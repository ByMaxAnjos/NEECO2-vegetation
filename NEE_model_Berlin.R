#=======================================================================================
#
# Title:       Net Ecosystem Exchange of CO2.
# Author:      Dr. Max Anjos (maxanjos@campus.ul.pt)
# Description: More details on the approach are available at:
#             https://github.com/ByMaxAnjos/CO2-vegetation-emissions
# Data: 28.06.2023
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
library(rdwd)
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

# Get study area polygon from OpenStreetMap data
city <- "Berlin"
mycrs <- "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"

shp_verify <- osmdata::getbb(city, format_out = "sf_polygon", limit = 1, featuretype = "city")
# Check if polygon was obtained successfully
if(!is.null(shp_verify$geometry) & !inherits(shp_verify, "list")) {
  study_area <- shp_verify$geometry
  study_area <- st_make_valid(study_area) %>%
    st_as_sf() %>% 
    st_transform(crs = mycrs)
} else {
  study_area <- shp_verify$multipolygon
  study_area <- st_make_valid(study_area) %>%
    st_as_sf() %>%
    st_transform(crs= mycrs)
}

qtm(study_area)# Plot map

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
## Soil temperature
#=================================================================
#Get and read climate data from DW-------------------------
# station_id <- as.data.frame(findID(city, exactmatch=FALSE)) |> # List of all stations within the city
#   magrittr::set_colnames(c('id')) |>
#   tibble::rownames_to_column("names")
# print(station_id)
# selected <- 5 #selected station
# metadata <- metaInfo(station_id$id[selected]) # Get metadata of station
# 
# data <- selectDWD(id=station_id$id[selected], res="hourly", var="soil_temperature", per="historical", current=TRUE) # since 2019
# file <- dataDWD(data, read=FALSE)
# soil_temperature <- readDWD(file, varnames=TRUE) %>%
#   dplyr::select(MESS_DATUM, V_TE002.Erdbodentemperatur_002cm, V_TE005.Erdbodentemperatur_005cm,V_TE010.Erdbodentemperatur_010cm, 
#                 V_TE020.Erdbodentemperatur_020cm, V_TE050.Erdbodentemperatur_050cm, V_TE100.Erdbodentemperatur_100cm) %>%
#   rename(date= MESS_DATUM, soil_temp_002cm = V_TE002.Erdbodentemperatur_002cm, soil_temp_005cm = V_TE005.Erdbodentemperatur_005cm,
#          soil_temp_010cm = V_TE010.Erdbodentemperatur_010cm, soil_temp_020cm = V_TE020.Erdbodentemperatur_020cm,
#          soil_temp_050cm = V_TE050.Erdbodentemperatur_050cm, soil_temp_100cm = V_TE100.Erdbodentemperatur_100cm) %>% 
#   openair::selectByDate(year = 2015:2023)

#soil_temperature <- fread("Data/soil_temperature_ber_2015_2022.csv")

#=================================================================
# Calculate NEE with LCZ air temperature interpolation
#=================================================================
#Define the period
#imonth <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
imonth <- "jul"
iyear <- c(2018)
check_day <- 1
imput <- expand.grid(imonth, iyear)

NeeUHI <- function(imput, roi = study_area, 
                   veg_df = raster_veg, 
                   rg_df = Rg_data, 
                   air_df = air_UCON,
                   lcz = lcz_map,
                   spRes = 500,
                   Kriging = FALSE,
                   IDW = TRUE) {

  imonth <- imput[1]
  iyear <- imput[2]

  #Select Rg data
  if(is.null(check_day)) {
    rg_model <- rg_df %>%
      selectByDate(year = iyear, month = imonth) %>%
      mutate(Rg = ifelse(Rg < 4, 0, Rg))
  } else {
    rg_model <- rg_df %>%
      selectByDate(year = iyear, month = imonth, day = check_day) %>%
      mutate(Rg = ifelse(Rg < 4, 0, Rg))
  }

  #Select air data
  if(is.null(check_day)) {
    air_model <- air_df %>%
      openair::selectByDate(year = iyear, month = imonth)
  } else {
    air_model <- air_df %>%
      openair::selectByDate(year = iyear, month = imonth, day = check_day)
  }

  meteo_model <- left_join(air_model, rg_model, by="date") 

  veg_model <- terra::as.data.frame(veg_df, xy=TRUE, na.rm = TRUE) %>% set_names(c("x", "y", "veg_fraction")) %>%
    #filter(veg_fraction >10) %>%
    mutate(FID=row_number())

    # resample iLCZ
    # Convert to spatial pixel
    grid_stations <- raster(raster::extent(roi), res = spRes)
    raster::values(grid_stations) <- 1:raster::ncell(grid_stations)
    raster::crs(grid_stations) <- mycrs

    iLCZ <- raster::resample(lcz, grid_stations)
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

        if(Kriging == TRUE) {
          
          ## [using ordinary kriging]
          krige_vgm <- autofitVariogram(airT ~ lcz, as(train_mod, "Spatial"))$var_model
          krige_mod = gstat(formula = airT ~ lcz, model = krige_vgm$var_model, data = train_mod)
          krige_map = predict(krige_mod, newdata=lcz_stars, debug.level = 1)
          krige_map = krige_map["var1.pred",,]
          krige_map = raster::raster(rast(krige_map))
          krige_map = raster::crop(krige_map, raster::extent(roi))
          krige_map = raster::mask(krige_map, roi)
          mydate <- rad_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
          mydate <- gsub("[: -]", "", mydate$date, perl=TRUE)
          names(krige_map) <- paste0("Krige_LCZ_", mydate)
          
          # Resample air map
          air_resample = raster::resample(krige_map, veg_df)
          #raster::writeRaster(air_resample, paste0("Building/outputs/maps/UHI/", mydate, "_UHI.TIF"), format="GTiff", overwrite = TRUE)
          #Merge building with air raster
          air_poi <- as_tibble(rasterToPoints(air_resample)) %>% set_names(c("x", "y", "airT"))
          #veg_poi <- as_tibble(rasterToPoints(veg_df)) %>% set_names(c("x", "y", "veg_fraction"))
          veg_cal <- inner_join(fcover_model %>% dplyr::select(x, y, Rg, veg_fraction), air_poi, by= c("x", "y")) %>%
            mutate(hour = paste0(myhour))
          
        }
        
        if(IDW==TRUE){
          
          #IDW
          idw_mod = gstat(formula = airT ~ 1, data = train_mod)
          idw_map = predict(idw_mod, lcz_stars)
          idw_map = idw_map["var1.pred",,]
          idw_map <- raster(rast(idw_map))
          mydate <- rad_model %>% distinct(date, .keep_all = FALSE) %>% as_tibble()
          mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
          names(idw_map) <- paste0("IDW_airT_",mydate)
          
          # Resample air map
          air_resample = raster::resample(idw_map, veg_df)
          #raster::writeRaster(air_resample, paste0("Building/outputs/maps/UHI/", mydate, "_UHI.TIF"), format="GTiff", overwrite = TRUE)
          #Merge building with air raster
          air_poi <- as_tibble(rasterToPoints(air_resample)) %>% set_names(c("x", "y", "airT"))
          #veg_poi <- as_tibble(rasterToPoints(veg_df)) %>% set_names(c("x", "y", "veg_fraction"))
          veg_cal <- inner_join(fcover_model %>% dplyr::select(x, y, Rg, veg_fraction), air_poi, by= c("x", "y")) %>%
            mutate(hour = paste0(myhour))
        }
        
        #Prepare the function
        ##Get Rre
          night_period <- meteo_model %>%
            filter(Rg <4)
          # Calculate seasonality of temperature with rolling window
          temp_seasonal <- zoo::rollapply(data = night_period$airT, width = 7, FUN = mean, by = 4)
          # Define Q10 function for respiration
          q10 <- function(T, Tref, Q10) {Q10 ^ ((T - Tref)/10)}
          # Set parameters for Q10 function
          Tref <- 10  # reference temperature in Celsius
          Q10 <- 2.0  # assumed change in respiration per 10 degrees Celsius increase in temperature
          # Calculate Rref for each time point using Q10 function
          Rref <- q10(temp_seasonal, Tref, Q10)
          Rref <- abs(mean(Rref))
          
        #Calculate the NEE CO2 fluxes
          
        NEEboost <- function(x) {
          x <- x %>% #filter(veg_fraction > 15) %>%
            mutate(B = abs(-8.25 + 0.35*veg_fraction)
                   ,a = 0.0002*veg_fraction + 0.0052 # 0.005 + 0.016*veg_fraction
                   #,yc = if_else(Rg < 20, REddyProc::fLloydTaylor(10, 330, airT + 273.15)[1], 0)
                   ,Reco = ifelse(Rg <4, Rref * (exp(230 * (1/(10 - -46.02) - 1/(airT - -46.02))))*(veg_fraction/100), 0)
                   #,Reco = 1.3 * (exp(230 * (1/(15 - -46.02) - 1/(airT - -46.02))))*(veg_fraction/100)
                   #,yc = 1.3 * (exp(230 * (1/(15 - -46.02) - 1/(airT - -46.02))))*(veg_fraction/100)
                   ,o = 0.96
                   ,PAR = 0.505*Rg
                   ,NEE = Reco - 1/2*o*(a*PAR+B-((a*PAR+B)^2 -4*a*B*o*PAR)^0.5)
                   ,GPP = -1/2*o*(a*PAR+B-((a*PAR+B)^2 -4*a*B*o*PAR)^0.5)
                   #,GPP = ifelse(GPP1>0, 0, GPP1)
                   #,NEE = GPP-Reco
                   )
          return(x)
        }

        NEE_data <- NEEboost(veg_cal)

          #Get NEE
          NEE_ras <- NEE_data %>% dplyr::select(x, y, NEE)
          NEE_ras = raster::rasterFromXYZ(xyz = NEE_ras,crs = mycrs)
          NEE_ras = raster::crop(NEE_ras, raster::extent(roi))
          NEE_ras= raster::mask(NEE_ras, roi)
          mydate <- rad_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
          mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
          names(NEE_ras) <- paste0("NEE_", mydate)
          #raster::writeRaster(NEE_ras, paste0("Vegetation/output/maps/", mydate, "_NEE.TIF"), format="GTiff", overwrite = TRUE)
          return(NEE_ras)
          
        # if(GPP == TRUE) {
        #   #Get GPP
        #   GPP_ras <- NEE_data %>% dplyr::select(x, y, GPP)
        #   GPP_ras = raster::rasterFromXYZ(xyz = GPP_ras,crs = mycrs)
        #   GPP_ras = raster::crop(GPP_ras, raster::extent(roi))
        #   GPP_ras= raster::mask(GPP_ras, roi)
        #   mydate <- rad_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
        #   mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
        #   names(GPP_ras) <- paste0("GPP_", mydate)
        #   #raster::writeRaster(GPP_ras, paste0("Vegetation/output/maps/", mydate, "_NEE.TIF"), format="GTiff", overwrite = TRUE)
        #   #return(GPP_ras)
        # }
        # if(Reco ==TRUE) {
        # 
        #   #Get Reco
        #   Reco_ras <- NEE_data %>% dplyr::select(x, y, Reco)
        #   Reco_ras = raster::rasterFromXYZ(xyz = Reco_ras,crs = mycrs)
        #   Reco_ras = raster::crop(Reco_ras, raster::extent(roi))
        #   Reco_ras= raster::mask(Reco_ras, roi)
        #   mydate <- rad_model %>% distinct(date, .keep_all = FALSE) %>% as.data.frame()
        #   mydate <- gsub("[: -]", "" , mydate$date, perl=TRUE)
        #   names(Reco_ras) <- paste0("Reco_", mydate)
        #   #raster::writeRaster(Reco_ras, paste0("Vegetation/output/maps/", mydate, "_NEE.TIF"), format="GTiff", overwrite = TRUE)
        #   #return(Reco_ras)
        # }

      }

      MapHour <-  pbapply::pbapply(ihour, 1,NEEhour)
      Myhour <- unlist(MapHour)
      return(Myhour)

    }

    Mapday <- apply(iday, 1, NEEday)
    Myday <- unlist(Mapday)
    return(Myday)
}

myNEE <- apply(imput, 1, NeeUHI)
NEE_max <- unlist(myNEE)
NEE_max <- raster::stack(NEE_max)

qtm(NEE_max[[10]])
#raster::writeRaster(NEE_max[[21]],"NEE_2018070112.TIF", format="GTiff", overwrite = TRUE)

tmap_mode("plot")
tm_shape(NEE_max[[12]])+
  tm_raster(n=10, palette = "Spectral", 
            title = "CO2 µmol/m2/s") +
  tm_layout(scale = 1.5) +
  qtm(study_area, fill = NULL)



