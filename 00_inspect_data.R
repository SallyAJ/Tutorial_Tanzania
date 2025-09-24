
# In this tutorial, we will primarily use data extracted from the Climate Data 
# Store of the Copernicus Climate Change Service, which has undergone 
# some preliminary processing steps. 
# Furthermore, we will use bias-corrected and downscaled climate projectios 
# from CMIP6 models (based on DBCCA, please refer to Jahn et al. 2025, 
# in preparation).
# All data are provided with a monthly temporal resolution and 
# a spatial resolution of 0.1°.

# We use country and administrative unit level data from GADM providing maps and 
# spatial data for all countries and their sub-divisions: https: //gadm.org/index.html


# Data is available at: 
# https://drive.google.com/drive/folders/15OHTdY9qsx4TqfRQu6yu4G6Y5BmP9vo6?usp=sharing

# Getting used to the (new) data format(s):

# 1) the NetCDF format (.nc), a common data format used for weather and climate data. 
# Most weather and climate datasets will be published primarily or additionally
# in the NetCDF format. 
# Further information: https://climatedataguide.ucar.edu/climate-tools/NetCDF

# 2) The Shapefile (.shp) format, a common data format used for storing 
# vector-based geographic data such as points, lines, and polygons. 
# Many spatial datasets, including administrative boundaries, are published 
# in the Shapefile format.
# A nice tool for quick assessments of shapefiles is https://mapshaper.org/ .

# Please note that there is a wide variety of packages available for working 
# with NetCDF data and shapefiles. In R, packages such as raster, terra or stars 
# (and in combination with sf) can be used to process these datasets. For more 
# advanced users, the ncdf4 library is also a good alternative. However, for 
# long-term projects, switching to Python and using packages like xarray is 
# recommended. 

# In this tutorial, we will primarily rely on the library terra, as it is 
# likely the package most users are familiar with. This choice may not always 
# provide the best performance or programming flexibility, but it serves well 
# for demonstrating key concepts. The aim of this tutorial is hence to provide a 
# practical introduction to working with spatial data, especially highlighting 
# the types of analysis that can be / should be performed.
################################################################################

# R Script: Inspect Data and Calculate Climatology

# Here, we will inspect observational near-surface 2m air temperature data 
# (tas) in Tanzania and analyze climatology.

################################################################################

path <- "C:/Users/sjahn/Desktop/Data/Tanzania" # your data path

# At first, please install and load libraries.
library(terra)
library(sp)
library(sf)
library(exactextractr)
library(lubridate)
library(dplyr)

################################################################################
# 1. Get to Know the Data
# 1) Read in the downloaded data.
path_file <- file.path(path, 
                       "Africa_tas_ERA5_Land_observation_ref_1985_2014_box_monthly.nc")
region_raster <- terra::rast(path_file)

# Data Quality Checks: Now that the data is loaded, let's perform some essential checks:
# 1. Longitude range: Is the data on a [-180, 180] longitude grid?
# 2. Coordinate Reference System (CRS): What is its current projection? 
#    Do we need to transform it for our analysis?
# 3. Units: Is the temperature recorded in Celsius?
# 4. Missing values: Are there any missing or invalid entries in the dataset?
# 5. Additional checks: Are the data dimensions, time steps, and spatial extent 
#    consistent with expectations?
# ...

print(region_raster) # familiarize yourself with the data
current_crs <- crs(region_raster, proj=TRUE)
print(current_crs)

# If CRS exists but is not EPSG:4326, reproject
if (!is.na(current_crs) && !grepl("WGS84", as.character(current_crs))) {
  region_raster <- terra:project(region_raster, crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  message("CRS transformed to EPSG:4326 (WGS84).")
# If CRS is undefined but CRS confirmed, you can simply assign EPSG:4326
} else if (is.na(current_crs)) {
  crs(region_raster) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  message("CRS was undefined. Assigned EPSG:4326 (WGS84).")
} else {
# Check passed.
  message("CRS is already EPSG:4326 (WGS84). No action needed.")
}

# etc. 

# create a first simple plot of data (first month)
plot(region_raster[[1]], main = "Temperature (first time step)")

# 2) Read in shapefile and clip 
# Now read in GADM shapefiles for Tanzania.
# The coordinate reference system is longitude/latitude and the WGS84 datum.
admin_sel <- 0
admin_layer <- paste0("ADM_", admin_sel)
gadm41_fil <- paste0("gadm41_TZA_", admin_sel, ".shp")

region_shape <- terra::vect(file.path(path, "gadm41_TZA_shp", 
                              gadm41_fil))  

# check CRS (always!)
if (!crs(region_raster, proj=TRUE) == crs(region_shape, proj=TRUE)) {
  region_shape <- terra::project(region_shape, crs(region_raster))
}

# add to plot
plot(region_shape, add = TRUE, border = "red", lwd = 2)

# To select a specific administrative unit (e.g. level 1), set admin_sel = 1
# and filter by the respective unit name: For example, selecting the region "Mkoa"
# Mkoa <- region_shape[region_shape@data$NAME_1 == "Mkoa", ]

# We can use terra package to create a cropped and masked version of our data.
cropped <- terra::crop(region_raster, region_shape)
r_clip  <- terra::mask(cropped, region_shape, touches = TRUE) # True: all intersecting pixels
# plot area
plot(r_clip[[1]], main = "Temperature (first time step) Tanzania")
plot(region_shape, add = TRUE, border = "red", lwd = 2)

# 3) Let’s plot a time series of the closest grid cell to our location
# Define coordinates
coords <- data.frame(lon = 36.68, lat = -3.37)

ts_arusha <- terra::extract(r_clip, coords)
values <- as.numeric(ts_arusha[ , -1]) # no ID
dates <- time(r_clip)
ts_df <- data.frame(
  date = dates,
  temperature = values
)
ts_df$year_month <- format(ts_df$date, "%Y-%m")

# Plot with monthly x-axis
plot(ts_df$temperature, type = "l",
     xaxt = "n",  
     main = "Temperature Time Series at Arusha", col="red",
     xlab = "Date (Year-Month)", ylab = "Temperature (°C)")
axis(1, at = seq(1, nrow(ts_df), by = 12), labels = ts_df$year_month[seq(1, nrow(ts_df), by = 12)])

################################################################################
# 2. Climate Normals and Anomalies

# 1) Standard reference period and climate normals.

# Anthropogenic activities and natural variations from years to decades shape 
# the Earth’s climate. In order to evaluate anomalous conditions of a specific 
# month or year, the World Meteorological Organization (WMO) defines standard 
# reference periods used to create climate normals. Climate normals can be 
# considered as the typical climate for the period the normals are based on. 
# Until 2020, the most current and widely used standard reference period was 
# the 30-year range of 1981-2010. With the start of 2021, the WMO recommended 
# updating the climate normal reference period to the range 1991-2020. For 
# other assessments, other periods are used, e.g., the IPCC report recently 
# used a 20-year reference period of 1995-2014 in his AR6 report. 

# Here, we base our analysis on the reference period of 1985-2014.
# First, let us calculate the tas climate normal 
# for the reference period based on a yearly time series. 
# For this, we will create for each year average tas. 

years <- format(dates, "%Y")
yearly_means <- tapp(r_clip, years, fun = mean, na.rm = TRUE)

# Please note you might have to subset the data. Here, it is already clipped to
# the period 1985-2014.
ref_mean_yearly <- app(yearly_means, fun = mean, na.rm = TRUE)

# 2) Now let's calculate the anomaly of a specific year with respect to the 
# climate normal. The term anomaly refers to the deviation of a value from the 
# long-term average. Positive or negative anomalies indicate that the average 
# temperatures of a particular year were respectively warmer or cooler than
# the reference value. Here, we now calculate the anomaly of a specific year 
# (e.g., 2010) with respect to the climate normal. 
t2m_2010 <- subset(yearly_means, "X2010")
anom_2010 <- t2m_2010 - ref_mean_yearly

# Plot results
cols <- colorRampPalette(c("blue", "white", "red"))(100)
plot(anom_2010, main = "Temperature Anomaly 2010 (Climatology 1985-2014)",
     col = cols, range = c(-2, 2))
plot(region_shape, add = TRUE, border = "red", lwd = 2)

# 3) Let's repeat parts of step 1) and 2) for calculate the monthly climatology 
# of tas. We will also view the anomalies with respect to the climatology for a 
# particular month (e.g., July). 
months <- as.numeric(format(dates, "%m"))
# climatological months (across all years)
monthly_clim <- tapp(r_clip, index = months, fun = mean, na.rm = TRUE)

anom_month <- r_clip
for (m in 1:12) {
  layer_idx <- which(months == m)               
  anom_month[[layer_idx]] <- r_clip[[layer_idx]] - monthly_clim[[m]]
}

aug_2010_idx <- which(years == 2010 & months == 8)
plot(anom_month[[aug_2010_idx]], 
     main = "Temperature Anomaly August 2010 (Climatology 1985-2014)", 
     col = cols, 
     range = c(-3, 3))  # fix color scale
plot(region_shape, add = TRUE, border = "red", lwd = 2)

# 4) Zonal statistics: Monthly climatology averaged over the entire Tanzanian 
#    region.
 
# To do this we need to average over the latitude and longitude dimensions. 
# A very important consideration however is that we will need to take into 
# account the varying size of the gridded data cells as a function of latitude.
# One way to do this e.g., is to use the cosine of the latitude as a proxy for the 
# varying sizes.

# There exist several packages to spatially average gridded data. Here we use
# R package exactextractr, as it performs faster other packages for many
# real-world applications. It quickly and accurately summarizes raster 
# values over polygonal areas, commonly referred to as zonal statistics.

# Then, we can use the "weighted_mean" argument. One application of this feature 
# is the calculation of zonal statistics on raster data in geographic coordinates.
# Without, the calculation of mean tas e.g., across Tanzanian admin units 
# would assume that each raster cell covered the same area, which is not correct
# for rasters in geographic coordinates (latitude/longitude).
# We can correct for varying cell areas by creating a weighting raster with the 
# area of each cell in the primary raster using the area function. 

region_shape_sf <- sf::st_as_sf(region_shape)  
# compute weighted mean
weighted_mean <- exact_extract(
  monthly_clim, 
  region_shape_sf, 
  fun = 'weighted_mean', weights = cellSize(monthly_clim))
months <- month.abb 
colnames(weighted_mean) <- months
values_weighted_mean <- as.numeric(weighted_mean)

ts_df_weighted_mean <- data.frame(
  date = months,
  temperature = values_weighted_mean
)

plot(ts_df_weighted_mean$temperature, type = "l",
     xaxt = "n", col="red",
     main = "Climatological Annual Cycle Tanzania (1985-2014)",
     xlab = "Month", ylab = "Temperature (°C)",  ylim = c(15, 30))
axis(1, at = seq(1, nrow(ts_df_weighted_mean), by = 1), labels = months)

################################################################################
# 3. Seasonal analysis of tas in Tanzania
# We will now compare seasonal trends in tas. To do this we return to our 
# monthly geographically averaged dataset, and we will downsample the monthly 
# averages to seasonal averages. 
# As an example, we split the data into consecutive three-month periods. 
# We calculate the average of each three-month period.

monthly_tas_tanzania <- exact_extract(
  r_clip, 
  region_shape_sf, 
  fun = 'weighted_mean', weights = cellSize(monthly_clim))

values_tanzania <- as.numeric(monthly_tas_tanzania)
ts_df_tanzania <- data.frame(
  date = dates,
  temperature = values_tanzania
)
ts_df_tanzania$year_month <- format(ts_df_tanzania$date, "%Y-%m")

ts_df_tanzania <- ts_df_tanzania %>%
  mutate(
    month = month(date),
    year = year(date),
    season = case_when(
      month %in% c(12, 1, 2)  ~ "Winter",
      month %in% c(3, 4, 5)   ~ "Spring",
      month %in% c(6, 7, 8)   ~ "Summer",
      month %in% c(9, 10, 11) ~ "Autumn"
    )
  )

seasonal_df <- ts_df_tanzania %>%
  group_by(season, year) %>%
  summarize(t2m = mean(temperature, na.rm = TRUE))

seasons <- c("Winter", "Spring", "Summer", "Autumn")
par(mfrow = c(4,1), mar = c(4,4,2,1))  # 4 rows, 1 column

for (s in seasons) {
  df <- seasonal_df %>% filter(season == s)
  mid_range <- round(mean(df$t2m, na.rm = TRUE))
  plot(df$year, df$t2m, type = "l", col = "red",
       ylim = c(15, 30),
       xlab = "Year", ylab = "Temp (°C)", main = s)
}

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

################################################################################
# 4. Own Exercise: 
# 4.1 Calculate and plot climate normals and anomalies (yearly and monthly) 
#     based on historical CanESM5 data.
# 4.2 Generate monthly climatology (climatological annual cycle) averaged 
#     over the entire Tanzanian region.
# 4.3 Check Model/Obs Consistency: Check e.g., the climatologcial monthly mean  
#     bias between the historical CanESM5 and ERA5-Land data. Plot the bias for 
#     one month (e.g., July).
