
# In this tutorial, we will primarily use data extracted from the Climate Data 
# Store of the Copernicus Climate Change Service, which has undergone 
# some preliminary processing steps. 
# Furthermore, we will use bias-corrected and downscaled climate projectios 
# from CMIP6 models (based on DBCCA, please refer to Jahn et al. 2025, 
# in preparation).
# All data are provided with a monthly temporal resolution and 
# a spatial resolution of 0.1°.

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

# R Script: Climate Change and Bioclimatic Variables

# Here, we will now also include the climate change signal by comparing the 
# historical tas from CanESM5 to the future CanESM5 data based on 30-year periods
# (climatology!). 

# We will also inspect bioclimatic variables derived from the climatological 
# monthly tas temperature values in order to generate more biologically 
# meaningful variables. These are often used in species distribution modeling 
# and related ecological modeling techniques. The bioclimatic variables represent
# e.g., annual trends (e.g., mean annual temperature), seasonality (e.g., 
# annual range in temperature) and extreme or limiting environmental factors 
# (e.g., temperature of the coldest and warmest month). A quarter is a period of
# three months (1/4 of the year).

# For other core climate extreme indices, see the ETCCDI Climdex indices:
# https://etccdi.pacificclimate.org/list_27_indices.shtml
# Respective R packages are available to calculate these indices, making it 
# straightforward to analyze temperature and precipitation extremes using 
# standard ETCCDI definitions.

path <- "Z:/Tanzania" # your data path

# At first, please install and load libraries.
library(terra)
library(dismo)
library(viridis)

################################################################################
# 1. Get the data. 
# 1) Read in all downloaded data.

# observational period (tas, tasmin, tasmax, pr)
path_file_tas_obs <- file.path(path, 
                       "Africa_tas_ERA5_Land_observation_ref_1985_2014_box_monthly.nc")
tas_obs <- terra::rast(path_file_tas_obs)

path_file_tasmin_obs <- file.path(path, 
                               "Africa_tasmin_ERA5_Land_observation_ref_1985_2014_box_monthly.nc")
tasmin_obs <- terra::rast(path_file_tasmin_obs)

path_file_tasmax_obs <- file.path(path, 
                                  "Africa_tasmax_ERA5_Land_observation_ref_1985_2014_box_monthly.nc")
tasmax_obs <- terra::rast(path_file_tasmax_obs)

path_file_tp_obs <- file.path(path, 
                                  "Africa_pr_CHIRPS_observation_ref_1985_2014_box_monthly.nc")
pr_obs <- terra::rast(path_file_tp_obs)

# CanESM5 (bias-corrected and downscaled projections) (tas)
path_file_tas_hist_gcm <- file.path(path, 
                                  "Africa_tas_DBCCA_CanESM5_1985_2014_r1i1p1f1_historical_box_monthly.nc")
tas_hist_gcm <- terra::rast(path_file_tas_hist_gcm)

path_file_tas_scen_gcm <- file.path(path, 
                                    "Africa_tas_DBCCA_CanESM5_2015_2100_r1i1p1f1_ssp245_box_monthly.nc")
tas_scen_gcm <- terra::rast(path_file_tas_scen_gcm)

# now, we should do all the data sanity checks....

# 2)  Check tas climate change signal Tanzania
# Generate Tanzania files
admin_sel <- 0
admin_layer <- paste0("ADM_", admin_sel)
gadm41_fil <- paste0("gadm41_TZA_", admin_sel, ".shp")

region_shape <- terra::vect(file.path(path, "gadm41_TZA_shp", 
                                      gadm41_fil))  
# check CRS 
if (!crs(tas_scen_gcm, proj=TRUE) == crs(tas_hist_gcm, proj=TRUE)) {
  tas_hist_gcm <- project(tas_hist_gcm, crs(tas_scen_gcm))
}
if (!crs(tas_hist_gcm, proj=TRUE) == crs(region_shape, proj=TRUE)) {
  region_shape <- project(region_shape, crs(tas_hist_gcm))
}
tas_scen_gcm_cropped <- terra::crop(tas_scen_gcm, region_shape)
tas_scen_gcm_clip  <- terra::mask(tas_scen_gcm_cropped, region_shape, touches = TRUE) 
tas_hist_gcm_cropped <- terra::crop(tas_hist_gcm, region_shape)
tas_hist_gcm_clip  <- terra::mask(tas_hist_gcm_cropped, region_shape, touches = TRUE) 

# Create climatological means based on the selected 30-year periods.
# Compare historical (1985-2014) with mid-century future conditions (2041-2070).
# The analysis uses the SSP2-4.5 scenario and CMIP6 climate model CanESM5.

dates_hist <- time(tas_hist_gcm_clip)
years_hist <- format(dates_hist, "%Y")
tas_hist_gcm_yearly_means <- tapp(tas_hist_gcm_clip, years_hist, fun = mean, na.rm = TRUE)
tas_hist_gcm_mean_yearly <- app(tas_hist_gcm_yearly_means, fun = mean, na.rm = TRUE)

dates_scen <- time(tas_scen_gcm_clip)
years_scen <- format(dates_scen, "%Y")
tas_scen_gcm_yearly_means <- tapp(tas_scen_gcm_clip, years_scen, fun = mean, na.rm = TRUE)
layer_years <- as.numeric(sub("X", "", names(tas_scen_gcm_yearly_means)))
subset_idx <- which(layer_years >= 2041 & layer_years <= 2070)
tas_scen_gcm_yearly_means_sub <- tas_scen_gcm_yearly_means[[subset_idx]]
tas_scen_gcm_mean_yearly <- app(tas_scen_gcm_yearly_means_sub, fun = mean, na.rm = TRUE)

signal <- tas_scen_gcm_mean_yearly - tas_hist_gcm_mean_yearly

cols <- viridis(100, option = "inferno") 
plot(signal, main = "Climate Change (2041-2070 vs. 1985-2014)",
     col = cols, range = c(0, 3.5))
plot(region_shape, add = TRUE, border = "red", lwd = 2)

# 2. Calculate the biovars with the dismo package.
# As an example, we will here focus on BIO6 = Min Temperature of Coldest Month 
# and also analyse observational, historical, and future data.
# The package needs 4 mandatory columns (year, ppt, tmin, and tmax), and 12 rows 
# (months) for each year sorted from Jan to Dec.
# Fore more information: https://github.com/rspatial/dismo/blob/master/R/biovars.R

# Notes / Steps to Consider (skipped in this demo)
# Make sure all input data are aligned:
#    - CRS (coordinate reference system)
#    - Spatial extent
#    - Time periods and resolution
#    - etc. 
# Optionally, clip the data to a region of interest (e.g., Tanzania)

# In this demo, we directly calculate climatological monthly means
# without performing the above preprocessing steps.

dates <- time(tasmin_obs)
months <- as.numeric(format(dates, "%m"))
monthly_clim_tasmin <- tapp(tasmin_obs, index = months, fun = mean, na.rm = TRUE)
monthly_clim_tasmax <- tapp(tasmax_obs, index = months, fun = mean, na.rm = TRUE)
monthly_clim_pr <- tapp(pr_obs, index = months, fun = mean, na.rm = TRUE)

# Convert terra raster to raster package RasterBrick (dismo requires this)
tmin_r <- raster::brick(monthly_clim_tasmin) # in Celsius
tmax_r <- raster::brick(monthly_clim_tasmax) # in Celsius
prec_r <- raster::brick(monthly_clim_pr)  # in mm

# Calculate all 19 bioclimatic variables
bio_vars <- dismo::biovars(prec_r, tmin_r, tmax_r)

# Extract BIO6 (Minimum Temperature of Coldest Month)
bio6 <- bio_vars[[6]]  # BIO6 is the 6th layer

# Plot result
plot(bio6, main = "BIO6: Min Temperature of Coldest Month (°C)", col = cols)
plot(region_shape, add = TRUE, border = "red", lwd = 2)

################################################################################
# 3. Own Exercise: 
# 3.1 Calculate and plot climate change signal for end-century conditions based
#     on SSP2-4.5 and CMIP6 model CanESM5 (e.g., period 2071-2100)
# 3.2 Evaluate further bioclimatic variables (e.g., BIO18: Precipitation of 
#     Warmest Quarter) in the observational period. 
# 3.3 ....
