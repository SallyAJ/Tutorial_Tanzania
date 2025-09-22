# Climate Data Tutorial

In this tutorial, we will primarily use data extracted from the Climate Data Store of the Copernicus Climate Change Service, which has undergone some preliminary processing steps.  
Furthermore, we will use bias-corrected and downscaled climate projections from CMIP6 models (based on DBCCA, please refer to Jahn et al. 2025, in preparation).  
All data are provided with a **monthly temporal resolution** and a **spatial resolution of 0.1°**.

## Data Availability

Data is available at: 
https://drive.google.com/drive/folders/15OHTdY9qsx4TqfRQu6yu4G6Y5BmP9vo6?usp=sharing

## Getting used to the (new) data format(s)

1. **NetCDF format (.nc)**  
   - A common data format used for weather and climate data.  
   - Most weather and climate datasets are published primarily or additionally in the NetCDF format.  
   - Further information: [https://climatedataguide.ucar.edu/climate-tools/NetCDF](https://climatedataguide.ucar.edu/climate-tools/NetCDF)

2. **Shapefile (.shp) format**  
   - A common format for storing vector-based geographic data such as points, lines, and polygons.  
   - Many spatial datasets, including administrative boundaries, are published in the Shapefile format.  
   - A nice tool for quick assessments of shapefiles: [https://mapshaper.org/](https://mapshaper.org/)

## Tools and Packages

Please note that there is a wide variety of packages available for working with NetCDF data and shapefiles: In **R**: `raster`, `terra`, or `stars` (and in combination with `sf`) can be used to process these datasets. 
For more advanced users, the `ncdf4` library is also a good alternative. For long-term projects, switching to **Python** and using packages like `xarray` is recommended.
In this tutorial, we will primarily rely on the library `terra`, as it is likely the package most users are familiar with. This choice may not always provide the best performance or programming flexibility, 
but it serves well for demonstrating key concepts. The aim of this tutorial is hence to provide a practical introduction to working with spatial data, especially highlighting the types of analysis that can be / should be performed.