# NDVI & Chlorophyll-a

We used the Normalized Difference Vegetation Index (NDVI) and the concentratin of chlorophyll-a as proxies for the net primary productivity for land and marine environments, respectively.

NDVI was obtained from the composite Landsat product (Landsat Collection 2 Tier 1 Level 2 Annual NDVI Composite) at 30 meters resolution (Masek et al., 2006; Vermote et al. 2016).
We masked this layer using the Landsat Global Land Cover Facility inland surface water layer (Feng et al., 2016) and retained only land pixels.
We calculated the average NDVI for the time period spanning the sampling of the foodwebs plus one year prior for the area inscribed by a circle of radius of 100 meters from the food web location.
For lake food webs, we extended the radius to 2,000 meters, to be sure to capture some land pixels.

Chlorophyll-a concentration was obtained from the Aqua satellite at ca. 4,600 meters resolution (NASA/JPL, 2020).
We calculated the average Chlorophyll-a for the time period spanning food web sampling plus one year prior for the area inscibed by a circle of radius of 5,000 meters.

Average for NDVI and Chlorophyll-a concentration was calculated across time and space.

Landsat Spectral Indices products courtesy of the U.S. Geological Survey Earth Resources Observation and Science Center.
Chlorophill-a Indices products courtesy of NASA/JPL.
Data access and calculations were performed on Google Earth Engine using the Python API (Google, 2025).

# References

Masek, J. G., Vermote, E. F., Saleous, N. E., Wolfe, R., Hall, F. G., Huemmrich, K. F., ... & Lim, T. K. (2006). A Landsat surface reflectance dataset for North America, 1990-2000. IEEE Geoscience and Remote sensing letters, 3(1), 68-72.

Vermote, E., Justice, C., Claverie, M., & Franch, B. (2016). Preliminary analysis of the performance of the Landsat 8/OLI land surface reflectance product. Remote sensing of environment, 185, 46-56.

Feng, M., Sexton, J. O., Channan, S., & Townshend, J. R. (2016). A global, high-resolution (30-m) inland water body dataset for 2000: First results of a topographic–spectral classification algorithm. International Journal of Digital Earth, 9(2), 113-133.

NASA/JPL. (2020). GHRSST Level 2P Global Sea Surface Skin Temperature from the Moderate Resolution Imaging Spectroradiometer (MODIS) on the NASA Aqua satellite (GDS2) [Data set]. NASA Physical Oceanography DAAC.

Google. (2025). earthengine-api (Version 1.5.22). Retrieved from https://anaconda.org/conda-forge/earthengine-api.
