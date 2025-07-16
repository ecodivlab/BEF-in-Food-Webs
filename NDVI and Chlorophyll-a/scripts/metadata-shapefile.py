import pandas as pd
import geopandas as gpd
from os.path import join as path

d = pd.read_csv(path('data', 'Master dataset_FuSED.csv'))
d = d[[
  'lat',
  'lon',
  'study_ID',
  'FW_name',
  'ecosystem.type',
  'sampling.start.year',
  'sampling.end.year'
]]
p = gpd.GeoDataFrame(
  d,
  geometry = gpd.points_from_xy(d.lon, d.lat),
  crs = "EPSG:4326"
)
p.to_file(path('data', 'foodwebs.shp'))
