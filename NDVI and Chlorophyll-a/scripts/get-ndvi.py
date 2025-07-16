from os.path import join as path
import ee
import geopandas as gpd
import datetime as dt
from scripts.ndvi import landsat_ndvi

# ee.Authenticate()                       # authenticate on new machine
ee.Initialize(project = 'ee-emilioberti') # initialize module

# load park shapefile and define boundaries
location = gpd.read_file(path('data', 'foodwebs.shp'))
location = location[location.study_ID != 'Intertidal rockpools']

# check that all food webs have sampling years
check_start = [x for x in location['sampling.s']]
if all([isinstance(x, int) for x in check_start]):
  print(" - Start years all good")
else:
  print(" - Problems with start years")

check_end = [x for x in location['sampling.e']]
if all([isinstance(x, int) for x in check_end]):
  print(" - End years all good")
else:
  print(" - Problems with end years")

# start NDVI tasks for all foodwebs
# i = 148 # HEW05
# i = 317 # Maggiore
for i in [x for x in location.index]:

  print(' - ', location.loc[i].FW_name)
  
  if 'lake' in location.loc[i].study_ID or 'Lake' in location.loc[i].study_ID:
    roi = location \
      .to_crs(3395) \
      .buffer(2000) \
      .to_crs(4326) \
      .loc[i] \
      .bounds
  else:
    roi = location \
      .to_crs(3395) \
      .buffer(100) \
      .to_crs(4326) \
      .loc[i] \
      .bounds

  bbox = {
    'xmin': roi[0],
    'xmax': roi[2],
    'ymin': roi[1],
    'ymax': roi[3]
  }
  
  start_year = int(location.loc[i]['sampling.s'])
  end_year = int(location.loc[i]['sampling.e'])

  task = landsat_ndvi(
    bbox,
    start_year,
    end_year,
    location.loc[i].FW_name,
    show_thumbnail
  )
  task.ndvi()
  task.initialize()
  task.export()
  # task.check()
