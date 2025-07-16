from os.path import join as path
import ee
import geopandas as gpd
import datetime as dt
from scripts.chla import aqua_chla

# ee.Authenticate()
ee.Initialize(project = 'ee-emilioberti')

# load park shapefile and define boundaries
location = gpd.read_file(path('data', 'foodwebs.shp'))
location = location[location.study_ID == 'Intertidal rockpools']

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

buffer = location.to_crs(3395).buffer(5000).to_crs(4326)

# start Chl-a tasks for all foodwebs
for i in [x for x in location.index]:

  print(' - ', location.loc[i].FW_name)
  roi = buffer.loc[i].bounds
  bbox = {
    'xmin': roi[0],
    'xmax': roi[2],
    'ymin': roi[1],
    'ymax': roi[3]
  }
  
  start_year = int(location.loc[i]['sampling.s'])
  end_year = int(location.loc[i]['sampling.e'])

  task = aqua_chla(
    bbox,
    start_year,
    end_year,
    location.loc[i].FW_name
  )
  task.chla()
  task.initialize()
  task.export()
  # task.check()
