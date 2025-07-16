import pandas as pd
import numpy as np
import ee
import datetime as dt
from IPython.display import Image

def error(err = ""):
    print("\033[0;39;31m " + err + "\033[0;49;39m")

def warning(warn = ""):
    print("\033[0;32;33m " + warn + "\033[0;49;39m")

def message(mssg = None):
    print("\033[0;32;34m " + mssg + "\033[0;49;39m")

class landsat_ndvi:
    """
    This class gets NDVI from Google Earth Engine
    using WGS84 long-lat bounding boxes.
    - box must be a dictionary with keys: xmin, xmax, ymin, ymax.
    - start_year must be an integer (e.g. 2010).
    - end_year must be an integer (e.g. 2010).
    - name must be a string (e.g. 'Lago Maggiore').
    - show_url must be a boolean.
    """
    def __init__(self, box, start_year, end_year, name, show_url = False):
        self.name = name
        self.dataset = ee.ImageCollection('LANDSAT/COMPOSITES/C02/T1_L2_8DAY_NDVI')
        self.start_year = start_year
        self.end_year = end_year
        self.start_date = '{}-01-01'.format(start_year - 1)
        self.end_date = '{}-12-31'.format(end_year)
        self.box = box
        self.show_url = show_url
        if type(self.box) is not dict:
            error("Error: box is not a dictionary")
        if not all([x in self.box.keys() for x in ['xmin', 'xmax', 'ymin', 'ymax']]):
            error("Error: box should have all keys: xmin, xmax, ymin, ymax")
    
    def ndvi(self):
        self.roi = ee.Geometry.Polygon([
            [self.box['xmin'], self.box['ymin']],
            [self.box['xmax'], self.box['ymin']],
            [self.box['xmax'], self.box['ymax']],
            [self.box['xmin'], self.box['ymax']],
            [self.box['xmin'], self.box['ymin']]
        ])

        # Get the land surface for masking
        self.land = ee.ImageCollection('GLCF/GLS_WATER') \
            .select('water') \
            .min()
        
        self.land = self.land.mask(self.land.eq(1))

        # Load imageCollection
        # select NDVI
        # average through time
        # keep only land pixels (mask)
        self.ndvi = self.dataset \
            .filter(ee.Filter.date(self.start_date, self.end_date)) \
            .select('NDVI') \
            .mean() \
            .mask(self.land)

        # Image Thumbnail
        if self.show_url:
            url = self.ndvi.getThumbUrl({
                'min': -1,
                'max': 1,
                'region': self.roi,
                'dimensions': 1000,
                'palette': ['80146E', '83409B', '6C6AB5', '478EC1', '2BABC2', '4AC3BD', '7ED5B8', 'B0E3B8', 'DAEDC2', 'F5F2D8']
            })
            print(url)

        # Reduce image to featureCollection
        self.ndvi = self.ndvi \
            .reduceRegions(**{
              'collection': self.roi,
              'reducer': ee.Reducer.mean(),
              'scale': 30
            })

    def initialize(self):
        self.task_ndvi = ee.batch.Export.table.toDrive(
          collection = self.ndvi,
          selectors = 'mean',
          description = self.name,
          folder = 'andrew-ndvi',
          fileNamePrefix = self.name,
          fileFormat = 'CSV'
        )

    def export(self):
        # Start the tasks in Earth Engine server
        self.task_ndvi.start()

    def check(self):
        # Check status of tasks in Earth Engine server
        print('--- NDVI', self.name, ':', self.start_date, 'to', self.end_date, '--- ')
        print(self.task_ndvi.status()['state'])

    def cancel(self):
        self.task_ndvi.cancel()
