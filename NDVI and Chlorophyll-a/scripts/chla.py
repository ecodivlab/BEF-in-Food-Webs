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

class aqua_chla:
    """
    This class gets Chlorophyll-a concetnration from Google Earth Engine
    using WGS84 long-lat bounding boxes.
    - box must be a dictionary with keys: xmin, xmax, ymin, ymax.
    - year must be an integer (e.g. 2010).
    - month must be an integer, 1 to 12 (e.g. 2 for February).
    - clouds must be an integer, 1 to 100, masking out pixels with cloud proportion >= clouds.
    """
    def __init__(self, box, start_year, end_year, name):
        self.name = name
        self.dataset = ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI')
        self.start_year = start_year
        self.end_year = end_year
        self.start_date = '{}-01-01'.format(start_year - 1)
        self.end_date = '{}-12-31'.format(end_year)
        self.box = box
        if type(self.box) is not dict:
            error("Error: box is not a dictionary")
        if not all([x in self.box.keys() for x in ['xmin', 'xmax', 'ymin', 'ymax']]):
            error("Error: box should have all keys: xmin, xmax, ymin, ymax")
    def chla(self):
        self.roi = ee.Geometry.Polygon([
            [self.box['xmin'], self.box['ymin']],
            [self.box['xmax'], self.box['ymin']],
            [self.box['xmax'], self.box['ymax']],
            [self.box['xmin'], self.box['ymax']],
            [self.box['xmin'], self.box['ymin']]
        ])
        self.chl = self.dataset \
            .filter(ee.Filter.date(self.start_date, self.end_date)) \
            .select('chlor_a') \
            .mean() \
            .reduceRegions(**{
              'collection': self.roi,
              'reducer': ee.Reducer.mean(),
              'scale': 4616
            })
    def initialize(self):
        self.task_chl = ee.batch.Export.table.toDrive(
          collection = self.chl,
          selectors = 'mean',
          description = self.name,
          folder = 'andrew-chla',
          fileNamePrefix = self.name,
          fileFormat = 'CSV'
        )
    def export(self):
        # Start the tasks in Earth Engine server
        self.task_chl.start()
    def check(self):
        # Check status of tasks in Earth Engine server
        print('--- NDVI', self.name, ':', self.start_date, 'to', self.end_date, '--- ')
        print(self.task_chl.status()['state'])
    def cancel(self):
        self.task_chl.cancel()
