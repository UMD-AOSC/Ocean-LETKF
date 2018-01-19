
import os
from struct import unpack
import numpy as np
from coamps_grid import COAMPSGrid

class ncodaField3D:
  """ structure to hold NFLUX VAR 3D field"""
  def __init__(self):
    self.m = 0
    self.n = 0
    self.l = 0
    self.nest = 0

  def grid(self,COAMPSGrid,nest):
    self.nest = nest
    self.l = COAMPSGrid.kko
    self.m = COAMPSGrid.m[nest-1]
    self.n = COAMPSGrid.n[nest-1]

  def read(self,filename):
      f = open(filename, "rb")  
      data = np.fromfile(f, dtype=">f4", count=-1)
      self.data = np.reshape(data, (self.m, self.n, self.l), order='F')
      f.close()

class ncodaField2D:
  """ structure to hold NFLUX VAR 2D field"""
  def __init__(self):
    self.m = 0
    self.n = 0
    self.nest = 0
    self.Grid = COAMPSGrid(' ')

  def grid(self,COAMPSGrid,nest):
    self.nest = nest
    self.m = COAMPSGrid.m[nest-1]
    self.n = COAMPSGrid.n[nest-1]

  def read(self,filename):
      f = open(filename, "rb")  
      data = np.fromfile(f, dtype=">f4", count=-1)
      self.data = np.reshape(data, (self.m, self.n), order='F')
      f.close()

