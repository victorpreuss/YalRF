import numpy as np

import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

from yarf.Devices import *

import logging

class Transient():

    def __init__(self, name):
        self.name = name


