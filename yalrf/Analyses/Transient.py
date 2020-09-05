import numpy as np

import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

from yalrf.Devices import *
from yalrf.Utils import tr_logger as logger

class Transient():

    def __init__(self, name):
        self.name = name


