import numpy as np

import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

from yarf.Devices import *
from yarf.Utils import tr_logger as logger

class Transient():

    def __init__(self, name):
        self.name = name


