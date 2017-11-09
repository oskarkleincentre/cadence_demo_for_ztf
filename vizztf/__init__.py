from __future__ import absolute_import
import os
from .version import __version__
from .ztf import *
here = __file__
basedir = os.path.split(here)[0]
example_data = os.path.join(basedir, 'example_data')

