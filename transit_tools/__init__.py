import os
PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

from .version import __version__
from .preplc import *
from .main import *
