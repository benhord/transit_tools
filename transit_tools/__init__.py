import os
PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

from .version import __version__
from .fetch_lc import *
from .main import *
from .search import *
from .helpers import *
