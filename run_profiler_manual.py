
%load_ext autoreload
%autoreload 2
import profiler
from utils import config_utils
from datetime import datetime

profiler = profiler.Profiler(lsh_path=config_utils.lsh_path, hashes_path=config_utils.hashes_path, mat_path = config_utils.mat_path)

date_string = "{:%B_%d_%Y_%H_%M}".format(datetime.now())

infams = './'
