import ujson as json
import h5py
import multiprocessing as mp
import pandas as pd

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import TermCounts

from utils import config_utils, hashutils
from validation import validation_semantic_similarity

from time import time

# data prep
    print('hogs without annotations {}'.format(len(hogs_without_annotations)))
    print('hogs with annotations {}'.format(len(hogs_with_annotations)))
    return hogs_with_annotations , hogs_without_annotations
