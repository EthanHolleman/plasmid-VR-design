import pytest

import os
import sys
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agg_seq_metrics import *

# class Snakemake():
#     def __init__(self **kwargs):
#         self.__dict__.update(kwargs)

# @pytest.fixture
# def snakemake():  # simulate snakemake object