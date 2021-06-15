import numpy as np

SEED = 12311997  # turn this into a parameter in snakemake somewhere in future

RAND_GEN =  np.random.default_rng(SEED)