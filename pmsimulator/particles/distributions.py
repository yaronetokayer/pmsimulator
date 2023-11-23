import numpy as np

ALLOWED_DISTRIBUTIONS = ['uniform', 'normal', 'allones']

def initialize_uniform(n, bounds=(0, 1)):
    return np.random.uniform(low=bounds[0], high=bounds[1], size=n)

def initialize_normal(n, mean=0, std=1):
    return np.random.normal(mean, std, n)

def initialize_all_ones(n, scale=1.0):
    return scale * np.ones(n)