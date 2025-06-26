import numpy as np

def round_sf(x, sig = 3):
    if x == 0:
        return 0
    return np.round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)