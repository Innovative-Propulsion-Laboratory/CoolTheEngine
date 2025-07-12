from itertools import product
import numpy as np

# 1) Define your parameters in an ordered mapping
params = {
    "use_TEOS_PDMS": [True, False],
    "channel_width_inj": [0.001, 0.002, 0.003],
    "channel_width_conv": [0.003, 0.004, 0.004],
    "channel_width_throat": [0.005, 0.006, 0.005],
    "channel_width_exit": [0.007, 0.008, 0.009, 0.001],
    "channel_height_inj": [0.010, 0.011, 0.012, 0.013],
    "channel_height_conv": [0.014, 1.0, 2.0],
    "channel_height_throat": [0.015, 0.016],
    "channel_height_exit": [0.017, 0.018],
    "channel_angle_inj": [0.0, 30.0, 45.0],
    "channel_angle_conv": [0.0, 30.0, 45.0],
    "channel_angle_throat": [0.0],
    "channel_angle_exit": [0.0, 30.0, 45.0],
}

# 2) Generate all combinations
names = list(params.keys())
value_lists = [params[n] for n in names]
combinations = list(product(*value_lists))

# 3a) As a NumPy array (dtype=object since mixed types)
arr = np.array(combinations, dtype=object)
print(arr.shape)
