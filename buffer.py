from itertools import product
import pandas as pd
import json
import numpy as np


json_path = "input/input_bulk.json"
with open(json_path, "r") as f:
    settings = json.load(f)


# 2) Generate all combinations
names = list(settings.keys())
value_lists = [v if isinstance(v, list) else [v] for v in settings.values()]
combinations = list(product(*value_lists))

# 3a) As a NumPy array (dtype=object since mixed types)
df = pd.DataFrame(combinations, columns=names)
print(df)
total = np.prod([len(lst) for lst in value_lists])
print(total)
