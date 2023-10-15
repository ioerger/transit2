# What is this?

A simple python tool for hashing objects that python normally treats as unhashable. There are also a few other tools such as ways of deeply hashing functions, and an system for extending this hashability to work with other existing classes, such as dataframes or neural networks.

# How do I use this?

`pip install super_hash`


```python
from super_hash import super_hash, function_hashers, FrozenDict, helpers

normally_unhashable = {
    frozenset({
        frozenset({
            "key-deep-deep": 10
        }.items()): "key-deep",
    }.items()): "first_value",
    "second_value": [
        {"a": 10},
    ]
}
a_hash = super_hash(normally_unhashable)

# 
# extend what can be hashed
# 

# example1:
import pandas as pd
# tell super_hash that pandas dataframes should be converted to csv, then hashed
super_hash.conversion_table[pd.DataFrame] = lambda data_frame : super_hash(data_frame.to_csv())

# example2:
import torch
# create a custom checker function
is_non_scalar_pytorch_tensor = lambda value: isinstance(value, torch.Tensor) and len(value.shape) > 0
# create a custom converter
super_hash.conversion_table[is_non_scalar_pytorch_tensor] = lambda non_scalar_tensor: super_hash(non_scalar_tensor.tolist())

# example3:
class Thing:
    def __super_hash__(self):
        return self.file_path
```
