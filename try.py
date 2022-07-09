#!/usr/bin/env python3
import numpy as np
data=np.array([1,2,3,4,5,6])
k = 2
print(data)
print(data.reshape(-1, k))
data1=np.mean(data.reshape(-1, k), axis=1)
print(data1)