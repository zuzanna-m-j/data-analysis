#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def fit(x,cs,cp,squiggle,l0):
    return cs + (cp - cs) * 1/2 * (np.tanh(x/squiggle - l0) + 1)
cs = 0.5
cp = 0.9
squiggle = 2
l0 = 1
x = np.arange(-40,40,1)
ydata = cs + (cp - cs) * 1/2 * (np.tanh(x/squiggle - l0) + 1)
popt, pcov = curve_fit(fit, x, ydata, p0 = [0.5,0.8, 2, 1])
plt.plot(x,ydata,'o')
plt.plot(x,fit(x, *popt))
plt.show()
print(popt)