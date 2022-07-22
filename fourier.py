#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


# gratings.py

x = np.arange(-10, 11.5, 0.5)

X, Y = np.meshgrid(x, x)


angle = 0
l = 5
w = 2 * np.pi/l
grating = np.sin(w * X )
print(w * np.pi)

plt.set_cmap("gray")

plt.subplot(121)
plt.imshow(grating)

# Calculate Fourier transform of grating
ft = np.fft.ifftshift(grating)
ft = np.fft.fft2(ft)
ft = np.fft.fftshift(ft)

plt.subplot(122)
plt.imshow(abs(ft))
# plt.xlim([480, 520])
# plt.ylim([520, 480])  # Note, order is reversed for y
plt.show()

# from skimage.draw import rectangle, disk
# img = np.zeros((8, 10), dtype=np.uint8)
# start = (1, 4)
# extent = (2, 2)
# rr, cc = rectangle(start, extent=extent, shape=img.shape)
# img[rr, cc] = 11
# print(img)
# plt.imshow(img, cmap = 'binary_r')
# plt.show()
# shape = (4, 4)
# img = np.zeros(shape, dtype=np.uint8)
# rr, cc = disk((0, 0), 2, shape=shape)
# img[rr, cc] = 1
# img = np.zeros(shape, dtype=np.uint8)
# # Negative coordinates in rr and cc perform a wraparound
# rr, cc = disk((0, 0), 2, shape=None)
# img[rr, cc] = 1
# img = np.zeros((10, 10), dtype=np.uint8)
# rr, cc = disk((4, 4), 5)
# img[rr, cc] = 1





