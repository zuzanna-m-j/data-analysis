#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,41,41)
X, Y = np.meshgrid(x, x)
L = 5
w = 2*np.pi/L

print(f"L = {L}")
print(f"Omega = {w}, {w * (2*np.pi)}")

grating = np.sin(w * X)
plt.set_cmap("gray")
plt.subplot(121)
plt.imshow(grating)

ft = np.fft.ifftshift(grating)
ft = np.fft.fft2(ft)
ft = np.fft.fftshift(ft)
ft = abs(ft)

plt.subplot(122)
plt.imshow(ft)
plt.xlim([369//2 - 20, 369//2 + 20])  # Note, order is reversed for y
plt.ylim([369//2 + 20, 369//2 - 20])  # Note, order is reversed for y
print(f"{369//2 + 1}")
plt.show()

X = np.zeros_like(X)
Y = np.zeros_like(X)
Z = np.zeros_like(X)

for r in range(len(ft)):
    for c in range(len(ft[r])):
        X[r,c] = c
        Y[r,c] = r
        Z[r,c] = ft[r,c]

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()

# mask = np.argwhere((ft.real) != 0)
# ft = abs(ft)
# print(mask)
# for c in mask:
#     ar = ft[c]
#     print(f"{ar[0] - (369//2 + 1)}, {ar[1] - (369//2 + 1)}")


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





