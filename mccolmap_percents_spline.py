#!/usr/bin/env python3
# program mccolmap.py 

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
#from scipy.ndimage import gaussian_filter
#from scipy.ndimage import median_filter
from scipy.interpolate import SmoothBivariateSpline

# read table
# columns:
# 0:F  1:M  2:rate1  3:rate2  4:Femp  5:Pemp  6:Pstd
data = np.loadtxt("mcsim_results_weighted_6_regular.dat")
data = np.loadtxt("mcsim_results.dat")

# choose plotted column (Python numbering)
j = 3     # column 4 as indexing starts with 0

F = data[:,0]
M = data[:,1]
Z = data[:,j]

# unique grid values
Fvals = np.unique(F)
Mvals = np.unique(M)

# reshape to 2D grid
Zgrid = np.full((len(Mvals), len(Fvals)), np.nan)

for i in range(len(data)):
    ix = np.where(Fvals == F[i])[0][0]
    iy = np.where(Mvals == M[i])[0][0]
    Zgrid[iy, ix] = Z[i]

# smooth MC noise on the F-M grid
#sigma_smooth = 0.5   # in grid-cell units; try 0.5, 1.0, 1.5
#Zplot = gaussian_filter(Zgrid, sigma=sigma_smooth, mode="nearest")
#Zplot = median_filter(Zgrid, size=3, mode="nearest")

spline = SmoothBivariateSpline(F, M, Z, s=250)

Ffine = np.linspace(Fvals.min(), Fvals.max(), 300)
Mfine = np.linspace(Mvals.min(), Mvals.max(), 300)
Zfine = spline(Ffine, Mfine)


# custom red -> yellow -> green colormap
cmap = LinearSegmentedColormap.from_list(
    "ryg",
    ["red", "yellow", "green"]
)

# plot
plt.ion()	 # set interactive mode, so fig.is redrawn every draw() commanfig = plt.figure(1,figsize=(10,5))
plt.figure(figsize=(7,5))




levels = np.arange(70, 105, 5)
levels_filled = np.linspace(70, 100, 200)

im = plt.contourf(
    Ffine,
    Mfine,
    Zfine.T,
    levels=np.linspace(70, 100, 200),
    cmap=cmap,
    extend="both"
)

linewidths = [
    2.0 if lev == 90 else 0.7
    for lev in np.arange(70, 105, 5)
]

cs = plt.contour(
    Ffine, Mfine, Zfine.T,
    levels=np.arange(70, 105, 5),
    colors="black",
    linewidths=linewidths
)

# function curve
Fcurve = np.linspace(2.5,7.5, 100)
Mcurve = 70 * Fcurve**(-2.44)
plt.plot(
    Fcurve,
    Mcurve,
    color="blue",
    linewidth=2.5,
    linestyle="--"
)

# mark locus of EMAP:
M0=4.2
for F0 in np.arange(3.1,3.5,0.1) :
   plt.plot(
        F0, M0,
        marker="x",
        markersize=5,
        markeredgecolor="black",
        markerfacecolor="white",
        markeredgewidth=2
    )

# optional contour labels
plt.clabel(cs, fmt="%d", fontsize=8)

plt.xlabel(r"Max/min contrast ratio $\, F$")
plt.ylabel(r"Annual event rate$\, M$")
plt.title(f"Hit rate [%]")

#cbar = plt.colorbar(im)
cbar = plt.colorbar(
    im,
    ticks=np.arange(70, 101, 5)
)
cbar.set_label(r"%")

# saturate outside range
im.set_clim(70, 100)

plt.tight_layout()
#plt.show()


input("Press [enter] to terminate.")

plotfilename="out.png"
#plotfilename = 'algebraic_' + str(tau) + '_' + str(lambdar) + '.png'
plt.savefig(plotfilename)

sys.exit("Bye!")






# plot
plt.figure(figsize=(7,5))

im = plt.imshow(
    Zgrid,
    origin="lower",
    aspect="auto",
    extent=[Fvals[0], Fvals[-1], Mvals[0], Mvals[-1]],
    cmap=cmap
)

plt.xlabel("F")
plt.ylabel("M")
plt.title(f"Heat map of column {j+1}")

cbar = plt.colorbar(im)
cbar.set_label(f"Column {j+1}")

plt.tight_layout()
plt.show()
