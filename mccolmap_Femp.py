#!/usr/bin/env python3
# program mccolmap.py 

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter
from scipy.ndimage import median_filter

# read table
# columns:
# 0:F  1:M  2:rate1  3:rate2  4:Femp  5:Pemp  6:Pstd
data = np.loadtxt("mcsim_results_weighted_6_regular.dat")
data = np.loadtxt("mcsim_results.dat")

# choose plotted column (Python numbering)
j = 4     # column 5 as indexing starts with 0

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
sigma_smooth = 2.0   # in grid-cell units; try 0.5, 1.0, 1.5
Zplot = gaussian_filter(Zgrid, sigma=sigma_smooth, mode="nearest")
#Zplot = median_filter(Zgrid, size=17, mode="nearest")

# custom red -> yellow -> green colormap
cmap = LinearSegmentedColormap.from_list(
    "ryg",
    ["red", "yellow", "green"]
)

# plot
plt.ion()	 # set interactive mode, so fig.is redrawn every draw() commanfig = plt.figure(1,figsize=(10,5))
plt.figure(figsize=(7,5))

minlev=0.5
maxlev=3.5
levels = np.arange(minlev, maxlev, 5)
levels_filled = np.linspace(minlev, maxlev, 200)

im = plt.contourf(
    Fvals,
    Mvals,
    #Zgrid,
    Zplot,
    levels=levels_filled,
    cmap=cmap,
    extend="both"
    #levels=200,          # many levels -> smooth appearance
    #cmap=cmap,
    #vmin=70,
    #vmax=100
)

# regular contour lines every 5%
cs = plt.contour(
    Fvals,
    Mvals,
    Zplot,
    levels=levels,
    colors="black",
    linewidths=0.7
)

# emphasize the Femp=1.87 contour
cs34 = plt.contour(
    Fvals,
    Mvals,
    Zplot,
    levels=[2.15],       # [3.28],
    colors="black",
    linewidths=1.0
)

# mark its 2sigma limits
cs34 = plt.contour(
    Fvals,
    Mvals,
    Zplot,
    levels=[2.09,2.21],       # [3.28],
    colors="black",
    linewidths=0.5
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

plt.xlabel("F")
plt.ylabel("M")
plt.title(r'$F_{emp}$')

#cbar = plt.colorbar(im)
cbar = plt.colorbar(
    im,
    ticks=np.arange(minlev,maxlev, 0.5)
)
cbar.set_label(r'$F_{emp}$')

# saturate outside range
im.set_clim(0.5, 3.5)

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

input("Press [enter] to terminate.")

plotfilename="out.png"
#plotfilename = 'algebraic_' + str(tau) + '_' + str(lambdar) + '.png'
plt.savefig(plotfilename)
