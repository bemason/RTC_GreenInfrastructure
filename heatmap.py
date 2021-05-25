from numpy import genfromtxt
from numpy import array
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

# Captured Heatmap - Traditional
z1 = array([[-42.04, -21.36, -0.677, 20.00, 40.68], 
            [-162.0, -120.7, -79.40, -38.08, 3.229],
            [-243.6, -179.3, -114.9, -50.59, 13.76]])  

# Captured Heatmap - Traditional with RTC
z2 = array([[8.842, 29.52, 50.21, 70.89, 91.57], 
            [-0.140, 41.18, 82.49, 123.8, 165.1],
            [-55.69, 8.665, 73.02, 137.4, 201.7]]) 

# Captured Heatmap - Amended
z3 = array([[9.076, 29.76, 50.44, 71.11, 91.79], 
            [17.13, 58.45, 99.76, 141.1, 182.4],
            [26.65, 91.00, 155.4, 219.7, 284.1]])

# SCS 6HR Design Storm (mm)
rows1 = array([12.7, 12.7, 12.7, 12.7, 12.7,
               25.4, 25.4, 25.4, 25.4, 25.4,  
               50.8, 50.8, 50.8, 50.8, 50.8])

# TP Influent Concentration (mg/L)
cols1 = array([0.2, 0.6, 1.0, 1.4, 1.8,
               0.2, 0.6, 1.0, 1.4, 1.8,
               0.2, 0.6, 1.0, 1.4, 1.8])

xi1 = array([12.7, 25.4, 50.8], dtype=np.float)
yi1 = array([0.2, 0.6, 1.0, 1.4, 1.8], dtype=np.float)

I21 = interpolate.interp2d(rows1, cols1, z1, kind='linear')
I21 = I21(xi1, yi1)

I22 = interpolate.interp2d(rows1, cols1, z2, kind='linear')
I22 = I22(xi1, yi1)

I23 = interpolate.interp2d(rows1, cols1, z3, kind='linear')
I23 = I23(xi1, yi1)

levels = [0]

fig, ax = plt.subplots(1,3, figsize=(13,5), constrained_layout=True)
ax[0].contourf(xi1, yi1, I21, np.linspace(-300, 300, num=121), cmap='RdBu')
ax[0].contour(xi1, yi1, I21, levels, colors=('k',), linewidths=(2))
#ax[0].set_xlabel('SCS 6HR Design Storm (mm)', size=12)
ax[0].set_ylabel('TP Influent Conc. (mg/L)', size=12)
ax[0].set_title('Non-Amended Cell', size=12)

ax[1].contourf(xi1, yi1, I22, np.linspace(-300, 300, num=121), cmap='RdBu')
CS1 = ax[1].contour(xi1, yi1, I22, levels, colors=('k',), linewidths=(2))
ax[1].set_xlabel('SCS 6HR Design Storm (mm)', size=12)
#ax[1].set_ylabel('TP Influent Conc. (mg/L)', size=12) 
ax[1].set_title('Non-Amended Cell with RTC', size=12)

CS2 = ax[2].contourf(xi1, yi1, I23, np.linspace(-300, 300, num=121), cmap='RdBu')
#ax[2].contour(xi1, yi1, I23, levels, colors=('k',), linewidths=(2))
#ax[2].set_xlabel('SCS 6HR Design Storm (mm)', size=12)
#ax[2].set_ylabel('TP Influent Conc. (mg/L)', size=12)
ax[2].set_title('Amended Cell', size=12)
cbar = plt.colorbar(CS2, ticks=range(-300, 300, 50))
cbar.add_lines(CS1)
fig.suptitle('TP Load Captured or Released (g)', size=15)
plt.savefig('heatmap.eps')
plt.show()