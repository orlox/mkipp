#!/usr/bin/env python
import mkipp
import kipp_data
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
import numpy as np

#simple example (saves to default filename Kippenhahn.png)
mkipp.kipp_plot(mkipp.Kipp_Args())

#plot of Helium abundance against time, independent decoration
fig = plt.figure()
axis = plt.gca()
contour_plot, histories, xlims = mkipp.kipp_plot(mkipp.Kipp_Args(
        xaxis = "star_age",
        time_units = "Myr",
        identifier = "y",
        log10_on_data = False,
        levels = np.arange(0.0,1.001,0.01),
        decorate_plot = False,
        save_file = False), axis = axis)
bar = plt.colorbar(contour_plot,pad=0.05)
bar.set_label("Helium abundance")
axis.set_xlabel("Time (Myr)")
axis.set_ylabel("Mass (solar masses)")
plt.savefig("Kippenhahn2.png")

#Reading out mixing regions and data, and plotting independently
kipp_args = mkipp.Kipp_Args()
fig = plt.figure()
axis = plt.gca()
profile_names = ["LOGS/profile"+str(int(i))+".data" \
        for i in np.loadtxt("LOGS/profiles.index", skiprows = 1, usecols = (2,))]
#if data is distributed among several history.data files, you can provide them
history_names = ["LOGS/history.data"]
#read profile data
min_x_coord, max_x_coord, X, Y, Z = \
        kipp_data.get_xyz_data(profile_names, kipp_args.xaxis_divide, kipp_args)
#read mixing regions 
zones, mix_types, x_coords, y_coords, histories = kipp_data.get_mixing_zones(\
        history_names, kipp_args.xaxis_divide, min_x_coord, max_x_coord, kipp_args)
# just plot convection, overshooting and semiconvection
for i,zone in enumerate(zones):
    color = ""
    #Convective mixing
    if mix_types[i] == 1: #convection
        color = '#332288'
    #Overshooting 
    elif mix_types[i] == 3: #overshooting
        color = '#117733'
    #Semiconvective mixing
    elif mix_types[i] == 4: #semiconvection
        color = '#CC6677'
    else:
        continue
    axis.add_patch(PathPatch(zone, color=color, alpha = 0.5, lw = 0))
CS = plt.contour(X, Y, Z, [0,4,8], colors='k')
plt.clabel(CS, inline=1, fontsize=10)
axis.plot(x_coords,y_coords,'k',lw=4)
axis.set_xlabel("model_number")
axis.set_ylabel("Mass (solar masses)")
axis.set_xlim(0,max(x_coords))
axis.set_ylim(0,max(y_coords))
plt.savefig("Kippenhahn3.png")
