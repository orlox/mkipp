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
#mkipp.kipp_plot returns an object containing
#   kipp_plot.contour_plot : the return value of matplotlibs contourf. Can be
#                            used to create a colorbar with plt.colorbar() 
#   kipp_plot.histories    : list of history files read. data can be accesed from this
#                            using the get("column_name") function
#   kipp_plot.xlims        : limits of data in x coordinate
kipp_plot = mkipp.kipp_plot(mkipp.Kipp_Args(
        xaxis = "star_age",
        time_units = "Myr",
        identifier = "y",
        log10_on_data = False,
        levels = np.arange(0.0,1.001,0.01),
        decorate_plot = False,
        save_file = False), axis = axis)
bar = plt.colorbar(kipp_plot.contour_plot,pad=0.05)
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
#kipp_data.get_xyz_data returns an object containing
#   xyz_data.xlims : limits of data in x coordinate
#   xyz_data.X     : 2D array of xaxis values of profile data
#   xyz_data.Y     : 2D array of xaxis values of profile data
#   xyz_data.Z     : 2D array of xaxis values of profile data
# the last three can be used as inputs for matplotlib contour or contourf
xyz_data = kipp_data.get_xyz_data(profile_names, kipp_args.xaxis_divide, kipp_args)
#read mixing regions 
#kipp_data.get_mixing_zones returns an object containing
#   mixing_zones.zones     : matplotlib Path objects for each mixing zone.
#                            These can be plotted using add_patch(...)
#   mixing_zones.mix_types : Integer array containing the type of each zone
#   mixing_zones.x_coords  : x coordinates for points at the surface
#   mixing_zones.y_coords  : y coordinates for points at the surface
#   mixing_zones.histories : mesa_data history files to access additional data
# the last three can be used as inputs for matplotlib contour or contourf
mixing_zones = kipp_data.get_mixing_zones(history_names, kipp_args.xaxis_divide, xyz_data.xlims, kipp_args)
# just plot convection, overshooting and semiconvection
for i,zone in enumerate(mixing_zones.zones):
    color = ""
    #Convective mixing
    if mixing_zones.mix_types[i] == 1: #convection
        color = '#332288'
    #Overshooting 
    elif mixing_zones.mix_types[i] == 3: #overshooting
        color = '#117733'
    #Semiconvective mixing
    elif mixing_zones.mix_types[i] == 4: #semiconvection
        color = '#CC6677'
    else:
        continue
    axis.add_patch(PathPatch(zone, color=color, alpha = 0.5, lw = 0))
CS = plt.contour(xyz_data.X, xyz_data.Y, xyz_data.Z, [0,4,8], colors='k')
plt.clabel(CS, inline=1, fontsize=10)
axis.plot(mixing_zones.x_coords, mixing_zones.y_coords, 'k', lw=4)
axis.set_xlabel("model_number")
axis.set_ylabel("Mass (solar masses)")
axis.set_xlim(0,max(mixing_zones.x_coords))
axis.set_ylim(0,max(mixing_zones.y_coords))
plt.savefig("Kippenhahn3.png")
