#!/usr/bin/env python
import mkipp
import kipp_data
import mesa_data
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
import numpy as np

#simple example (saves to default filename Kippenhahn.png)
mkipp.kipp_plot(mkipp.Kipp_Args())

#specify xrange and filename
mkipp.kipp_plot(mkipp.Kipp_Args(save_filename = "Kippenhahn2.png", xlims = [300,600]))

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
axis.set_xlim(kipp_plot.xlims)
plt.savefig("Kippenhahn3.png")

#read out max age of star first, then create a log(tf-t) plot
fig = plt.figure()
axis = plt.gca()
#only need to read star_age column first
history = mesa_data.mesa_data("LOGS/history.data", read_data_cols = ["star_age"])
max_age = max(history.get("star_age"))
kipp_plot = mkipp.kipp_plot(mkipp.Kipp_Args(
        xaxis = "star_age",
        time_units = "yr",
        function_on_xaxis = lambda x: np.log10(max_age+0.01 - x),
        decorate_plot = False,
        save_file = False), axis = axis)
bar = plt.colorbar(kipp_plot.contour_plot,pad=0.05)
bar.set_label("log |eps_nuc|")
axis.set_xlabel("log (tf-t) [yrs]")
axis.set_ylabel("Mass (solar masses)")
plt.savefig("Kippenhahn4.png")

#add custom extractor to compute g
fig = plt.figure()
axis = plt.gca()
msun = 1.99e33
rsun = 6.96e10
cgrav = 6.67e-8
def g_extractor(identifier, log10_on_data, prof, return_data_columns = False):
    if return_data_columns:
        return ["radius", "mass"]
    data = cgrav*prof.get("mass")*msun/(prof.get("radius")*rsun)**2
    if log10_on_data:
        return np.log10(data)
    else:
        return data
kipp_plot = mkipp.kipp_plot(mkipp.Kipp_Args(
        extractor = g_extractor,
        decorate_plot = False,
        save_file = False), axis = axis)
bar = plt.colorbar(kipp_plot.contour_plot,pad=0.05)
bar.set_label("log g")
axis.set_xlabel("model_number")
axis.set_ylabel("Mass (solar masses)")
plt.savefig("Kippenhahn5.png")

#Reading out mixing regions and data, and plotting independently
kipp_args = mkipp.Kipp_Args()
fig = plt.figure()
axis = plt.gca()
profile_paths = mesa_data.get_profile_paths(["LOGS"])
#if data is distributed among several history.data files, you can provide them
history_paths = ["LOGS/history.data"]
#read profile data
#kipp_data.get_xyz_data returns an object containing
#   xyz_data.xlims : limits of data in x coordinate
#   xyz_data.X     : 2D array of xaxis values of profile data
#   xyz_data.Y     : 2D array of xaxis values of profile data
#   xyz_data.Z     : 2D array of xaxis values of profile data
# the last three can be used as inputs for matplotlib contour or contourf
xyz_data = kipp_data.get_xyz_data(profile_paths, kipp_args)
#read mixing regions 
#kipp_data.get_mixing_zones returns an object containing
#   mixing_zones.zones     : matplotlib Path objects for each mixing zone.
#                            These can be plotted using add_patch(...)
#   mixing_zones.mix_types : Integer array containing the type of each zone
#   mixing_zones.x_coords  : x coordinates for points at the surface
#   mixing_zones.y_coords  : y coordinates for points at the surface
#   mixing_zones.histories : mesa_data history files to access additional data
# the last three can be used as inputs for matplotlib contour or contourf
mixing_zones = kipp_data.get_mixing_zones(history_paths, kipp_args, xlims = xyz_data.xlims)
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
if (xyz_data.Z.size > 0):
    CS = plt.contour(xyz_data.X, xyz_data.Y, xyz_data.Z, [0,4,8], colors='k')
    plt.clabel(CS, inline=1, fontsize=10)
axis.plot(mixing_zones.x_coords, mixing_zones.y_coords, 'k', lw=4)
axis.set_xlabel("model_number")
axis.set_ylabel("Mass (solar masses)")
axis.set_xlim(0,max(mixing_zones.x_coords))
axis.set_ylim(0,max(mixing_zones.y_coords))
plt.savefig("Kippenhahn6.png")
