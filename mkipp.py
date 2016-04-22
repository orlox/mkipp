"""
Kippenhahn plotter based on SimpleMesaKippy.py and Kippy.py
Prehistory: SimpleMesaKippy.py, Author: Alfred Gautschy /23.VIII.2012

            Kippy.py, Author: Pablo Marchant /13.V.2013

            MKippi.py, Author: Matteo Cantiello /V.2013 (Added LOSSES. Color/Line Style inspired by convplot.pro IDL routine of A.Heger) 
                                                /VI.2013 (Added RADIUS, TIME and NORMALIZE Options. OMEGA Allows to plot w-contours)   
                                                /VII.2013 (Tentatively added Rotational diff. coefficients,conv. velocities and Equipartition B-field)

                                                (most of these were removed and replaced with flexible plotting options)

            mkipp.py, Author: Pablo Marchant    /II.2014 (full rewrite, code cleanup and now works as a module, not a script)
                                                /III.2014 (enhancement of mixing regions ploting to deal with holes and merging regions)
                                                /IV.2015 stopped listing enhancements

Requirements: history.data and profiles.data containing 
              History (star_age,model_number,star_mass,photosphere_r,mixing_regions,mix_relr_regions)
              Profile (mass,radius,eps_nuc)
"""
import numpy as np
from math import log10, pi
from mesa_data import *
from kipp_data import *
from collections import namedtuple

#matplotlib specifics
import matplotlib.pyplot as plt
from matplotlib.path import Path

# Create contour levels array (TBD: Need to be improved)
def get_levels_linear(min,max,n_levels):
    delta = (max-min)/n_levels
    levels = np.arange(min-delta,max+delta,delta)
    return levels;
    
def get_levels_log(min,max,n_levels):
    max = int(max)+1
    min = int(round(min,2))-1
    levels = range(0,max+1)
    return levels;

def default_extractor(identifier, log10_on_data, prof, return_data_columns = False):
    if return_data_columns:
        return [identifier]
    if log10_on_data:
        return np.log10(abs(prof.get(identifier))+1e-99)
    else:
        return prof.get(identifier)

#properties of the plotter
class Kipp_Args:
    def __init__(self,
            logs_dirs = ['LOGS'],
            profile_names = [],
            history_names = [],
            clean_data = True,
            extra_history_cols = [],
            identifier = "eps_nuc",
            extractor = default_extractor,
            log10_on_data = True,
            contour_colormap = plt.get_cmap("Blues"),
            levels = [],
            log_levels = True,
            num_levels = 8,
            xaxis = "model_number",
            xaxis_divide = 1.0,
            time_units = "Myr",
            yaxis = "mass",
            yaxis_normalize = False,
            show_conv = True, show_therm = True, show_semi = True, show_over = True, show_rot = False,
            core_masses = ["He","C","O"],
            yresolution = 1000,
            mass_tolerance = 0.0000001,
            radius_tolerance = 0.0000001,
            decorate_plot = True,
            show_plot = False,
            save_file = True,
            save_filename = "Kippenhahn.png"):
        """Initializes properties for a Kippenhahn plot

        Note:
            All arguments are optional, if not provided defaults are assigned

        Args:
            logs_dir (List[str]): List of paths to MESA LOGS directories. If profile_names and
                history_names are not provided, they are automatically generated from logs_dir.
            profile_names (List[str]): List of paths to MESA profile files.
            history_names (List[str]): List of paths to MESA history files.
            clean_data (bool): Clean history files in case of redos. If data has no redos,
                a small performance gain can be had by setting this to False.
            extra_history_cols (List[str]): Additional column names to be extracted from history files.
            identifier (str): String used as identifier of data to plot. If not using any custom
                extractors this is simply the column name in the profile file that will be extracted.
                Default uses eps_nuc.
            extractor : TODO, explain
            log10_on_data (bool): Determines if log10(abs()) is applied to profile data
            contour_colormap (matplotlib.cm): matplotlib color map used to plot contours
            levels (List): List of fixed levels for contour plot (int or float)
            log_levels (bool): if levels is an empty list then they are auto-generated. This
                variable specifies if the date is the log of a quantity or not, in order to
                produce discrete integer levels.
            num_levels (int): Number of automatically generated levels
            xaxis (str): variable for the xaxis, either "model_number" or "star_age"
            xaxis_divide (float): divide xaxis by this value
            time_units (str): When using xaxis = "star_age" this specifies the unit of time
                and sets the value of xaxis_divide. Options are "yr", "Myr" and "Gyr"
            yaxis (str): Quantity plotted in the yaxis. Either "mass" or "radius"
            yaxis_normalize (bool): If True Normalize yaxis at each point using total mass/total radius
            show_conv, show_therm, show_semi, show_over, show_rot (bool): Specifies whether or
                certain mixing regions are displayed.
            core_masses (List(str)): Strings with core masses to plot. Options are "He", "C" and "O".
                Only for yaxis=mass
            yresolution (int): resolution for contour plotting
            mass_tolerance (float): ignore mixing regions smaller than this in solar masses. Ignored
                if yaxis="radius"
            radius_tolerance (float): ignore mixing regions smaller than this in solar radii. Ignored
                if yaxis="mass"
            decorate_plot (bool): If True, then axis labels are included.
            show_plot (bool): If True, pyplot.show() is ran at the end
            save_file (bool): If True, plot is saved after pyplot.show()
            save_filename (str): Filename to save plot. Extension determines filetype.

        """

        self.logs_dirs = logs_dirs
        self.profile_names = profile_names
        self.history_names = history_names
        self.clean_data = clean_data
        self.extra_history_cols = extra_history_cols
        self.identifier = identifier
        self.extractor = extractor
        self.log10_on_data = log10_on_data
        self.contour_colormap = contour_colormap
        self.levels = levels
        self.log_levels = log_levels
        self.num_levels = num_levels
        self.xaxis = xaxis
        self.xaxis_divide = xaxis_divide
        self.time_units = time_units
        self.yaxis = yaxis
        self.yaxis_normalize = yaxis_normalize
        self.show_conv = show_conv
        self.show_therm = show_therm
        self.show_semi = show_semi
        self.show_over = show_over
        self.show_rot = show_rot
        self.core_masses = core_masses
        self.yresolution = yresolution
        self.mass_tolerance = mass_tolerance
        self.radius_tolerance = radius_tolerance
        self.decorate_plot = decorate_plot
        self.show_plot = show_plot
        self.save_file = save_file
        self.save_filename = save_filename
                

#kipp_plot: Plots a Kippenhahn diagram into the matplotlib axis given. No decoration
#           done (i.e. axis labeling or colorbars). Returns
def kipp_plot(kipp_args, axis=None):
    if axis == None:
        fig = plt.figure()
        axis = fig.gca()

    #Fill profile and history names if unspecified
    profile_names = kipp_args.profile_names
    if len(profile_names) == 0:
        profile_names = []
        for log_dir in kipp_args.logs_dirs:
            profile_names.extend(\
                    [log_dir+"/profile"+str(int(i))+".data" for i in np.loadtxt(log_dir+"/profiles.index", skiprows = 1, usecols = (2,))])
    history_names = kipp_args.history_names
    if len(history_names) == 0:
        history_names = []
        for log_dir in kipp_args.logs_dirs:
            history_names.append(log_dir+"/history.data")

    xaxis_divide = kipp_args.xaxis_divide
    if kipp_args.xaxis == "star_age":
        if kipp_args.time_units == "1000 yr":
            xaxis_divide = 1000
        elif kipp_args.time_units == "Myr":
            xaxis_divide = 1e6
        elif kipp_args.time_units == "Gyr":
            xaxis_divide = 1e9

    xyz_data = get_xyz_data(profile_names, xaxis_divide, kipp_args)

    #Get levels if undefined
    levels = kipp_args.levels
    if len(levels) == 0:
        if kipp_args.log_levels:
            levels = get_levels_log(np.nanmin(xyz_data.Z[:,:]), np.nanmax(xyz_data.Z[:,:]), \
                    kipp_args.num_levels)
        else:
            levels = get_levels_linear(np.nanmin(xyz_data.Z[:,:]), np.nanmax(xyz_data.Z[:,:]), \
                    kipp_args.num_levels)
    #make plot
    contour_plot = axis.contourf(xyz_data.X, xyz_data.Y, xyz_data.Z[:,:], \
                cmap=kipp_args.contour_colormap, levels=levels, antialiased = False)

    #zones, mix_types, x_coords, y_coords, histories = get_mixing_zones(\
    #        history_names, xaxis_divide, xyz_data.xlims, kipp_args)
    mixing_zones = get_mixing_zones(history_names, xaxis_divide, xyz_data.xlims, kipp_args)


    for i,zone in enumerate(mixing_zones.zones):
        color = ""
        #Convective mixing
        if mixing_zones.mix_types[i] == 1 and kipp_args.show_conv:
            color = "Chartreuse"
            hatch = "//"
            line  = 1
        #Overshooting 
        elif mixing_zones.mix_types[i] == 3 and kipp_args.show_over:
            color = "purple"
            hatch = "x"
            line  = 1
        #Semiconvective mixing
        elif mixing_zones.mix_types[i] == 4 and kipp_args.show_semi:
            color = "red"
            hatch = "\\\\"
            line  = 1
        #Thermohaline mixing
        elif mixing_zones.mix_types[i] == 5 and kipp_args.show_therm:
            color = "Gold" #Salmon
            hatch = "||"
            line  = 1
        #Rotational mixing
        elif mixing_zones.mix_types[i] == 6 and kipp_args.show_rot:
            color = "brown"
            hatch = "*"
            line  = 1
        #Anonymous mixing
        else: 
            color = "white"
            hatch = " "
            line = 0
        axis.add_patch(PathPatch(zone, fill=False, hatch = hatch, edgecolor=color, linewidth=line))

    #limit x_coords to data of contours and add line at stellar surface
    for i, x_coord in enumerate(mixing_zones.x_coords):
        if x_coord > xyz_data.xlims[1]:
            break
    #I also fill with white above the plot to cover any remainer of the contour plot
    axis.fill_between(mixing_zones.x_coords, max(mixing_zones.y_coords), mixing_zones.y_coords, color = "w")
    axis.plot(mixing_zones.x_coords[:i], mixing_zones.y_coords[:i], "k-")
    #add core masses
    if kipp_args.yaxis == "mass":
        for core_mass in kipp_args.core_masses:
            if core_mass == "He":
                field_name = "he_core_mass"
                color = "b:"
            elif core_mass == "C":
                field_name = "c_core_mass"
                color = "r:"
            elif core_mass == "O":
                field_name = "o_core_mass"
                color = "g:"
            for history in mixing_zones.histories:
                axis.plot(history.get(kipp_args.xaxis) / xaxis_divide, history.get(field_name), color)

    if kipp_args.decorate_plot:
        #add colorbar
        bar = plt.colorbar(contour_plot,pad=0.05)
        bar.set_label(kipp_args.identifier)

        #add axis labels
        if kipp_args.xaxis == "star_age":
            axis.set_xlabel(r'$t$ ('+kipp_args.time_units+')')
        else:
            axis.set_xlabel(r'$Model Number$')
        if kipp_args.yaxis == "radius":
            if kipp_args.yaxis_normalize:
                axis.set_ylabel(r'$r/R$')
            else:   
                axis.set_ylabel(r'$r/R_\odot$')
        else:
            if kipp_args.yaxis_normalize:
                axis.set_ylabel(r'$m/M$')
            else:   
                axis.set_ylabel(r'$m/M_\odot$')

    if kipp_args.show_plot:
        plt.show()
    if kipp_args.save_file:
        plt.savefig(kipp_args.save_filename)

    Kipp_Plot = namedtuple('Kipp_Plot', 'contour_plot histories xlims')

    return Kipp_Plot(contour_plot, mixing_zones.histories, xyz_data.xlims)

