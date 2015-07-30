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
from mixing_zones import *

#matplotlib specifics
import matplotlib.pyplot as plt
from matplotlib.path import Path

# Create contour levels array (TBD: Need to be improved)
def get_levels_linear(min,max,n_levels):
    delta = (max-min)/n_levels
    print min, max, delta, n_levels
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

    # Initialize interpolation grids
    Z_data_array = np.zeros((kipp_args.yresolution,len(profile_names)))

    # XY coordinates for data
    X_data_array = np.zeros((kipp_args.yresolution,len(profile_names)))
    Y_data_array = np.zeros((kipp_args.yresolution,len(profile_names)))

    # Extract data from profiles
    max_x_coord = -1
    min_x_coord = 1e99
    #first read all headers to determine max value on yaxis
    max_y = 0
    print "Reading profile data"
    if kipp_args.yaxis_normalize:
        max_y = star_mass = star_radius = 1.0
    else:
        for i,profile_name in enumerate(profile_names):
            try:
                prof = Mesa_Data(profile_name, only_read_header = True)
                if kipp_args.yaxis == "mass":
                    max_y = max(prof.header['star_mass'],max_y)
                elif kipp_args.yaxis == "radius":
                    max_y = max(prof.header['photosphere_r'],max_y)
            except Exception as e:
                print "Couldn't read profile " + profile_name
    #array to interpolate data in the yaxis
    y_interp = np.array([max_y * j / (kipp_args.yresolution-1) for j in range(kipp_args.yresolution)])
    #now read the data
    #"columns" is a list with the neccesary columns that will be parsed from the profile*.data files
    columns = []
    if kipp_args.yaxis == "mass":
        columns.append('mass')
    elif kipp_args.yaxis == "radius":
        columns.append('radius')
    columns.extend(kipp_args.extractor(\
            kipp_args.identifier, kipp_args.log10_on_data, prof, return_data_columns = True))
    for i,profile_name in enumerate(profile_names):
        try:
            prof = Mesa_Data(profile_name, read_data = False)
        except Exception as e:
            print "Couldn't read profile " + profile_name
        x_coord = prof.header[kipp_args.xaxis] / xaxis_divide
        if x_coord < max_x_coord:
            print "Profiles are not ordered in X coordinate!!!"
        max_x_coord = max(max_x_coord, x_coord)
        min_x_coord = min(min_x_coord, x_coord)

        #fill up positions
        for j in range(kipp_args.yresolution):
            X_data_array[j,i] = x_coord
            Y_data_array[j,i] = max_y * j / (kipp_args.yresolution-1)

        prof.read_data(columns)
        
        #read and interpolate data
        y_data = prof.get(kipp_args.yaxis)
        #reverse y_data and z_data for np.interp
        y_data = y_data[::-1]
        z_data = kipp_args.extractor(kipp_args.identifier, kipp_args.log10_on_data, prof)
        z_data = z_data[::-1]
        interp_z_data = np.interp(y_interp, y_data, z_data)
        #for j in range(kipp_args.numy):
        Z_data_array[:,i] = interp_z_data[:]

    #Get levels if undefined
    levels = kipp_args.levels
    if len(levels) == 0:
        if kipp_args.log_levels:
            levels = get_levels_log(np.min(Z_data_array[:,:]), np.max(Z_data_array[:,:]), \
                    kipp_args.num_levels)
        else:
            levels = get_levels_linear(np.min(Z_data_array[:,:]), np.max(Z_data_array[:,:]), \
                    kipp_args.num_levels)
    #make plot
    contour_plot = axis.contourf(X_data_array, Y_data_array, Z_data_array[:,:], \
                cmap=kipp_args.contour_colormap, levels=levels, antialiased = False)

    zones, mix_types, x_coords, y_coords, histories = get_mixing_zones(\
            kipp_args.logs_dirs, history_names, kipp_args.yaxis_normalize, kipp_args.xaxis, kipp_args.yaxis, \
            xaxis_divide, kipp_args.radius_tolerance, kipp_args.mass_tolerance, min_x_coord, max_x_coord)

    for i,zone in enumerate(zones):
        color = ""
        #Convective mixing
        if mix_types[i] == 1 and kipp_args.show_conv:
            color = "Chartreuse"
            hatch = "//"
            line  = 1
        #Overshooting 
        elif mix_types[i] == 3 and kipp_args.show_over:
            color = "purple"
            hatch = "x"
            line  = 1
        #Semiconvective mixing
        elif mix_types[i] == 4 and kipp_args.show_semi:
            color = "red"
            hatch = "\\\\"
            line  = 1
        #Thermohaline mixing
        elif mix_types[i] == 5 and kipp_args.show_therm:
            color = "Gold" #Salmon
            hatch = "||"
            line  = 1
        #Rotational mixing
        elif mix_types[i] == 6 and kipp_args.show_rot:
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
    for i, x_coord in enumerate(x_coords):
        if x_coord > max_x_coord:
            break
    axis.plot(x_coords[:i], y_coords[:i], "k-")
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
            for history in histories:
                axis.plot(history.get(kipp_args.xaxis) / xaxis_divide, history.get(field_name), color)

    if kipp_args.decorate_plot:
        #add colorbar
        bar = plt.colorbar(contour_plot,pad=0.05)
        bar.set_label('$\epsilon_{nuc}-\epsilon_{\\nu}$,  Log (erg/g/s)')

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

    return contour_plot, histories, [min_x_coord, max_x_coord]

