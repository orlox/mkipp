from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mesa_data import *
import re

# Someday, if I remember, I might explain how this works
# Class used to store mixing zones
class Zone_Block:
    def __init__(self, min_x, max_x, min_y, max_y):
        #vertices for this block, in the order:
        #-lower left
        #-upper left
        #-upper right
        #-lower right
        self.vertices = [
                Zone_Vertex((min_x,min_y)),
                Zone_Vertex((min_x,max_y)),
                Zone_Vertex((max_x,max_y)),
                Zone_Vertex((max_x,min_y))]
        for i, vertex in enumerate(self.vertices):
            vertex.prev_vertex = self.vertices[i-1]
            vertex.next_vertex = self.vertices[(i+1)%4]

class Zone_Vertex:
    def __init__(self, coords):
        self.coords = coords
        self.prev_vertex = None
        self.next_vertex = None
        self.checked = False

class Zone:
    def __init__(self, block, mix_type):
        self.old_blocks = []
        self.last_blocks = []
        self.new_blocks = [block]
        self.mix_type = mix_type
        self.open = True

    def extend(self, upper_vertex, lower_vertex, mix_type):
        if mix_type != self.mix_type:
            return False
        found = False
        for block in self.last_blocks:
            old_upper = block.vertices[2]
            old_lower = old_upper.next_vertex
            old_max_y = old_upper.coords[1]
            old_min_y = old_lower.coords[1]
            new_max_y = upper_vertex.coords[1]
            new_min_y = lower_vertex.coords[1]
            if (new_max_y <= old_max_y and new_max_y >= old_min_y) or \
               (new_min_y <= old_max_y and new_min_y >= old_min_y) or \
               (new_max_y > old_max_y and new_min_y < old_min_y):
               #associate vertices
               old_upper.next_vertex = upper_vertex
               upper_vertex.prev_vertex = old_upper
               lower_vertex.next_vertex = old_lower
               old_lower.prev_vertex = lower_vertex
               #if upper vertex is contained in an older block, no need to
               #continue
               if (new_max_y <= old_max_y and new_max_y >= old_min_y) and \
                  (new_min_y <= old_max_y and new_min_y >= old_min_y):
                  return True
               #otherwise, switch lower vertex
               lower_vertex = old_upper
               found = True
        return found

    def switch_new_blocks(self):
        self.old_blocks.extend(self.last_blocks)
        self.last_blocks = self.new_blocks
        self.new_blocks = []

    def get_path(self):
        self.old_blocks.extend(self.last_blocks)
        vertices = []
        coords = []
        codes = []
        for block in self.old_blocks:
            vertices.extend(block.vertices)

        for starting_vertex in vertices:
            if starting_vertex.checked:
                continue
            starting_vertex.checked = True
            coords.append(starting_vertex.coords)
            codes.append(Path.MOVETO)
            current_vertex = starting_vertex
            while current_vertex.next_vertex != starting_vertex and not current_vertex.next_vertex.checked:
                current_vertex = current_vertex.next_vertex
                current_vertex.checked = True
                coords.append(current_vertex.coords)
                codes.append(Path.LINETO)
            coords.append((0,0))
            codes.append(Path.CLOSEPOLY)

        return Path(coords, codes)

    def merge_zone(self, zone2):
        self.old_blocks.extend(zone2.old_blocks)
        self.last_blocks.extend(zone2.last_blocks)
        self.new_blocks.extend(zone2.new_blocks)

#returns a list of matplotlib path objects with the mixing regions
def get_mixing_zones(history_names, xaxis_divide, min_x_coord, max_x_coord, kipp_args):
   # Get mixing zones and mass borders from history.data files
   print "Reading history data"
   mix_data = []
   histories = []
   for history_name in history_names:
       history = Mesa_Data(history_name, read_data = False)
       columns = []
       for key in history.columns.keys():
           search_regex = "log_R|star_mass|model_number|.*core_mass|"
           if kipp_args.yaxis == "radius":
               search_regex = search_regex + "mix_relr_top.*"
           else:
               search_regex = search_regex + "mix_qtop.*"
           if re.search(search_regex,key):
               columns.append(key)
       history.read_data(columns)
       histories.append(Mesa_Data(history_name))
   x_coords = []
   for history in histories:
       x_coords.extend(history.get(kipp_args.xaxis) / xaxis_divide)
   y_coords = []
   if kipp_args.yaxis_normalize:
       y_coords = [1.0]*len(x_coords)
   elif kipp_args.yaxis == "radius":
       for history in histories:
           y_coords.extend(np.power(10,history.get('log_R')))
   else:
       for history in histories:
           y_coords.extend(history.get('star_mass'))

   print "Constructing mixing regions"
   mesa_mix_zones = 0
   while True:
       try:
           mesa_mix_zones = mesa_mix_zones + 1
           mix_type = []
           mix_top = []
           for history in histories:
               mix_type.extend(history.get('mix_type_'+str(mesa_mix_zones)))
               if kipp_args.yaxis == "radius":
                    mix_top.extend(history.get('mix_relr_top_'+str(mesa_mix_zones)))
               else:
                    mix_top.extend(history.get('mix_qtop_'+str(mesa_mix_zones)))
           mix_data.append([mix_type, mix_top])
       except Exception, e:
           #reached all mix zones included
           mesa_mix_zones = mesa_mix_zones - 1
           break

   if kipp_args.yaxis == "radius":
       tolerance = kipp_args.radius_tolerance
   else:
       tolerance = kipp_args.mass_tolerance

   zones = []
   mix_types = []
   open_zones = []
   new_zones = []
   for i in range(1,len(x_coords)):
       current_x = x_coords[i]
       if current_x > max_x_coord:
           break
       if current_x < min_x_coord:
           continue
       previous_x = x_coords[i-1]
       for j in range(0,mesa_mix_zones):
           mix_type = mix_data[j][0][i]
           if mix_type == 0 or mix_type == -1:
               continue 
           max_y_coord = mix_data[j][1][i]*y_coords[i]
           min_y_coord = 0
           if j > 0:
               min_y_coord = mix_data[j-1][1][i]*y_coords[i]
           #ignore too small regions
           if max_y_coord - min_y_coord < tolerance*y_coords[i]:
               continue
           zone_block = Zone_Block(previous_x, current_x, min_y_coord, max_y_coord)
           exists = False
           zones_to_merge = []
           for z in open_zones:
               if z.extend(zone_block.vertices[1], zone_block.vertices[1].prev_vertex, mix_type):
                   exists = True
                   z.new_blocks.append(zone_block)
                   zones_to_merge.append(z)
           #merge zones as needed
           for k in range(1,len(zones_to_merge)):
               zones_to_merge[0].merge_zone(zones_to_merge[k])
               open_zones.remove(zones_to_merge[k])
           #create zone if it has no predecesor
           if not exists:
               z = Zone(zone_block, mix_type)
               new_zones.append(z)
       open_zones.extend(new_zones)
       new_zones = []
       #separate zones which didn't continue here so we don't need to check them all the time
       for z in open_zones:
           z.switch_new_blocks()
       for z in open_zones:
           if len(z.last_blocks) == 0:
               zones.append(z.get_path())
               mix_types.append(z.mix_type)
               open_zones.remove(z)

   for z in open_zones:
       z.switch_new_blocks()
       zones.append(z.get_path())
       mix_types.append(z.mix_type)

   return zones, mix_types, x_coords, y_coords, histories

def get_xyz_data(profile_names, xaxis_divide, kipp_args):
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
    # Initialize interpolation grid
    Z_data_array = np.zeros((kipp_args.yresolution,len(profile_names)))
    # XY coordinates for data
    X_data_array = np.zeros((kipp_args.yresolution,len(profile_names)))
    Y_data_array = np.zeros((kipp_args.yresolution,len(profile_names)))
    for j in range(kipp_args.yresolution):
        Y_data_array[j,:] = max_y * j / (kipp_args.yresolution-1)
    for i,profile_name in enumerate(profile_names):
        try:
            prof = Mesa_Data(profile_name, read_data = False)
            prof.read_data(columns)
        except Exception as e:
            print "Couldn't read profile " + profile_name
        x_coord = prof.header[kipp_args.xaxis] / xaxis_divide
        if x_coord < max_x_coord:
            print "Profiles are not ordered in X coordinate!!!"
        max_x_coord = max(max_x_coord, x_coord)
        min_x_coord = min(min_x_coord, x_coord)
        if kipp_args.yaxis == "mass":
            prof_y = prof.header['star_mass']
        elif kipp_args.yaxis == "radius":
            prof_y = prof.header['photosphere_r']

        #fill up X positions
        X_data_array[:,i] = x_coord

        #read and interpolate data
        if kipp_args.yaxis == "mass":
            y_data = prof.get('mass')
        elif kipp_args.yaxis == "radius":
            y_data = prof.get('radius')
        #reverse y_data and z_data for np.interp
        y_data = y_data[::-1]
        z_data = kipp_args.extractor(kipp_args.identifier, kipp_args.log10_on_data, prof)
        z_data = z_data[::-1]
        interp_z_data = np.interp(y_interp, y_data, z_data)
        #set nans outside of range so there is no plotting
        #for j in range(kipp_args.yresolution):
        #    if (Y_data_array[j,i] > prof_y):
        #        Z_data_array[j,i] = np.nan
        #    else:
        Z_data_array[:,i] = interp_z_data

    return min_x_coord, max_x_coord, X_data_array, Y_data_array, Z_data_array
