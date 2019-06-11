import numpy as np
import numpy.ma as ma
import h5py
from collections import OrderedDict

#class to extract a mesa data file.
class mesa_data:
    def __init__(self, history_file, only_read_header = False, read_data = True, read_data_cols = [],
            clean_data = True, sample_every_n = 1, is_hdf5 = False):
        self.filename = history_file
        self.is_hdf5 = is_hdf5
        #header is a dictionary with the general info from the second and third line of file
        self.header = {}
        self.header_num = {}
        #columns is a dictionary which gives the column number (minus 1) corresponding to the key
        self.columns = {}
        columns = []
        if not is_hdf5:
            file = open(self.filename, "r")
            #first line is not used
            file.readline()
            #following two lines have header data
            header_names = file.readline().split()
            header_vals = file.readline().split()
            i = 0
            for i, header_name in enumerate(header_names):
                #need to properly account for these new columns
                if header_name in ["compiler", "build", "MESA_SDK_version", "date"]:
                    continue
                self.header[header_name] = float(header_vals[i])
                self.header_num[header_name] = i
                i+=1
            if only_read_header:
                file.close()
                return
            #next line is empty
            file.readline()
            #following two lines have column data
            nums = file.readline().split()
            names = file.readline().split()
            for i, name in enumerate(names):
                self.columns[name] = int(nums[i])-1
                columns.append(name)
            file.close()
        else:
            file = h5py.File(self.filename, "r")
            header_names = file['header_names'][:]
            header_vals = file['header_vals'][:]
            for i in range(len(header_names)):
                key = header_names[i].decode('utf-8')
                self.header[key] = header_vals[i]
                self.header_num[key] = i
            columns = file['data_names'][:].tolist()
            for i, col in enumerate(columns):
                self.columns[col.decode('utf-8')] = i
                columns[i] = col.decode('utf-8')
            file.close()

        if not read_data:
            return

        if len(read_data_cols) == 0:
            read_data_cols = columns
        self.read_data(read_data_cols, clean_data = clean_data)


    def read_data(self, column_names, clean_data = True, sample_every_n = 1):
        #always include model_number if its part of the data
        if "model_number" not in column_names and "model_number" in self.columns:
            column_names.append("model_number")

        #be sure there are no repeated column names
        #(could use set but that breaks the order of the columns, which is needed if I want to save the file)
        column_names = list(OrderedDict.fromkeys(column_names))

        self.read_columns = column_names

        if not self.is_hdf5:
            #read data
            data = np.loadtxt(self.filename, skiprows = 6, \
                usecols = tuple([self.columns[k] for k in column_names]), unpack = True)
        else:
            file = h5py.File(self.filename, "r")
            data = file['data_vals'][:,sorted([self.columns[k] for k in column_names])]
            data = data.transpose()
            file.close()

        self.data = {}
        #Be careful in case only one column is required
        if len(data.shape) > 1:
            for i, column in enumerate(column_names):
                self.data[column] = data[i]
        else:
            self.data[column_names[0]] = data

        #clean redos
        if clean_data and "model_number" in self.columns and len(self.data["model_number"]) > 1:
            #create a mask
            model_number = self.data["model_number"]
            mask = np.zeros(len(model_number))
            max_model_number = model_number[-1]
            #last entry is valid, start from there and remove repeats
            for i in range(len(model_number)-2,-1,-1):
                if model_number[i] >= max_model_number:
                    #exclude this point
                    mask[i] = 1
                else:
                    max_model_number = model_number[i]

            if sum(mask) > 0:
                for column in column_names:
                    self.data[column] = ma.masked_array(self.data[column], mask=mask).compressed()

        #subsample points
        if sample_every_n > 1 and "model_number" in self.columns and len(self.data["model_number"]) > 2:
            #keep first and last entry
            #create a mask
            model_number = self.data["model_number"]
            mask = np.zeros(len(model_number))
            for i in range(1,len(model_number)-1):
                if (i+1)%sample_every_n != 0:
                    #exclude this point
                    mask[i] = 1

            if sum(mask) > 0:
                for column in column_names:
                    self.data[column] = ma.masked_array(self.data[column], mask=mask).compressed()

        #count number of points using first entry in dict
        self.num_points = len(self.data[list(self.read_columns)[0]])

    def get(self,key):
        return self.data[key]

    def save_as_hdf5(self, filename, header_str_dtype="S28", data_str_dtype="S40", compression_opts=4):
        f = h5py.File(filename, "w")
        dset_header_names = f.create_dataset("header_names", (len(self.header),), dtype=header_str_dtype)
        dset_header_vals = f.create_dataset("header_vals", (len(self.header),), dtype="d")
        for key in self.header:
            dset_header_names[self.header_num[key]] = np.string_(key)
            dset_header_vals[self.header_num[key]] = self.header[key]
        dset_column_names = f.create_dataset("data_names", (len(self.read_columns),), dtype=data_str_dtype)
        dset_column_vals = f.create_dataset("data_vals", (self.num_points,len(self.read_columns)), dtype="d",
                compression='gzip',compression_opts=compression_opts)
        for k, key in enumerate(self.read_columns):
            dset_column_names[k] = np.string_(key)
            dset_column_vals[:,k] = self.data[key]
        f.close()

    #creates a mesa look-alike output file
    #prints all integers as doubles
    #not the most efficient code but I don't care
    def save_as_ascii(self, filename, header_str_format="{0:>28}", header_double_format="{0:>28.16e}",
            data_str_format="{0:>40}", data_double_format="{0:>40.16e}"):
        f = open(filename, "w")
        for i in range(len(list(self.header))):
            f.write(header_str_format.format(i+1))
        f.write("\n")
        #create an ordered list of keys
        header_keys = []
        for i in range(len(list(self.header))):
            for key in self.header:
                if self.header_num[key] == i:
                    header_keys.append(key)
                    break
        for i, key in enumerate(header_keys):
            f.write(header_str_format.format(key))
        f.write("\n")
        for i, key in enumerate(header_keys):
            f.write(header_double_format.format(self.header[key]))
        f.write("\n")
        f.write("\n")

        for i in range(len(list(self.read_columns))):
            f.write(data_str_format.format(i+1))
        f.write("\n")

        for i, key in enumerate(self.read_columns):
            f.write(data_str_format.format(key))
        for k in range(self.num_points):
            f.write("\n")
            for i, key in enumerate(self.read_columns):
                f.write(data_double_format.format(self.data[key][k]))

        f.close()



#reads the profiles.index files in the folders specified by the logs_dirs array and returns
#an array containing paths to the individual profile files, after cleaning up redos and backups
def get_profile_paths(logs_dirs = ["LOGS"]):
    profile_paths = []
    for log_dir in logs_dirs:
        print(log_dir, logs_dirs)
        model_number, paths = np.loadtxt(log_dir+"/profiles.index", skiprows = 1, usecols = (0,2), unpack = True)
        mask = np.zeros(len(paths))
        max_model_number = model_number[-1]
        #last entry is valid, start from there and remove repeats
        for i in range(len(model_number)-2,-1,-1):
            if model_number[i] >= max_model_number:
                mask[i] = 1
            else:
                max_model_number = model_number[i]

        if sum(mask) > 0:
            paths = ma.masked_array(paths, mask=mask).compressed()
        profile_paths.extend([log_dir+"/profile"+str(int(i))+".data" for i in paths])
    return profile_paths

