import numpy as np

#class to extract a mesa data file.
class Mesa_Data:
    def __init__(self, history_file, only_read_header = False, read_data = True, read_data_cols = []):
        self.filename = history_file
        #header is a dictionary with the general info from the second and third line of file
        self.header = {}
        #columns is a dictionary which gives the column number (minus 1) corresponding to the key
        self.columns = {}
        file = open(self.filename, "r")
        #first line is not used
        file.readline()
        #following two lines have header data
        header_names = file.readline().split()
        header_vals = file.readline().split()
        for i, header_name in enumerate(header_names):
            self.header[header_name] = float(header_vals[i])
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
        file.close()

        if not read_data:
            return

        if len(read_data_cols) == 0:
            read_data_cols = self.columns.keys()
        self.read_data(self.columns.keys())

    def read_data(self, column_names):
        #read data, TODO: ignore inexistant columns
        data = np.loadtxt(self.filename, skiprows = 6, \
            usecols = ([self.columns[k] for k in column_names]), unpack = True)

        self.data = {}
        for i, column in enumerate(column_names):
            self.data[column] = data[i]

    def get(self,key):
        return self.data[key]
