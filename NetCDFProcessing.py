# -*- coding: utf-8 -*-
import os
import re
import urllib2
import math
import errno
from netCDF4 import Dataset, num2date, date2num
from datetime import datetime, timedelta
import numpy as np
from glob import glob
from mpl_toolkits.basemap import Basemap, interp

from variablemetadata import metadata, remotebasepath, bounding_box
from hydrounits.timeseriesconversion import ConvertSeriesTo as CST

LIMIT = 2016
START = 1948

class NCEPDataDownload(object):
    
    def __init__(self):
        super(NCEPDataDownload, self).__init__()
        self.BasePath = remotebasepath
        self.resolution = 1
        
    def setupDataset(self, filepath, variable, resolution = 1):
        if os.path.isfile(filepath):
            self.dataset = Dataset(filepath, "r")
            self.lats = self.dataset.variables['lat']  # extract/copy the data
            self.lons = self.dataset.variables['lon']
            self.times = self.dataset.variables['time']
            self.data = self.dataset.variables[variable]
            print self.lats[:]
        else :
            print ("File does not exist.") 
    
    def createDirectories(self, subdirs, filename = None):
        """
        This function uses the connection result of the remoteConnect method and open an sftp connection to create remote directories
        """
        tempbasepath = self.BasePath
        for each_subdir in subdirs.split("/"):
            tempbasepath = tempbasepath + each_subdir + "/"
            if not os.path.isdir(tempbasepath):
                os.mkdir(tempbasepath)
                
        self.datadirectory = tempbasepath
        return self.datadirectory
    
    def generateDirectoriesAndUrls(self):
        """
        This functions generates the needed data urls and both generate and create the directory where the raw data is to be saved
        """ 
        #/AirTemperature/NCEP/60min/Daily" 
        for key in metadata.keys():
            for subkey in metadata[key].keys():
                #Use the self.remotesubdirectory value to rename processed file???               
#                 if subkey == "Monthly":
#                     filesurls = [metadata[key][subkey]["path"] + metadata[key][subkey]["file"]]

                if subkey == "Daily":
                    self.raw_subdirectory = key + "/NCEP/RAW/" + subkey
                    self.raw_directory = self.createDirectories(self.raw_subdirectory)
                    
                    self.processed_subdirectory = key + "/NCEP/{}min/".format(60*self.resolution) + subkey  
                    self.processed_directory = self.createDirectories(self.processed_subdirectory)
                    
                    filesurls = [metadata[key][subkey]["path"] + metadata[key][subkey]["file"]+str(i)+".nc" for i in range(START, LIMIT, 1)]

                    for i, url in enumerate(filesurls):
                        self.year = re.findall("[(\d+)]{4}", url)
                        remote_con = urllib2.urlopen(url)
                        #try exception???
                        
                        self.raw_file_path = self.raw_directory + url.split("/")[-1]
                        
                        #while not ...:
                        remotefile = open(self.raw_file_path,"wb")
                        remotefile.write(remote_con.read())
                        remotefile.close()
                         
                        self.variable = metadata[key][subkey]["netcdfvariable"]
                        self.processNCEPData(self.raw_file_path, self.variable, key, subkey)
                        
                        #Create the new folder for storing the processed data
                         
                        #self.createDirectories(self.processed_subdirectory)
                        #print self.datadirectory
                        self.processed_file_path = self.processed_directory + "_".join(self.processed_subdirectory.split("/")) + "_" + str(self.year)+ ".nc" 
                        self.createNewNetCDF(self.processed_file_path, self.raw_file_path, self.variable)
                    
    
    def processNCEPData(self, filepath, variable, key, subkey):
        #Set up the dataset  
        if os.path.isfile(filepath):
            self.setupDataset(filepath, variable, 1)                   
        #Change the time axis   
            if len(self.year) != 0:
                self.year = int(self.year[0])
                #Change time axis for daily datasets
                self.changeTimeAxis(self.times, self.year)
            self.unit_abbrev = None                     
            #Convert Units if necessarry
            if metadata[key][subkey]["new_unitcv"] != None :
                self.unit_abbrev = metadata[key][subkey]["abbreviation"]
                self.convertnetCDFUnit(self.data[:], metadata[key][subkey]["unitcv"], metadata[key][subkey]["new_unitcv"])
                
                #Regrid all the data
                self.regrid(self.new_data_values, bounding_box, self.resolution)
            else:
                self.regrid(self.data[:], bounding_box, self.resolution)
        else :
            raise ("File {} does not exist").format(os.path.basename(filepath))
        
    def regrid(self, data, bounding_box, resolution):
        #Setup Mercator basemap
        m = Basemap(projection='merc', llcrnrlat=-90, urcrnrlat=90,\
            llcrnrlon=0, urcrnrlon=360, resolution='c', lon_0=0)
        new_data = []  
        lats = np.flipud(self.lats[:])
        #Create the regularly spaced grid 
        new_lat = np.arange(bounding_box[0][0],bounding_box[0][1], self.resolution)
        new_lon =  np.arange(bounding_box[1][0],bounding_box[1][1], self.resolution)
        
        lon2d, lat2d = np.meshgrid(new_lon, new_lat)
        #Shift the data to draw them correctly according to the map origin (-180 to 180) and also perform interpolation 
        for i in range(len(data)):
            lons, shift_data = m.shiftdata(self.lons[:], datain = data[i,:,:], lon_0=0)
            #The data has to be flipped because interp wants both the lats and lons to be in increasing order
            interp_data = interp(np.flipud(shift_data), lons, lats, lon2d, lat2d, order = 1)
            #we reflip the data again to be in the previous order it was before interpolation
            new_data.append(np.flipud(interp_data)) 
        #As we reflip the data, the latitude most be flipped also   
        self.reg_lon, self.reg_lat = new_lon, np.flipud(new_lat)
        self.new_data_values = np.array(new_data)                  
        
    def changeTimeAxis(self, times, year):
        self.date_time = num2date(times[:],units = times.units, calendar = "gregorian")
        #Need this for the new file
        self.new_times = [delt.days for delt in (self.date_time-datetime(year,1,1))]
        #print self.new_times
        
    def convertnetCDFUnit(self, dataset, origin_unit, destination_unit):
        converted_dataset = []
        data_shape = dataset.shape
        
        if "month" in origin_unit or "month" in destination_unit:
            date_time = num2date(self.times[:],units = self.times.units, calendar = "gregorian")
            for i in range(data_shape[0]):
                values = dataset[i].flatten()
                year = (date_time[i].timetuple()).tm_year
                month = (date_time[i].timetuple()).tm_mon
                new_values = CST([values], [origin_unit], [destination_unit], year = year, month = month)
                converted_dataset.append((np.array(new_values.convert_unit()).reshape(data_shape[1], data_shape[2])))
            valid_range = CST([self.data.valid_range], [origin_unit], [destination_unit], year = year, month = month)
            
        else: 
            values = dataset.flatten()
            #print values, "values"
            new_values = CST([values], [origin_unit], [destination_unit])
            self.new_data_values = np.array(new_values.convert_unit()).reshape(data_shape[0], data_shape[1], data_shape[2])
            valid_range = CST([self.data.valid_range], [origin_unit], [destination_unit] )
        self.valid_range = np.array(valid_range.convert_unit()[0]).astype(np.float64, copy=False) 
        #need the conversion factor to compute the actual range for the converted dataset    
        #self.new_data_values = converted_dataset
        #print self.new_data_values.shape, "convertunit","\n"
    #To be called when all files are downloaded                    
    def mergeWindSpeeds(self, file_path_one, file_path_two, variable_one, variable_two):
        #Check if datafile exists
        if os.path.isfile(file_path_one) == True and os.path.isfile(file_path_two) == True:
            file_buffer_one = Dataset(file_path_one, "r")
            file_buffer_two = Dataset(file_path_two, "r")
            dataset_one = file_buffer_one.variables[variable_one][:]
            dataset_two = file_buffer_two.variables[variable_two][:]
            merged_data = zip(dataset_one.flatten(),dataset_two.flatten())
            merged_dataset = np.array([math.sqrt(x**2+y**2) for (x,y) in merged_data])
            merged_dataset = merged_dataset.reshape(dataset_one.shape)
            self.new_data_values = np.array(merged_dataset)

        else :
            raise ("at least one of the given file names does not exist") 

    def processWindSpeed(self, dir_1, dir_2, var_1, var_2): #To modify and the mergeAndProcessWindSpeed method
        #merging datasets
        file_list_1 = dir_1 + "/*.nc"
        file_list_2 = dir_2 + "/*.nc"
        #
        for i, filepath in enumerate(glob(file_list_1)):
            self.year = int(re.findall("[(\d+)]{4}", filepath)[0])
            self.processed_directory = "TotalWindSpeed/NCEP/{0}min/Daily/TotalWindSpeed_NCEP_{1}min_Daily_{2}.nc".format(60*self.resolution,60*self.resolution,self.year)
            self.createDirectories(os.path.dirname(self.processed_directory))
            self.mergeWindSpeeds(filepath, glob(file_list_2)[i], var_1, var_2)
            self.regrid(self.new_data_values, bounding_box, self.resolution)
            self.setupDataset(filepath, var_1, self.resolution)
            self.changeTimeAxis(self.times, self.year)
            self.unit_abbrev = "m/s"            
            self.createNewNetCDF(remotebasepath + self.processed_directory,filepath ,var_1, "twnd")
            
    def createNewNetCDF(self, new_dataset_path, original_dataset_path, original_variable, new_variable = None):
        #Check if the path exists
        newdataset = Dataset(new_dataset_path, "w", format = self.dataset.file_format)
        #The following 4 lines allow one to send either a file name to this method or a netCDF4.Dataset Object
        if "netCDF4.Dataset" in str(type(original_dataset_path)):
            original_dataset = original_dataset_path        
        elif os.path.isfile(original_dataset_path):
            original_dataset = Dataset(original_dataset_path, "r")

        #Set attributes for the dataset to be created   
        newdataset.Conventions = original_dataset.Conventions
        newdataset.title = original_dataset.title
        newdataset.description = original_dataset.description
        newdataset.platform = original_dataset.platform
        newdataset.references = original_dataset.references
        newdataset.history = original_dataset.history
        
        data = original_dataset.variables[original_variable]

        #Create dimensions
        lon = newdataset.createDimension(unicode("lon"), len(self.reg_lon[:]))
        lat = newdataset.createDimension(unicode("lat"), len(self.reg_lat[:]))
        time = newdataset.createDimension(unicode("time"), None)

        #create latitude dimension and variable
        lat_var = newdataset.createVariable('lat', 'f8',\
                                       ('lat',))
        for ncattr in original_dataset.variables['lat'].ncattrs():
            lat_var.setncattr(ncattr, original_dataset.variables['lat'].getncattr(ncattr))
        #lat_var.actual_range = [min(lat_var[:]), max(self.data[:])]
        newdataset.variables['lat'][:] = np.array(self.reg_lat[:]).astype(np.float64, copy=False)
        
        #create longitude variable
        lon_var = newdataset.createVariable('lon', 'f8',\
                                       ('lon',))
        newdataset.variables['lon'][:] = self.reg_lon[:]
        for ncattr in original_dataset.variables['lon'].ncattrs():
            lon_var.setncattr(ncattr, original_dataset.variables['lon'].getncattr(ncattr))        
        #lon_var.actual_range = [min(lat_var[:]), max(self.data[:])]
                        
        #create time variables
        time_var = newdataset.createVariable('time', 'i',\
                                       ('time',))
        for ncattr in original_dataset.variables['time'].ncattrs():
            time_var.setncattr(ncattr, original_dataset.variables['time'].getncattr(ncattr))
        
        if len([self.year]) != 0:#for the daily datasets
            time_var.units = "days since %s-01-01 00:00:00.0" %self.year
            time_var.actual_range = [min(self.new_times), max(self.new_times)]
            newdataset.variables['time'][:] = np.array(self.new_times[:]).astype(np.int, copy=False)
            
        else:
            newdataset.variables['time'][:] = np.array(self.times[:]).astype(np.int, copy=False) #for the monthly datasets
        
        #create data/value variable
        if new_variable != None:
            data_var = newdataset.createVariable(new_variable, 'f8', ('time','lat','lon',),\
                                                  zlib=False, least_significant_digit= self.data.least_significant_digit)
            newdataset.variables[new_variable][:] = self.new_data_values.astype(np.float64, copy=False)
            
        else : 
            data_var = newdataset.createVariable(original_variable, 'f8', ('time','lat','lon',), \
                                                 zlib=False, least_significant_digit= self.data.least_significant_digit)
            newdataset.variables[original_variable][:] = self.new_data_values.astype(np.float64, copy=False)
             
        for ncattr in original_dataset.variables[original_variable].ncattrs():
            data_var.setncattr(ncattr, original_dataset.variables[original_variable].getncattr(ncattr))
        #Reset valid range for the new dataset
        data_var.valid_range = [(self.valid_range[0]), int(self.valid_range[1])]
        print data_var.valid_range
        #Set actual range attribute 
        if new_variable != "twnd":
            data_var.actual_range = [min((newdataset.variables[original_variable][:]).flatten()), \
                                 max((newdataset.variables[original_variable][:]).flatten())]         
        elif new_variable == "twnd":
            data_var.long_name = "mean Daily total wind at 10 m"
            data_var.actual_range = [min((newdataset.variables[new_variable][:]).flatten()), \
                                 max((newdataset.variables[new_variable][:]).flatten())]  
        #Set the unit abbreviation attribute      
        if self.unit_abbrev != None:
            data_var.units = self.unit_abbrev
        
        newdataset.close()
   
if __name__ == "__main__":
    ncep = NCEPDataDownload()
    filepath = "air.2m.gauss.2013.nc"
    ncep.setupDataset(filepath, "air", 1)
    #print ncep.data[:]
    ncep.convertnetCDFUnit(ncep.data[:], "degree kelvin", "degree celsius")
    ncep.regrid(ncep.new_data_values, bounding_box, 1)
    ncep.changeTimeAxis(ncep.times, 2013)
    ncep.year = 2013
    ncep.unit_abbrev = "DegC"
    ncep.createNewNetCDF("AirTemperature_NCEP_60min_Daily_2013_last.nc", filepath, "air")
    print "ok done"
    #ncep.generateDirectoriesAndUrls()
    #ncep.processWindSpeed()
    #merging datasets
#     file_list_1 = "C:/data/NETCDFArchive/Global/NorthSouthWindSpeed/NCEP/RAW/Daily/*.nc"
#     file_list_2 = "C:/data/NETCDFArchive/Global/WestEastWindSpeed/NCEP/RAW/Daily/*.nc"
#     for i, filepath in enumerate(glob(file_list_1)):
#         ncep.year = int(re.findall("[(\d+)]{4}", filepath)[0])
#         ncep.processed_directory = "TotalWindSpeed/NCEP/{0}min/Daily/TotalWindSpeed_NCEP_{1}min_Daily_{2}.nc".format(60*ncep.resolution,60*ncep.resolution,ncep.year)
#         ncep.createDirectories(os.path.dirname(ncep.processed_directory))
#         #print remotebasepath+"TotalWindSpeed/NCEP/{0}min/Daily/TotalWindSpeed_NCEP_{1}min_Daily_{2}.nc".format(60*ncep.resolution,60*ncep.resolution,ncep.year)
#         ncep.mergeAndprocessWindSpeed(filepath, glob(file_list_2)[i], "vwnd", "uwnd")
#         ncep.regrid(ncep.new_data_values, bounding_box, ncep.resolution)
#         ncep.createNewNetCDF(remotebasepath+ncep.processed_directory,filepath ,"vwnd", "twnd")
    
#     file_path_one = "uwnd.10m.gauss.1948.nc"
#     file_path_two = "vwnd.10m.gauss.1948.nc"
#     ncep.mergeValues(file_path_one, file_path_two, "uwnd", "vwnd")
#"/data/RGISarchive/Global/AirTemperature/NCEP/60min/Daily"    
