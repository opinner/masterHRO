import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np
import geopy.distance as geo #https://pypi.org/project/geopy/
import math
from scipy.integrate import simps

#searches for the array item, that is nearest to "search"
def which_arrayindex(search,array):
    deviation = 0
    while deviation < 1:
        
        for i in np.arange(array.size):
            if np.abs(array[i]-search)<=deviation: 
                return i
            
        deviation += 0.0001
        #deviation = deviation*2
    
#returns the Index of the first nonNan item in a array
def firstIndexNonNan(array):
  for index in np.arange(array.size):
    if math.isnan(array[index]) == False:
      return index

#returns the first nonNaN item of array 
def firstNonNan(array):
  for item in array:
    if math.isnan(item) == False:
      return item


def save_mat_txt(topo,filename="warnowdata.txt"):
    mat = np.matrix(topo)
    with open(filename,'wb') as f:
        for line in mat:
            np.savetxt(f, line, fmt='%.2f')


#https://stackoverflow.com/questions/33282368/plotting-a-2d-heatmap-with-matplotlib
#https://matplotlib.org/api/_as_gen/matplotlib.pyplot.imshow.html


def problem1(lon,lat,topo):
    """Plot the topography of the Warnow River!"""
    plt.imshow(topo, cmap='Blues_r', interpolation='nearest') #np.flip(topo, 0)

    plt.xticks(np.arange(0,lon.size,60), np.around(lon[::60],2))
    plt.yticks(np.arange(0,lat.size,50), np.around(lat[::50],2))

    plt.xlabel("longitude /degrees")
    plt.ylabel("latitude /degrees")
    plt.title("topography of the warnow river mouth")

    cbar = plt.colorbar()
    cbar.set_label('river depth / m')


def problem2(lon,lat,topo):
    """What is the maximum and mean depth of the river?"""
    #We cut off the baltic sea and only look at data points in the river.  
    print("maximum depth:",np.nanmin(topo[0:120,:]))
    print("mean depth:",np.nanmean(topo[0:120,:]))


def problem3(lon,lat,topo):
    """Calculate the mean horizontal resolution of the data set!"""
    
    #TODO Kurze Erklärung
    
    difference_degrees = np.ones(lon.size-1)
    difference_meter = np.ones(lon.size-1)

    mean_lat = np.mean(lat)
    mean_lon = np.mean(lon)

    for i in np.arange(difference_degrees.size):
        difference_degrees[i]= np.abs(lon[i]-lon[i+1])
        difference_meter[i] = geo.distance((mean_lat,lon[i]),(mean_lat,lon[i+1])).km*1000

    mean_diff_degrees = np.mean(difference_degrees)
    mean_diff_meter = np.mean(difference_meter)
    print("mean resolution in degrees:",np.around(mean_diff_degrees,7))
    print("mean resolution in meters:",np.around(mean_diff_meter,2),"m")
    
    
def problem4(lon,lat,topo,zonal_cut = 54+7/60):
    """Plot the water depth against distance along a zonal section (54°7'N)!"""
    
    #To find the corresponding line in the data     
    i = which_arrayindex(zonal_cut,lat)
    
    topo_zonal_cut = topo[i,:]
    
    index_west_shore = firstIndexNonNan(topo_zonal_cut)
    index_east_shore = (topo_zonal_cut.size-1) - firstIndexNonNan(topo_zonal_cut[::-1])
    
    lon_in_meter = np.zeros(index_east_shore-index_west_shore+1)
    difference_meter = np.zeros(index_east_shore-index_west_shore)

    for i in np.arange(difference_meter.size):
        difference_meter[i] = geo.distance((zonal_cut,lon[index_west_shore+i]),(zonal_cut,lon[index_west_shore+i+1])).km*1000
        lon_in_meter[i+1] = lon_in_meter[i] +  difference_meter[i]
    
    #Integrating the topography with Simpson's rule. 
    cross_section = simps(topo_zonal_cut[index_west_shore:index_east_shore+1],lon_in_meter)
    print("Cross section:",np.abs(cross_section))
    
    plt.plot(lon_in_meter,topo_zonal_cut[index_west_shore:index_east_shore+1])
    
    plt.xlabel("distance / meter")
    plt.ylabel("river depth / meter")
    
    
def problem5(lon,lat,topo,zonal_cut = 54+7/60):
    """Calculate the width of the River (in m ) 
    and its cross section (in m 2 ) at that latitude!"""
    
    i = which_arrayindex(zonal_cut,lat)
    
    topo_zonal_cut = topo[i,:]
    
    
    index_west_shore = firstIndexNonNan(topo_zonal_cut)
    index_east_shore = (topo_zonal_cut.size-1) - firstIndexNonNan(topo_zonal_cut[::-1])
    
    lon_west_shore = lon[index_west_shore]
    lon_east_shore = lon[index_east_shore]
    
    lon_in_meter = np.zeros(index_east_shore-index_west_shore+1)
    difference_meter = np.zeros(index_east_shore-index_west_shore)
    
    #river width as distance between west and east shore
    river_width = geo.distance((zonal_cut,lon_west_shore),(zonal_cut,lon_east_shore)).km*1000
    
    print("river_width at 54 7':",np.around(river_width,2))
    
    for i in np.arange(difference_meter.size):
        difference_meter[i] = geo.distance((zonal_cut,lon[index_west_shore+i]),(zonal_cut,lon[index_west_shore+i+1])).km*1000
        lon_in_meter[i+1] = lon_in_meter[i] +  difference_meter[i]
    
    cross_section = np.abs(simps(topo_zonal_cut[index_west_shore:index_east_shore+1],lon_in_meter))
    print("Cross section:",cross_section)
    
    
def problem6(lon,lat,topo,zonal_cut = 54+7/60): 
    """Estimate the transport through this cross section"""
    
    mean_v = 0.013 # m/s
    
    #For the sake of completeness: Again calculating the cross section
    i = which_arrayindex(zonal_cut,lat)
    
    topo_zonal_cut = topo[i,:]
    
    index_west_shore = firstIndexNonNan(topo_zonal_cut)
    index_east_shore = (topo_zonal_cut.size-1) - firstIndexNonNan(topo_zonal_cut[::-1])
    
    lon_west_shore = lon[index_west_shore]
    lon_east_shore = lon[index_east_shore]
    
    lon_in_meter = np.zeros(index_east_shore-index_west_shore+1)
    difference_meter = np.zeros(index_east_shore-index_west_shore)
    
    river_width = geo.distance((zonal_cut,lon_west_shore),(zonal_cut,lon_east_shore)).km*1000
    
    print("river_width at 54 7':",np.around(river_width,2),"m")
    
    for i in np.arange(difference_meter.size):
        difference_meter[i] = geo.distance((zonal_cut,lon[index_west_shore+i]),(zonal_cut,lon[index_west_shore+i+1])).km*1000
        lon_in_meter[i+1] = lon_in_meter[i] +  difference_meter[i]
    
    cross_section = np.abs(simps(topo_zonal_cut[index_west_shore:index_east_shore+1],lon_in_meter))
    print("Cross section:",cross_section)
  
  
    transport_per_s = cross_section*mean_v
    print("Transport per second:",np.around(transport_per_s,3),"m^3/s")
    print("Transport per year:",np.around(transport_per_s*3.154*10**7,2),"m^3/a")
    print("Transport per year:",np.around(transport_per_s*3.154*10**7*10**(-9),2),"km^3/a")
  
    
if __name__ == "__main__":

    #open .mat files in python:
    #https://docs.scipy.org/doc/scipy/reference/tutorial/io.html
    warnow_data = sio.loadmat("topo_warnow.mat")
    #print(sio.whosmat("topo_warnow.mat"))

    topo = warnow_data["topo"]
    lat = warnow_data["lat"][0]
    lon = warnow_data["lon"][0]

    #I am not really sure why, but to get the correct display of the topography
    #I have to flip the array
    topo = -1*np.flip(topo, 0)
    lat = np.flip(lat)


    #List of the problems. All of them function indenpendently of each other

    #problem1(lon,lat,topo)
    problem2(lon,lat,topo)
    #problem3(lon,lat,topo)
    #problem4(lon,lat,topo)    
    #problem5(lon,lat,topo)
    #problem6(lon,lat,topo)
    plt.show()


#All numerical output and answers
"""
maximum depth: -15.23506
mean depth: -6.32878514026439
mean resolution in degrees: 0.0003308
mean resolution in meters: 21.62 m
Cross section: 4699.678158801872
river_width at 54 7': 713.91 m
Transport per second: 61.096 m^3/s
Transport per year: 1926962038.67 m^3/a
Transport per year: 1.93 km^3/a
"""
