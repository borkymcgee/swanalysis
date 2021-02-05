

#!/usr/bin/env python3

####################### Swanalysis version 3.7 ###############################
##################### Created by Juno Presken ############################
# Downloads SWAN data and uses it to generate estimates of solar activity

import cartopy.img_transform as cit
import cartopy.crs as ccrs
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from scipy.signal import convolve2d
from scipy import interpolate
import requests
import resource
from mpl_toolkits.mplot3d import Axes3D
import datetime
import numpy as np
import sys
import argparse
import tempfile
import os
from multiprocessing.pool import ThreadPool, Pool


# plt.xkcd()
# Set the theme to one that matches Juno's desktop
matplotlib.rcParams.update({#'toolbar':'None',
                            'axes.facecolor':'0d0f18',
                            'figure.facecolor':'#1a1f26',
                            'savefig.facecolor':'#1a1f26',    
                            'text.color':'white',
                            'axes.labelcolor':'white',
                            'xtick.color':'white',
                            'ytick.color':'white'})

first_reading = datetime.datetime.fromisoformat("1997-01-02")

# Used to calculate progress in various tasks
program_start_time = datetime.datetime.now()
hackybefore = datetime.datetime.now()

# Where the raw data is stored
swan_array = None

version = 3.7
args = None
animate = True

###############LINES SHOULD BE NO MORE THAN THIS LENGTH########################

def main():
    global swan_array
    
    get_args()
    if args.update:
        Vprint("Updating cache file")
        updateCache(63)

    #get swan array from storage, but only the section requested
    start_num = (args.startDate - first_reading).days
    end_num = (args.endDate - first_reading).days
    swan_array = fits.getdata(args.input_file)[start_num:end_num]

    # really rudimentary code to try different masking parameters, needs to be hardcoded
    # showMasks()
    # return True
    
    if args.mask:
        maskArray()
        if args.interpolate:
            interpArray()
    if args.calculate:
        calculateArray()
    # If we've modified the array and haven't specified the -n flag, save to output_file
    if (args.mask or args.calculate) and not args.no_save:
        Vprint("Writing data to file")
        fits.writeto(args.output_file, swan_array)
    if args.display:
        dislider()

# function for testing out various options for masking and displaying them simultaneously. barely works, don't use
def showMasks():
    
    fig, axs = plt.subplots(9,9)
    
   #adjust vmin and vmax to display data well:
    maxvalue = np.max(maskSwanDay(swan_array[0]))
    Vprint("maximum in array: {}".format(maxvalue))

    testarray = swan_array[140].copy()
 
    for row, diffnum in enumerate(range(1,10)):
        for column, floodnum in enumerate(range(0,41,5)):
            Vprint("{},{}".format(diffnum, floodnum))
            plt.subplots_adjust(left = 0, right = 1, top = 1, bottom = 0,hspace = .001, wspace = .001)
            axs[row,column].set_xticks([])
            axs[row,column].set_yticks([])
            if args.interpolate:
                axs[row,column].imshow(interpSwanDay(maskSwanDay(testarray,50,diffnum,floodnum)), cmap=args.color, animated = False, aspect = 'auto', vmax=maxvalue)
            else:
                axs[row,column].imshow(maskSwanDay(testarray,50,diffnum,floodnum), cmap=args.color, animated = False, aspect = 'auto', vmax=maxvalue)
            testarray = swan_array[140].copy()
    plt.show()
            
        

# Combines data from LISIRD with the whole swan array to get an actual irradiance map for all the days
def calculateArray():
    global swan_array
    start_num = (args.startDate - first_reading).days
    end_num = (args.endDate - first_reading).days
    
    Vprint("Calculating irradiance")
    LISIRDarray = getLisirdAverage(52)
    before = datetime.datetime.now()

    # nonfunctional code to average the entire array and take the difference between that and each day of the array
    # swan_average = np.average(swan_array, 0)
    # swan_array = swan_array - swan_average

    # Run calculation in parallel for all specified dates   
    for day, array in enumerate(Pool(4).imap(calculateDay, zip(swan_array, LISIRDarray[start_num:end_num]))):
        # Code to calculate progress of operation and estimate time to completion
        swan_array[day] = array
        percent_complete = (day/(len(swan_array))*100) + 0.01
        estimated_completion = ((datetime.datetime.now() - before) / percent_complete) * (100 - percent_complete)
        Vprint("{} remaining, {}% completed".format(str(estimated_completion)[0:7], str(percent_complete)[0:4]), nl=True)
    Vprint("")
    Vprint("Done, calculating took {}".format(str(datetime.datetime.now() - before)[0:7]))

    Vprint("Calculation done")


# Calculates irradiance map for a single day, argument is a tuple of the format [swan_day, lisird irradiance]
def calculateDay(indata):
    day_array = indata[0]
    irr = indata[1]
    # Skip days where no reading was taken
    if day_array[0,0] < 0:
        return day_array
    # Code to flip array vertically to attempt to remove more noise, doesn't work at all
    # inarray[np.flip(inarray,0) <= 0] = 0
    # inarray -= np.flip(inarray,0)
    # inarray = abs(inarray)
    return irr * (1 + day_array / 0.903)


# Identify stars in array and replace all areas where stars are detected with 0s
def maskArray():
    global swan_array
    
    Vprint("Masking data")
    before = datetime.datetime.now()

    # Run masking function on all days in parallel
    for day, array in enumerate(Pool(4).imap(maskSwanDay, swan_array)):
        # calculate progress and ETA
        swan_array[day] = array
        percent_complete = (day/(len(swan_array))*100) + 0.01
        estimated_completion = ((datetime.datetime.now() - before) / percent_complete) * (100 - percent_complete)
        Vprint("{} remaining, {}% completed".format(str(estimated_completion)[0:7], str(percent_complete)[0:4]), nl=True)
    Vprint("")
    Vprint("Done, masking took {}".format(str(datetime.datetime.now() - before)[0:7]))


# Same as maskArray but for a single day. mask by...
# 1. Taking the difference of the array over both axes, and masking any values with a difference higher than threshold.
# 2. Doing this diffnum times to quickly inflate the borders of the masked region to better mask the edges of stars.
# 3. Perform a floodfill operation floodnum times, masking any pixel with more than 4 masked neighbors. This helps remove unmasked bits in the middles of stars
def maskSwanDay(day_array, threshold = 50, diffnum = 2, floodnum = 10):
    if day_array[0,0] < 0:
        return day_array

    for i in range(diffnum):
        # make a destination array for the mask
        fits_diff_bool_array = np.ones_like(day_array, bool)

        #take the difference of each axis. offset each axis in the direction of the differentiation to mask the bigger of the values
        fits_diff_bool_array[1:,:] = np.diff(day_array, axis=0) > threshold
        fits_diff_bool_array[:-1,:] = np.logical_or(fits_diff_bool_array[:-1,:], np.diff(day_array, axis=0) < -threshold)
        fits_diff_bool_array[:,1:] = np.logical_or(fits_diff_bool_array[:,1:], np.diff(day_array, axis=1) > threshold)
        fits_diff_bool_array[:,:-1] = np.logical_or(fits_diff_bool_array[:,:-1], np.diff(day_array, axis=1) < -threshold)

        for i in range(floodnum):
            fits_diff_bool_array = np.logical_or(fits_diff_bool_array, convolve2d(fits_diff_bool_array, np.ones((3,3)), 'same') > 4)
        # apply the mask to the array, replacing masked values with 0
        day_array[fits_diff_bool_array] = 0
    return day_array


# Linearly interpolate the masked out areas of the whole array
def interpArray():
    global swan_array
    
    Vprint("Interpolating data")
    before = datetime.datetime.now()
    # Interpolate all the days, utilizing multithreading
    for day, array in enumerate(Pool(4).imap(interpSwanDay, swan_array)):
        #calculate the percent complete
        swan_array[day] = array
        percent_complete = (day/(len(swan_array))*100) + 0.01
        estimated_completion = ((datetime.datetime.now() - before) / percent_complete) * (100 - percent_complete)
        Vprint("{} remaining, {}% completed".format(str(estimated_completion)[0:7], str(percent_complete)[0:4]), nl=True)
    Vprint("")
    Vprint("Done, Interpolating took {}".format(str(datetime.datetime.now() - before)[0:7]))


# Linearly interpolate the masked out areas of a single day    
def interpSwanDay(day_array):
    if day_array[0,0] < 0:
        return day_array

    # the interpolating function can't operate on the edges of the image, so we vertically flip, move down, and translate by 1/2 the image the bottom 22 rows of pixels so the interpolation algorithm knows what's on the other side

    day_array_longer = np.empty([202, 360])
    day_array_longer[22:,:] = day_array
    day_array_longer[0:22,:] = np.roll(np.flip(day_array[0:22,:],0),180,1)  

    #Voodoo magic i lifted from stackoverflow. it works great as is so i haven't really looked into how it works
    day_array_longer[day_array_longer == 0] = np.nan
    day_array_longer = np.ma.masked_invalid(day_array_longer)

    x = np.arange(0, day_array_longer.shape[1])
    y = np.arange(0, day_array_longer.shape[0])

    xx, yy = np.meshgrid(x, y)

    x1 = xx[~day_array_longer.mask]
    y1 = yy[~day_array_longer.mask]
    newarr = day_array_longer[~day_array_longer.mask]

    day_array = interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method='linear', fill_value = 0)[22:,:]    
    return day_array    
    

#Make sure cache file contains all dates needed, if no cache file is provided, generates a new one
def updateCache(threads):
    tempdir = tempfile.TemporaryDirectory()
    if (not os.path.exists(args.input_file)):
        Vprint("Generating new cache file")

        #make cache file that goes from the date of the first reading to today
        cache_start_date = first_reading
        cache_end_date = datetime.datetime.now()

        #create new array to use for cache
        new_cache_array = np.empty([(cache_end_date - cache_start_date).days + (365 * 5),180,360])
        new_cache_array[:] = np.nan
        Vprint("Array full of nans created of shape: {}".format(new_cache_array.shape))

        # Vprint("Downloading irradiance data, adding as secondaryHDU")
        # irradiance_hdu = fits.BinTableHDU.from_columns(
            # [fits.Column(name='earth irradiance', format='E', array=getLisirdArray())])

        # print(type(irradiance_hdu))
        
        #turn cache into an HDU to use with astropy
        
        cacheHDUL = fits.HDUList([fits.PrimaryHDU(new_cache_array)])
        
        Vprint("HDUList made from empty array")
    else:
        Vprint("Cache file extant, opening")
        cacheHDUL = fits.open(args.input_file, mode='update')

    cache_array = cacheHDUL[0].data
    
    if len(cache_array) < (args.endDate - first_reading).days:
        print("Hello, future astrophysicist! If you're reading this message, it means that your cache file is too small. If you want to try to fix this code, go ahead. Otherwise you can simply generate an entire new cache. Good luck from the past! -Juno")

    # iterate through the important section of the cache, replacing missing dates (which are filled with nan) with either data (if there's data)
    # or -1 (if there's no data)
    start_index = (args.startDate - first_reading).days
    end_index = (args.endDate - first_reading).days
    
    cache_changed = False

    # Generate an array of dates to download
    new_data_dates = []
    for index, (static_array, date) in enumerate(zip(cache_array[start_index:end_index], datesArray(args.startDate, args.endDate))):
        if np.isnan(static_array[0,0]):
            if args.verbose:
                sys.stdout.write("\rNo data cached for {}".format(date.date()))
                sys.stdout.flush()
            cache_changed = True
            new_data_dates.append(date)
        else:
            if args.verbose:
                sys.stdout.write("\rData cached for {}".format(date.date()))
                sys.stdout.flush()
    Vprint("")
            
    if len(new_data_dates) > 0:
        before = datetime.datetime.now()
        for day, date_and_array in enumerate(ThreadPool(threads).imap(getSwanArray, new_data_dates)):
            #add the new data to the appropriate date in the cache
            cache_array[(date_and_array[0] - first_reading).days] = date_and_array[1]
            percent_complete = (day/(len(new_data_dates))*100) + 0.01
            estimated_completion = ((datetime.datetime.now() - before) / percent_complete) * (100 - percent_complete)
            Vprint("{} remaining, {}% completed".format(str(estimated_completion)[0:7], str(percent_complete)[0:4]), nl=True)
        Vprint("")
        Vprint("Done, Downloading took {}".format(str(datetime.datetime.now() - before)[0:7]))
        
        if not os.path.exists(args.input_file):
            Vprint("Writing data to new cache file")
            cacheHDUL.writeto(args.input_file)
        else:
            Vprint("Writing data to extant cache file")
        cacheHDUL.close()
        Vprint("Cache file updated!")
    else:
        Vprint("Took {} to check cache".format(datetime.datetime.now() - startTime))
        Vprint("Cache already up-to-date!")
    del cacheHDUL[0].data

    
# Download array of lyman-alpha readings from LISIRD, from 26 days before SWAN's first reading to the present day, and return as a numpy array
def getLisirdArray():
    
    # Calculate begining and end of gliding average window (add 1/2 window size to date)
    window_begin = args.startDate - datetime.timedelta(days=26)
    window_end = args.endDate + datetime.timedelta(days=25)


    # URL request to LISIRD for the daily earth lyman-alpha irradiance average
    # for the $earth_averagin_window days preceding $window_ending_date
    lisird_url = ("http://lasp.colorado.edu/lisird/latis/dap/composite_lyman_a"
                  "lpha.tab?time,irradiance&time>{}&time<{}&format_time"
                 "(G)".format(window_begin.strftime('%Y-%m-%d'),window_end.strftime('%Y-%m-%d')))
    
    Vprint("lisird API request URL:\n{}".format(lisird_url), 2)
        
    lisird_request = requests.get(lisird_url)
    lisird_text_file = lisird_request.text 

    # Strip time, whitespace, and convert into array
    earth_irradiance_array = (lisird_text_file
                             .replace("AD\t","")
                             .strip().split())
    
    #array is still made up of strings. cast to floats
    earth_irradiance_array = [float(i) for i in earth_irradiance_array]

    return earth_irradiance_array

    
# Return an array conaining all dates between the given days
def datesArray(startDate, endDate, period = 1):
    dates = []
    thisDate = startDate
    while thisDate <= endDate:
        dates.append(thisDate)
        thisDate += datetime.timedelta(days=period)
    return dates
    

# Establish cli arguments
def get_args():
    global args
    parser = argparse.ArgumentParser(description="Generates a fits file containing an estimate of solar data for the given date(s) based on SWAN data")

    parser.add_argument("startDate",
                        nargs='?',
                        default = first_reading.strftime("%Y-%m-%d"),
                        help = "Beginning of date range. Defaults to 1997-01-02, the date of SWAN's first reading")
    parser.add_argument("endDate",
                        nargs='?',                        
                        default = (datetime.datetime.now() - datetime.timedelta(days=5)).strftime("%Y-%m-%d"),
                        help = "End of the date range. Defaults to 5 days before the present")
    parser.add_argument("-i", "--input-file",
                        dest = "input_file",
                        help = "file to perform operations on, default is a raw cache file named SWANcache.fits in the working directory")
    parser.add_argument("-o", "--output-file",
                        help = "Location to write the combined data to. Defaults to swan_processed.fits for the name and the working directory for the location if not given.")
    parser.add_argument("-u", "--update",
                        action = "store_true",
                        help = "Update cache file")
    parser.add_argument("-m", "--mask",
                        action = "store_true",
                        help = "Mask input file")
    parser.add_argument("-t", "--interpolate",
                        action = "store_true",
                        help = "interpolate input file")
    parser.add_argument("-c", "--calculate",
                        action = "store_true",
                        help = "Mask input file")
    parser.add_argument("-d", "--display",
                        action = "store_true",
                        help = "Display input file")
    parser.add_argument("-p", "--projection",
                        default = "Mercator",
                        help = "Map projection displayed, format: PROJECTION,CENTERED_LONG,CENTERED_LAT. projection options are: PlateCarree, AlbersEqualArea, AzimuthalEquidistant, EquidistantConic, LambertConformal, LambertCylindrical, Mercator, Miller, Mollweide, Orthographic, Robinson, Sinusoidal, Stereographic, TransverseMercator, UTM, InterruptedGoodeHomolosine, RotatedPole, OSGB, EuroPP, Geostationary, NearsidePerspective, EckertI, EckertII, EckertIII, EckertIV, EckertV, EckertVI, EqualEarth, Gnomonic, LambertAzimuthalEqualArea, NorthPolarStereo, OSNI, SouthPolarStereo")
    parser.add_argument("-n", "--no-save",
                        dest = "no_save",
                        action = "store_true",
                        help = "don't save the output of the calculations as a new fits file")
    parser.add_argument("-f", "--date-format",
                        dest = "date_format",
                        default = "%Y-%m-%d",
                        help = "Format code of the dates given. defaults to %%Y-%%m-%%d")
    parser.add_argument("-w", "--window",
                        type = int,
                        default = 52,
                        help = "Window size in days used to calculate the "
                               "average irradiance measured from earth. "
                               "(default 52)")
    parser.add_argument("-O", "--overwite-file",
                        dest = "overwrite",
                        action = "store_true",
                        help = "If output file already exists, delete it.")
    parser.add_argument("-C", "--color",
                        dest = "color",
                        default = "inferno",
                        help = "colormap to use for data, default is inferno")
    parser.add_argument("-v", "--verbose",
                        action = "store_true",
                        help = "increase verbosity")
    parser.add_argument("-vv", "--more-verbose",
                        action = "store_true",
                        help = "increase verbosity even more")
    parser.add_argument("-V", "--save-video",
                        dest = "save_video",
                        help = "Saves video to specified filename, implies display")



    args = parser.parse_args()

    # turn "projection" into an actual projection.
    plist = args.projection.split(',')
    lon = 0
    lat = 0
    if len(plist) > 2:
        lon = int(plist[2])
    if len(plist) > 1:
        lat = int(plist[1])

    projectionary = {"PlateCarree":ccrs.PlateCarree(central_longitude = lon),"EquidistantConic":ccrs.EquidistantConic(central_longitude = lon, central_latitude = lat),"AlbersEqualArea":ccrs.AlbersEqualArea(central_longitude = lon, central_latitude = lat),"LambertConformal":ccrs.LambertConformal(),"AzimuthalEquidistant":ccrs.AzimuthalEquidistant(central_longitude = lon, central_latitude = lat),"LambertCylindrical":ccrs.LambertCylindrical(central_longitude = lon),"Mercator":ccrs.Mercator(central_longitude = lon),"Robinson":ccrs.Robinson(central_longitude = lon),"EckertIII":ccrs.EckertIII(central_longitude = lon),"InterruptedGoodeHomolosine":ccrs.InterruptedGoodeHomolosine(central_longitude = lon),"Gnomonic":ccrs.Gnomonic(central_longitude = lon, central_latitude = lat),"Stereographic":ccrs.Stereographic(central_latitude = lat, central_longitude = lon),"EckertV":ccrs.EckertV(central_longitude = lon),"NorthPolarStereo":ccrs.NorthPolarStereo(central_longitude = lon),"EqualEarth":ccrs.EqualEarth(central_longitude = lon),"Geostationary":ccrs.Geostationary(central_longitude = lon),"SouthPolarStereo":ccrs.SouthPolarStereo(central_longitude = lon),"Sinusoidal":ccrs.Sinusoidal(central_longitude = lon),"EckertIV":ccrs.EckertIV(central_longitude = lon),"RotatedPole":ccrs.RotatedPole(),"LambertAzimuthalEqualArea":ccrs.LambertAzimuthalEqualArea(central_longitude = lon, central_latitude = lat),"TransverseMercator":ccrs.TransverseMercator(central_longitude = lon, central_latitude = lat, approx = False),"EckertVI":ccrs.EckertVI(central_longitude = lon),"OSNI":ccrs.OSNI(approx = False),"NearsidePerspectivex":ccrs.NearsidePerspective(central_longitude = lon, central_latitude = lat),"EckertI":ccrs.EckertI(central_longitude = lon),"EckertII":ccrs.EckertII(central_longitude = lon),"Miller":ccrs.Miller(central_longitude = lon),"Mollweide":ccrs.Mollweide(central_longitude = lon),"Orthographic":ccrs.Orthographic(central_longitude = lon, central_latitude = lat)}
    
    args.projection = projectionary[plist[0]]

    # If specified input_file is not provided, or is a directory, put it in the default location and/or name it with the default filename
    if args.input_file is None:
        args.input_file = os.getcwd()
    if os.path.isdir(args.input_file):
        args.input_file = os.path.join(args.input_file,"SWANcache.fits")
    Vprint("Cache file path: ", args.input_file)


    ow = False
    ow = args.overwrite
    # Do the same for output_file
    if args.output_file is None:
        args.output_file = os.getcwd()
    if os.path.isdir(args.output_file):
        fname = "SWAN"
        if args.calculate:
            fname = fname + "Calculated"
        if args.mask:
            fname = fname + "Masked"
        fname = fname + "{}_{}.fits".format(args.startDate, args.endDate)
        if args.calculate or args.mask:
            args.output_file = os.path.join(args.output_file, fname)
            if os.path.exists(args.output_file) and ow:
                os.remove(args.output_file)
                Vprint("Output filename already existed, overwriting...") 
        else:
            args.output_file = args.input_file
    else:
        if os.path.exists(args.output_file) and ow:
                os.remove(args.output_file)
                Vprint("Output filename already existed, overwriting...")       
    Vprint("Output file: {}".format(args.output_file))

        
    # turn dates into datetime objects
    args.startDate = datetime.datetime.strptime(args.startDate, args.date_format)
    args.endDate = datetime.datetime.strptime(args.endDate, args.date_format)

    #this is just here bc functions ending with comments break code folding
    pass

    
#Downloads a swan file for the specified date and returns a tuple with the date it was given and the array of the first HDU. if no file is found, returns an array of the right size filled with -1
def getSwanArray(swan_date):
    Vprint("getting data for date: {}".format(swan_date), 2)
    url = "http://swan.projet.latmos.ipsl.fr/data/L2//FSKMAPS/swan_bckgn_sm_{}_0001.fts".format(swan_date.strftime('%Y%m%d'))
    Vprint("url: {}".format(url), 2)
    # if swan_date.day == 1:
        # sys.stdout.write("\rdownloading data for {}".format(swan_date.strftime("%Y-%m")))
    swan_request = requests.get(url, stream = True)
    if swan_request.status_code != 200:
        # if args.verbose:
            # sys.stdout.write("\rNo data available for {}".format(swan_date.strftime(args.date_format)))
        empty_array = np.empty((180,360))
        empty_array.fill(-1)
        return (swan_date, empty_array)
    else:
        # if args.verbose:
            # sys.stdout.write("\rDownloading data for  {}".format(swan_date.strftime(args.date_format)))
        temp = tempfile.TemporaryFile()
        temp.write(swan_request.content)
        with fits.open(temp, ignore_missing_end=True, mode='update') as fitsFile:
            return (swan_date, fitsFile[0].data)

        
# returns an array containing a running average of irradiance data from LISIRD
def getLisirdAverage(window):

    return np.convolve(getLisirdArray(), np.ones(window)/window, mode='valid')
    for i, e in enumerate(inarray):
        outarray.append(np.mean(inarray[int(i-window/2):int(i+window/2)]))
    return outarray
 

# Display data with a slider
def dislider():
    global swan_array
    global hackybefore
    Vprint("displaying data")


    # Initialize the plot area and add gridlines
    fig, ax = plt.subplots(subplot_kw={'projection': args.projection})
    ax.set_global()
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

    before = datetime.datetime.now()
    Vprint("Processing data for display:")
    start_num = (args.startDate - first_reading).days
    end_num = (args.endDate - first_reading).days

    # We're going to remove layers that don't have a reading from the data,
    # so we need a list of dates to keep track of when each layer is from
    # Get an array of all dates between the start and end
    dates_array = np.array(datesArray(args.startDate, args.endDate))[:len(swan_array)]

    Vprint("\t- Removing days with no readings", True)
    # remove days for which no reading exists
    dates_array = dates_array[swan_array[:,0,0] >= 0]
    swan_array = swan_array[swan_array[:,0,0] >= 0]
    Vprint(" (took {})".format(str(datetime.datetime.now() - before)[0:7]))

    Vprint("Reprojecting data")
    date_title = plt.title(dates_array[0].date())

    before = datetime.datetime.now()
    #using multiple processes, reproject each day in swan_array, and also pass the function the percent done
    for day, array in enumerate(Pool(4).imap(reproject, swan_array)):
        swan_array[day] = array
        percent_complete = (day/len(swan_array)*100) + 0.01
        estimated_completion = ((datetime.datetime.now() - before) / percent_complete) * (100 - percent_complete)
        Vprint("{} remaining, {}% completed".format(str(estimated_completion)[0:7], str(percent_complete)[0:4]), nl=True)
    Vprint("")
    Vprint("Done, reprojection took {}".format(str(datetime.datetime.now() - before)[0:7]))


    #adjust vmin and vmax to display data well:
    maxvalue = np.max(maskSwanDay(swan_array[0]))
    Vprint("maximum in array: {}".format(maxvalue))
    
    # plot holding the data
    frame_image = ax.imshow(swan_array[0], cmap=args.color, animated = True, aspect = 'auto', transform=ccrs.Mercator(), vmax=maxvalue)
    
    slider_axis = plt.axes([0.125, 0.03, 0.775, 0.02])
    day_slider = Slider(slider_axis, 'Date', 0, len(swan_array) - 1, valinit=0, valstep=1, valfmt = '%.0s')

        
    # Handle keypresses
    def on_key_press(event):
        global animate
        if event.key == 'right':
            day_slider.set_val((day_slider.val + 1) % len(swan_array))
        elif event.key == 'left':
            day_slider.set_val((day_slider.val - 1) % len(swan_array))
        elif event.key == ',':
            day_slider.set_val((day_slider.val - 10) % len(swan_array))
        elif event.key == '.':
            day_slider.set_val((day_slider.val + 10) % len(swan_array))
        elif event.key == ' ' or args.save_video is not None:
            animate = not animate
            

    # Redraw image when slider is changed
    def sliderUpdate(day):
        global args
        global hackybefore
        if args.save_video or animate:
            if day_slider.val != day:
                day_slider.set_val(day % len(swan_array))
        else:
            day = day_slider.val
        if args.save_video is not None and day >= len(swan_array):
            Vprint("displaying data took {}".format(str(datetime.datetime.now() - before)[0:7]))
            hackybefore = datetime.datetime.now()
            Vprint("saving video...")
            ani.save(args.save_video, savefig_kwargs={'facecolor':fig.get_facecolor()})
            Vprint("saved! saving took {}".format(str(datetime.datetime.now() - hackybefore)[0:7]))
            Vprint("program done, took {} to complete".format(str(datetime.datetime.now() - program_start_time)[0:7]))            
            args.save_video = None
        elif hackybefore is not None and day >= len(swan_array):
            Vprint("program done, took {} to complete".format(str(datetime.datetime.now() - program_start_time)[0:7]))
            hackybefore = None
        day = day % len(swan_array)
        date_title.set_text("SWAN data for {}".format(dates_array[int(day)].date()))        
        frame_image.set_data(swan_array[int(day)])
        fig.canvas.draw_idle
        return 1,

    day_slider.on_changed(sliderUpdate)
    fig.canvas.mpl_connect('key_press_event', on_key_press)
    ani = animation.FuncAnimation(fig, sliderUpdate, interval = 50, save_count = len(swan_array))

    os.system("mpg123 -q --no-control notify.mp3")
    hackybefore = datetime.datetime.now()
    plt.show()


# Reproject the given array into the projection specified by args.projection
def reproject(source_array):
    if isinstance(args.projection, ccrs.Mercator):
        return source_array
    else:
        return cit.warp_array(source_array, args.projection, ccrs.Mercator(central_longitude = 180), target_res=(360,180), mask_extrapolated=True)[0]


#Vprint if verbosity is enabled
def Vprint(text = "", importance = 1, nl = False):
    if (args.verbose and importance == 1) or (args.more_verbose and importance == 2):
        if nl:
            sys.stdout.write("\r" + str(text))
            sys.stdout.flush()
        else:
            print(text)
        

# Limit memory to percent of free memory
def memory_limit(percent):
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    print("memory limited to {} / {}".format(int(get_memory() * (percent / 100)) / 1048576, get_memory() / 1048576))
    resource.setrlimit(resource.RLIMIT_AS, (int(get_memory() * 1024 * 0.75), hard))

# get amount of free memory, probably only works on POSIX systems
def get_memory():
    with open('/proc/meminfo', 'r') as mem:
        free_memory = 0
        for i in mem:
            sline = i.split()
            if str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
                free_memory += int(sline[1])
    return free_memory

# memory_limit(99)
try:
    main()
except MemoryError:
    sys.stderr.write('\n\nERROR: Memory Exception\n')
    sys.exit(1)
