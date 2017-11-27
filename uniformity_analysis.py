# work on python 2.7 and python 3.6

# test data provided at https://drive.google.com/file/d/0B99lmHE3zUdWQXZPb1laY0Rsbnc/view?usp=sharing

import numpy as np
import re
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import json
import sys

def main():

    print (sys.argv[1])
    parameters_file = sys.argv[1]

    parameters = {}
    with open(parameters_file) as json_file:
        parameters = json.load(json_file)
    #print(json.dumps(parameters, sort_keys=True, indent=4))


    scale = parameters["scale_factor"]*parameters["arcsec_per_pixel"]
    parameters["title_gaussian_fit"] = 'my gaussian plot '+str(parameters["gaussian_start_pixel"])+ ' to ' +str(parameters["gaussian_end_pixel"])

    if parameters["perform_slice_analysis"] == True and parameters["scaled"] == True:
        parameters["scaled"] = False
        print("scaled set to false, cant perform slice analysis with scaled=True if desired rerun separatly with set perform_slice_analysis=False scaled=True")

    # start of analysis
    if parameters["find_center_of_intensity"]:
        if parameters["mask_center"]:
            print('cannot mask center while finding center of mass')
            return 0

    f = open(parameters['data_file'],'r')

    keys = f.readline()

    if (parameters['pix_widthX']%2==0):
        parameters['pix_widthX'] = parameters['pix_widthX'] -1
    if (parameters['pix_widthY']%2==0):
        parameters['pix_widthY'] = parameters['pix_widthY'] -1

    printable_info = ''
    printable_info = printable_info + 'Bin size ' + str(parameters['bin_size']) + ' pixels; map pixel center ' + str(parameters['centerX']) + ', ' + str(parameters['centerY']) + '; width in pixels '+ str(parameters['pix_widthX']) + ' height in pixels ' + str(parameters['pix_widthY'])

    print('Bin size: ' + str(parameters['bin_size']))

    print ("Center (%d, %d) width %d height %d " %(parameters['centerX'],parameters['centerY'],parameters['pix_widthX'],parameters['pix_widthY']))

    half_widthX = (parameters['pix_widthX'] -1 ) / 2
    half_widthY = (parameters['pix_widthY'] -1 ) / 2

    max_size = half_widthX+1 # the lesser of x or y size
    if half_widthY < half_widthX:
        max_size = half_widthY+1

    parameters['circle_levels'] = []
    for index in range(int(max_size / parameters['bin_size'])):
        parameters['circle_levels'].append(index*parameters['bin_size']+parameters['bin_size'])

    startX = parameters['centerX']-half_widthX
    endX = parameters['centerX']+half_widthX
    startY = parameters['centerY']-half_widthY
    endY = parameters['centerY']+half_widthY

    parameters["scaled_startX"] = startX
    parameters["scaled_endX"] = endX
    parameters["scaled_startY"] = startY
    parameters["scaled_endY"] =endY
    if  parameters["scaled"]:
        parameters["scaled_startX"] = -half_widthX * scale
        parameters["scaled_endX"] = half_widthX * scale
        parameters["scaled_startY"] = -half_widthY * scale
        parameters["scaled_endY"] = half_widthY * scale

    x = []
    y = []
    intensities = []
    miniZ = []
    fullZ = []
    count = 1

    parameters['high_values'] = []
    total_pixels_available = 0

    for line in f:
        array = parseLine(line)
        x_value = float(array[0])
        y_value = float(array[1])
        intensity_value = float(array[2])
        if x_value >= startX and x_value <= endX and y_value >= startY and y_value <= endY:
            x.append(x_value)
            y.append(y_value)
            if intensity_value == 0:
                print('intensity_value = 0')
            if parameters['mask_center'] and ((x_value-parameters['centerX'])**2+(y_value-parameters['centerY'])**2)/max_size**2 < parameters['mask_inner_percentage']:
                intensity_value = 0
            angle = np.arctan2(y_value - parameters['centerY'], x_value - parameters['centerX']) # returns angle -Pi to -Pi in radians , y_value must come first
            if angle < 0:
                angle = angle + 2*np.pi
            if parameters['scoped'] and (angle < scoped_angle_start*np.pi/180 or angle >= scoped_angle_end*np.pi/180):
                intensity_value = 0
            if intensity_value > 0:
                total_pixels_available = total_pixels_available + 1
            if (intensity_value > parameters['high_threshold']):
                parameters['high_values'].append(intensity_value)
            intensities.append(intensity_value)
            miniZ.append(intensity_value)

            if (count%parameters['pix_widthX'] == 0):
                fullZ.append(miniZ)
                miniZ = []
            count = count + 1

    if parameters['find_center_of_intensity']:
        find_center_of_intensity_function(x, y, intensities, parameters['centerX'], parameters['centerY'], parameters['inner_percentage'], max_size);
        print ('center of intensity determined with ' + str(int(max_size *parameters['inner_percentage'])) + ' pixels')
        print ("throw_negatives " + str(parameters['throw_negatives']))
        print ("throw_high_values " + str(parameters['throw_high_values']))

    if parameters['throw_high_values']:
        print (parameters['high_values'])
        print ("%d high pixel values thrown out. %2.2f%% of available pixels" %(len(parameters['high_values']), len(parameters['high_values'])/total_pixels_available*100))
    else:
        print ('No high pixel values thrown')

    if parameters['print_contour']:
        xlist = np.linspace(parameters["scaled_startX"], parameters["scaled_endX"], parameters['pix_widthX'])
        ylist = np.linspace(parameters["scaled_startY"], parameters["scaled_endY"], parameters['pix_widthY'])
        X, Y = np.meshgrid(xlist, ylist)
        Z=fullZ

        fig = plt.figure()
        if parameters['print_custom_intensity_levels']:
            contour = plt.contour(X, Y, Z, parameters['levels'], colors='k')
            plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
            contour_filled = plt.contourf(X, Y, Z, parameters['levels'])
        else:
            contour = plt.contour(X, Y, Z, colors='k')
            plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
            contour_filled = plt.contourf(X, Y, Z)

        plt.colorbar(contour_filled)

        x_circles = np.linspace(parameters["scaled_startX"], parameters["scaled_endX"], parameters['pix_widthX'])
        y_circles = np.linspace(parameters["scaled_startY"], parameters["scaled_endY"], parameters['pix_widthY'])
        for index in range(len(parameters['circle_levels'])):
            parameters['circle_levels'][index] = parameters['circle_levels'][index]*parameters['circle_levels'][index]
        X, Y = np.meshgrid(x_circles, y_circles)
        a = (X - parameters['centerX']) ** 2 + (Y - parameters['centerY']) ** 2
        if parameters['scaled']:
            a = X ** 2 + Y ** 2
        if parameters['print_circles']:
            c = plt.contour(x_circles, y_circles, a, parameters['circle_levels'], colors='k')

        plt.title(parameters['title_of_contour'])
        plt.xlabel('x (pixel)')
        plt.ylabel('y (pixel)')
        if parameters['scaled']:
            plt.xlabel('delta mas')
            plt.ylabel('delta mas')
        plt.savefig('contour.png')

    # histagram
    max_size = half_widthX+1
    if half_widthY < half_widthX:
        max_size = half_widthY+1

    master_array = []
    area_array = []
    mas_array = []
    pixel_array = []

    number_of_arrays = int(max_size / parameters['bin_size']+2) # we leave on zero so bin_radius can be the key in master_array, then an extra one to push partial boxes to a higher bin
    for index in range(number_of_arrays):
        master_array.append(0)
        area_array.append(0)
        mas_array.append(parameters['bin_size']*index*parameters['arcsec_per_pixel']*1000)
        pixel_array.append(parameters['bin_size']*index)

    for index in range(len(intensities)):
        if np.sqrt((x[index] - parameters['centerX']) ** 2 + (y[index] - parameters['centerY']) ** 2) > max_size:
            continue
        if parameters['throw_negatives'] and intensities[index] < 0:
            continue
        if (parameters['throw_high_values'] and intensities[index] > parameters['high_threshold']):
            continue
        if intensities[index] == 0: # from masking the center or scoping with angles, includes naturally zero points
            continue
        radius = np.sqrt((x[index] - parameters['centerX']) ** 2 + (y[index] - parameters['centerY']) ** 2) + 1 #plus one is to say the first pixel is of size one
        bin_radius = int(radius / parameters['bin_size'] + .99999) # decides to push partial pixels to a higher bin
        area = 1
        if parameters['scaled']:
            area = area *(parameters['arcsec_per_pixel']*1000)**2 # pixel^2 to mas^2
        area_array[bin_radius] = area_array[bin_radius] +area
        master_array[bin_radius] = master_array[bin_radius] + intensities[index]
        # master_array[bin_radius] = master_array[bin_radius] + 1 # for testing
    for index in range(len(area_array)):
        if area_array[index] > 0:
            master_array[index] = master_array[index]/area_array[index]

    if parameters['throw_zero_bin']:
        master_array = master_array[1:]
        mas_array = mas_array[1:]
        pixel_array = pixel_array[1:]

    plt.figure()
    # plt.plot(range(number_of_arrays)[1:],master_array) #prints bin number #depends on parameters['throw_zero_bin'] # if thrown zero bin = true there is a off by one error for the x axis, maybe I fixed it
    # plt.xlabel('bin number')

    if parameters['print_radial']:
        plt.title(parameters['title_of_radial_plot'])
        plt.xlabel('pixel')
        plt.ylabel('millijy/beam/area (mJy/beam/pix^2)')
        if parameters['scaled']:
            plt.plot(mas_array,master_array)
            plt.ylabel('millijy/beam/area (mJy/beam/mas^2)')
            plt.xlabel('mas')
        else:
            plt.plot(pixel_array,master_array)
        plt.savefig('radial.png')

    parameters['high_value'] = 0
    high_bin = 0
    for index in range(len(master_array)):
        if master_array[index] > parameters['high_value']:
            parameters['high_value'] = master_array[index]
            high_bin = index

    print ("High bin: " + str(high_bin))
    print ("High bin value: " +str(master_array[high_bin]))
    print ("To pixel " + str(high_bin * parameters['bin_size']))
    print ("To mas %2.2f"% (parameters['arcsec_per_pixel']*high_bin * parameters['bin_size']*1000))
    # print ("physical size" , parameters['arcsec_per_pixel']*high_bin * parameters['bin_size']/3600 *np.pi/180 * 2.4*1000*206264.807497, "a.u." )
    # print ("physical size" , parameters['arcsec_per_pixel']*high_bin * parameters['bin_size']/3600 *np.pi/180 * 2.4*1000, "pc" )
    print ("Physical size %2.3f*10^16 cm" %(parameters['arcsec_per_pixel']*high_bin * parameters['bin_size']/3600*np.pi/180 * parameters['radius_from_source']*1000 *3085677600000000000 /10000000000000000) )

    if parameters['perform_gaussian_analysis']:
        start = int(parameters['gaussian_start_pixel'] / parameters['bin_size']) - 1
        end = int(parameters['gaussian_end_pixel'] / parameters['bin_size'])
        sub_array_pixels = pixel_array[start:end]
        sub_array_intensities = master_array[start:end]

        mean = 0
        length = len(sub_array_pixels)
        for index in range(length):
            mean = mean + sub_array_pixels[index]/length
        my_sum = 0
        for index in range(length):
            my_sum = my_sum + (sub_array_pixels[index] - mean)**2 / length
        stddev = np.sqrt(my_sum)


        # plot the data
        if parameters["print_gaussian_fit"]:
            plt.figure()
            plt.plot(sub_array_pixels,sub_array_intensities, 'b+:', label='data')

        # fit the data
        popt, pcov = curve_fit(gaus2, sub_array_pixels,sub_array_intensities, p0=[300, mean, stddev, 1])
        # print ("popt ", popt)

        # guassian statistic are contained in the popt array
        a = popt[0]
        x0 = popt[1]
        sigma = math.fabs(popt[2])
        c = popt[3]

        if popt[2] < 0:
            print("Warning negative sigma value", popt[2])

        # approx solution to FWHM
        fwhm = 2.3548 * sigma

        # increase resolution to plot the fit
        high_resolution_array = np.linspace(parameters['gaussian_start_pixel'], parameters['gaussian_end_pixel'], 100)

        #plot the higher resolution data
        if parameters["print_gaussian_fit"]:
            plt.plot(high_resolution_array, gaus2(high_resolution_array, *popt))

        # floor the gaussian fit
        popt[3] = 0

        # plot the higher resolution data as a floored function
        if parameters["print_gaussian_fit"]:
            plt.plot(high_resolution_array, gaus2(high_resolution_array, *popt))
            plt.legend()
            plt.title(parameters['title_gaussian_fit'])
            plt.xlabel('pixel')
            plt.ylabel('milliJy/beam/area (mJy/beam/mas^2)')
            plt.savefig('gaussian.png')

        # move the guassian statistic to common variable names
        parameters['high_value'] = gaus2(x0, *popt)
        parameters['high_index'] = x0
        parameters['width'] = fwhm

        print('high value: '+ str(parameters['high_value']) + ' at ' + str(parameters['high_index']))
        print ("Physical size %2.3f*10^16 cm" %(parameters['arcsec_per_pixel']*parameters['high_index'] /3600*np.pi/180 * parameters['radius_from_source']*1000 *3085677600000000000 /10000000000000000) )
        print('Pixel width ' + str(parameters['width']))
        print('mas width ' + str(parameters['width']*parameters['arcsec_per_pixel']*1000))
        print ("Gaussian Physical width %2.3f*10^16 cm" %(parameters['arcsec_per_pixel'] * parameters['width'] / 3600 * np.pi / 180 * parameters['radius_from_source'] * 1000 * 3085677600000000000 / 10000000000000000) )

        # gather all the statistics to print them to screen and file
        if parameters['throw_high_values']:
            printable_info = printable_info + str(len(parameters['high_values'])) + ' high pixel values thrown out. ' + str(len(parameters['high_values'])/total_pixels_available*100) +' of available pixels '
        else:
            printable_info = printable_info + '\nNo high pixel values thrown'

        printable_info = printable_info + "\nHigh bin " + str(high_bin) + '\nAngular radius ' + '%3.3f' %(parameters['arcsec_per_pixel'] * high_bin * parameters['bin_size']*1000) + ' mas; Physical radius '+ '%3.3f' %(parameters['arcsec_per_pixel'] * high_bin * parameters['bin_size'] / 3600 * np.pi / 180 * parameters['radius_from_source']*1000 * 3085677600000000000 / 10000000000000000) + '*10^16 cm'
        printable_info = printable_info + "\nGaussian high pixel " + str(parameters['high_index']) + '\nAngular radius ' + '%3.3f' %(parameters['arcsec_per_pixel'] * parameters['high_index'] *1000) + ' mas; Physical radius '+ '%3.3f' %(parameters['arcsec_per_pixel'] * parameters['high_index'] / 3600 * np.pi / 180 * parameters['radius_from_source']*1000 * 3085677600000000000 / 10000000000000000) + '*10^16 cm'
        printable_info = printable_info + '\nAngular thickness ' + '%3.3f' %(parameters['width'] * parameters['arcsec_per_pixel'] * 1000) + ' mas; Physical thickness ' + '%3.3f' % (parameters['arcsec_per_pixel']*parameters['width']/3600*np.pi/180 * parameters['radius_from_source'] * 1000 * 3085677600000000000 / 10000000000000000) + '*10^16 cm'
        print (printable_info)

        if parameters['perform_slice_analysis']:
            # These values will be over written
            # parameters['scoped'] = scoped
            # parameters['scoped_angle_start'] = scoped_angle_start
            # parameters['scoped_angle_end'] = scoped_angle_end

            # needed for execution for contour slices
            parameters["array_results"] = []

            # number of nonsymmetrical slices
            bad_sym = 0.0

            number_of_cuts = int(360 / parameters['angle_slice'])

            if parameters["scoped"]:
                number_of_cuts = int( (scoped_angle_end - scoped_angle_start) / parameters['angle_slice'])
            else:
                scoped_angle_start = 0
                scoped_angle_end = 360

            parameters['scoped'] = True

            print ('number of cuts ' + str(number_of_cuts))

            plt.figure()

            for ind in range(number_of_cuts):
                parameters['scoped_angle_start'] = scoped_angle_start + ind*parameters['angle_slice']
                parameters['scoped_angle_end'] = parameters['scoped_angle_start'] + parameters['angle_slice']

                # function that runs the same analysis on an angular slice
                print ("angle: " +str(parameters['scoped_angle_start']))
                results = pie_slice(parameters)

                # uses for contour slices
                results['scoped_angle_start'] = parameters['scoped_angle_start']
                results['scoped_angle_end'] = parameters['scoped_angle_end']
                parameters["array_results"].append(results)

                # extract data into variables
                slice_fwhm = results['fwhm']
                slice_peak_pixel = results['peak_index']

                # change data based on bad recieved variables
                if results['fwhm'] > 3 * parameters['width']:
                    print ('throwing ' + str(results['fwhm']))
                    slice_fwhm = 3 * parameters['width']

                print ("pixels separation from overall fit to slice fit: " + str(math.fabs(slice_peak_pixel - parameters['high_index'])))
                print ("Gaussian peak value: " + str(results['a']))
                print ("fwhm: " + str(results['fwhm']))

                if math.fabs(slice_peak_pixel - parameters['high_index']) > parameters['width'] / 2.0 or results['a'] < 0 or results['fwhm'] > 3 * parameters['width']:
                    bad_sym = bad_sym + 1

            if number_of_cuts == 0:
                print("no slcies")
                number_of_cuts=1
            uniformity = (1.0-bad_sym/number_of_cuts)*100
            print ("Uniformity: " + str(uniformity) + '%')

            slice_contour(parameters, X, Y, Z)

            printable_info = printable_info + "\nAngle Slice Thickness " + str(parameters['angle_slice']) + "\nNumber of Slices " + str(number_of_cuts) + "\nUniformity " + str(uniformity) + "%"

        doc = 'info.txt'
        print (doc)
        g = open(doc,'w')
        g.write(printable_info)

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def gaus2(x,a,x0,sigma,c):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+c

def slice_contour(parameters, X, Y, Z):
    fig = plt.figure()
    if parameters["print_custom_intensity_levels"]:
        contour = plt.contour(X, Y, Z, parameters["levels"], colors='k')
        plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
        contour_filled = plt.contourf(X, Y, Z, parameters["levels"])
    else:
        contour = plt.contour(X, Y, Z, colors='k')
        plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
        contour_filled = plt.contourf(X, Y, Z)

    plt.colorbar(contour_filled)

    circle_levels = []
    circle_levels.append(parameters["high_index"]-parameters["width"])
    circle_levels.append(parameters["high_index"])
    circle_levels.append(parameters["high_index"]+parameters["width"])

    x_circles = np.linspace(parameters["scaled_startX"], parameters["scaled_endX"], parameters["pix_widthX"])
    y_circles = np.linspace(parameters["scaled_startY"], parameters["scaled_endY"], parameters["pix_widthY"])
    for index in range(len(circle_levels)):
        circle_levels[index] = circle_levels[index]*circle_levels[index] * circle_levels[index]/math.fabs(circle_levels[index])

    X, Y = np.meshgrid(x_circles, y_circles)
    a = (X - parameters["centerX"]) ** 2 + (Y - parameters["centerY"]) ** 2
    if parameters["scaled"]:
        a = X ** 2 + Y ** 2

    c = plt.contour(x_circles, y_circles, a, circle_levels, colors='y')

    for results in parameters["array_results"]:
        scoped_angle_start = results['scoped_angle_start']
        scoped_angle_end = results['scoped_angle_end']
        fwhm = results['fwhm']
        peak_index = results['peak_index']
        a = results['a']
        width = parameters["width"] # the fwhm of the entire image
        # peak_value = results['peak_value']
        slice_fwhm = fwhm # changes if too large to draw
        slice_peak_pixel = results['peak_index']

        color = "#66ff00" # green

        if peak_index == 0:
            print ('peak_index is zero')
            color = "red"
        elif  fwhm == 0:
            print ('fwhm is zero')
            color = "red"
        elif math.fabs(slice_peak_pixel - parameters["high_index"]) > 5 * width / 2.0:
            print ('outside 2.5 fwhm')
            color = "red"
        elif math.fabs(slice_peak_pixel - parameters["high_index"]) > width / 2.0 or results['a'] < 0 or results['fwhm'] > 3 * width:
            color = "red"

        # plots the lines that are dividers for the angular slices
        r = parameters['pix_widthX']/2.0 - 1
        simple_angle = to_radians(get_simple_angle(scoped_angle_start))

        x_first_offset = parameters["centerX"]
        x_second_offset = x_sign(scoped_angle_start) * r * math.cos(simple_angle)
        y_first_offset = parameters["centerY"]
        y_second_offset = y_sign(scoped_angle_start) * r * math.sin(simple_angle)

        x1 = x_first_offset + x_sign(scoped_angle_start) * (r-7) * math.cos(simple_angle)
        y1 = y_first_offset + y_sign(scoped_angle_start) * (r-7) * math.sin(simple_angle)

        x2 = x_first_offset + x_second_offset
        y2 = y_first_offset + y_second_offset

        plt.plot([x1, x2],[y1, y2], color)

        #plot the end angle in the range
        offset_angle = .5
        simple_angle = to_radians(get_simple_angle(scoped_angle_end-offset_angle))

        x_first_offset = parameters["centerX"]
        x_second_offset = x_sign(scoped_angle_end-offset_angle) * r * math.cos(simple_angle)

        y_first_offset = parameters["centerY"]
        y_second_offset = y_sign(scoped_angle_end-offset_angle) * r * math.sin(simple_angle)

        x1 = x_first_offset + x_sign(scoped_angle_end-offset_angle) * (r-7) * math.cos(simple_angle)
        y1 = y_first_offset + y_sign(scoped_angle_end-offset_angle) * (r-7) * math.sin(simple_angle)

        x2 = x_first_offset + x_second_offset
        y2 = y_first_offset + y_second_offset

        plt.plot([x1, x2],[y1, y2], color)

        if color == "red":
            continue

        # plots for the gaussian fits on the contour plot
        middle_angle = (scoped_angle_start + scoped_angle_end) / 2.0
        simple_angle = to_radians(get_simple_angle(middle_angle))

        x_first_offset = parameters["centerX"]
        x_second_offset = x_sign(middle_angle) * slice_peak_pixel * math.cos(simple_angle)
        x_third_offset = x_sign(middle_angle) * slice_fwhm * math.cos(simple_angle)

        y_first_offset = parameters["centerY"]
        y_second_offset = y_sign(middle_angle) * slice_peak_pixel * math.sin(simple_angle)
        y_third_offset = y_sign(middle_angle) * slice_fwhm * math.sin(simple_angle)

        x1 = x_first_offset + x_second_offset - x_third_offset
        y1 = y_first_offset + y_second_offset - y_third_offset

        x2 = x_first_offset + x_second_offset + x_third_offset
        y2 = y_first_offset + y_second_offset + y_third_offset

        # correct lengths if too long
        if math.fabs(y_second_offset + y_third_offset) > r or math.fabs(x_second_offset + x_third_offset) > r:
            y_first_offset = parameters["centerY"]
            y_second_offset = y_sign(middle_angle) * (r+0) * math.sin(simple_angle)
            y2 = y_first_offset + y_second_offset

            x_first_offset = parameters["centerX"]
            x_second_offset = x_sign(middle_angle) * (r+0) * math.cos(simple_angle)
            x2 = x_first_offset + x_second_offset

            y1 = y_first_offset + y_sign(middle_angle) * 10 * math.sin(simple_angle)
            x1 = x_first_offset + x_sign(middle_angle) * 10 * math.cos(simple_angle)

        plt.plot([x1, x2],[y1, y2], color)

        new_angle = middle_angle + 90
        simple_angle = to_radians(get_simple_angle(new_angle))
        x_forth_offset = 2*math.cos(simple_angle) * parameters["arcsec_per_pixel"] *1000 * x_sign(new_angle) #  first direction at 90 angle doesn't matter
        y_forth_offset = 2*math.sin(simple_angle) * parameters["arcsec_per_pixel"] *1000 * y_sign(new_angle)

        # close error looking cap
        x1b = x1 - x_forth_offset
        y1b = y1 - y_forth_offset

        x2b = x1 + x_forth_offset
        y2b = y1 + y_forth_offset

        plt.plot([x1b,x2b], [y1b,y2b], color)

        # far error looking cap
        x1b = x2 - x_forth_offset
        y1b = y2 - y_forth_offset

        x2b = x2 + x_forth_offset
        y2b = y2 + y_forth_offset

        plt.plot([x1b,x2b], [y1b,y2b], color)

    plt.title(parameters["title_of_uniformity_plot"])
    plt.xlabel('x (pixel)')
    plt.ylabel('y (pixel)')
    if parameters["scaled"]:
        plt.xlabel('delta mas')
        plt.ylabel('delta mas')
    plt.savefig('contour_slices.png')

def get_simple_angle(angle):
    # handles picking the right angle 0 < angle < 90 for the trig functions
    if angle < 0:
        angle = angle + 360
    if angle > 360:
        angle = angle - 360
    simple_angle = angle
    if angle >= 90 and angle <= 180:
        simple_angle = 180 - angle
    elif angle >= 180 and angle <= 270:
        simple_angle = angle - 180
    elif angle >= 270 and angle <= 360:
        simple_angle = 360 - angle

    return simple_angle

def to_radians(angle):
    return angle*np.pi/180

def x_sign (angle):
    if angle >= 90 and angle <= 270:
        return -1
    return 1

def y_sign (angle):
    if angle >= 180 and angle <= 360:
        return -1
    return 1

def pie_slice(parameters):

    master_high_index = parameters['high_index']
    master_width = parameters['width']

    bin_size = parameters['bin_size']
    pix_widthX = parameters['pix_widthX']
    pix_widthY = parameters['pix_widthY']
    centerX = parameters['centerX']
    centerY = parameters['centerY']
    radius_from_source = parameters['radius_from_source']
    scaled = parameters['scaled']
    throw_zero_bin = parameters['throw_zero_bin']
    mask_center = parameters['mask_center']
    mask_inner_percentage = parameters['mask_inner_percentage']
    throw_negatives = parameters['throw_negatives']
    scoped = parameters['scoped'] # this is forced to True

    scoped_angle_start = parameters['scoped_angle_start']
    scoped_angle_end = parameters['scoped_angle_end']

    throw_high_values = parameters['throw_high_values']
    high_threshold = parameters['high_threshold']
    perform_gaussian_analysis = parameters['perform_gaussian_analysis']
    gaussian_start_pixel = parameters['gaussian_start_pixel']
    gaussian_end_pixel = parameters['gaussian_end_pixel']
    title_gaussian_fit = parameters['title_gaussian_fit']
    arcsec_per_pixel = parameters['arcsec_per_pixel']
    scale = parameters['scale']

    title_of_contour = parameters['title_of_contour']
    title_of_radial_plot = parameters['title_of_radial_plot']
    print_circles = parameters['print_circles']
    print_custom_intensity_levels = parameters['print_custom_intensity_levels']
    levels = parameters['levels']
    print_contour = parameters['print_contour']
    print_radial = parameters['print_radial']
    find_center_of_intensity = parameters['find_center_of_intensity']
    data_file = parameters['data_file']

    f = open(data_file,'r')

    keys = f.readline()

    scoped_angle_start = scoped_angle_start *np.pi/180
    scoped_angle_end = scoped_angle_end *np.pi/180

    if (pix_widthX%2==0):
        pix_widthX = pix_widthX -1
    if (pix_widthY%2==0):
        pix_widthY = pix_widthY -1

    half_widthX = (pix_widthX -1 ) / 2
    half_widthY = (pix_widthY -1 ) / 2

    max_size = half_widthX+1 # the lesser of x or y size
    if half_widthY < half_widthX:
        max_size = half_widthY+1

    startX = centerX-half_widthX
    endX = centerX+half_widthX
    startY = centerY-half_widthY
    endY = centerY+half_widthY

    scaled_startX = startX
    scaled_endX = endX
    scaled_startY = startY
    scaled_endY =endY
    if scaled:
        scaled_startX = -half_widthX * scale
        scaled_endX = half_widthX * scale
        scaled_startY = -half_widthY * scale
        scaled_endY = half_widthY * scale

    x = []
    y = []
    intensities = []
    miniZ = []
    count = 1

    high_values = []
    total_pixels_available = 0

    for line in f:
        array = parseLine(line)
        x_value = float(array[0])
        y_value = float(array[1])
        intensity_value = float(array[2])
        if x_value >= startX and x_value <= endX and y_value >= startY and y_value <= endY:
            x.append(x_value)
            y.append(y_value)
            if intensity_value == 0:
                print('intensity_value = 0')
            if mask_center and ((x_value-centerX)**2+(y_value-centerY)**2)/max_size**2 < mask_inner_percentage:
                intensity_value = 0
            angle = np.arctan2(y_value - centerY, x_value - centerX) # returns angle -Pi to -Pi in radians , y_value must come first
            if angle < 0:
                angle = angle + 2*np.pi
            if scoped and (angle < scoped_angle_start or angle >= scoped_angle_end):
                intensity_value = 0
            if intensity_value > 0:
                total_pixels_available = total_pixels_available + 1
            if (intensity_value > high_threshold):
                high_values.append(intensity_value)
            intensities.append(intensity_value)

            count = count + 1

    max_size = half_widthX+1
    if half_widthY < half_widthX:
        max_size = half_widthY+1

    master_array = []
    area_array = []
    mas_array = []
    pixel_array = []
    number_of_circle_pixels = 0

    number_of_arrays = int(max_size / bin_size+2) # we leave on zero so bin_radius can be the key in master_array, then an extra one to push partial boxes to a higher bin
    for index in range(number_of_arrays):
        master_array.append(0)
        area_array.append(0)
        mas_array.append(bin_size*index*arcsec_per_pixel*1000)
        pixel_array.append(bin_size*index)

    for index in range(len(intensities)):
        if np.sqrt((x[index] - centerX) ** 2 + (y[index] - centerY) ** 2) > max_size:
            continue
        if throw_negatives and intensities[index] < 0:
            continue
        if (throw_high_values and intensities[index] > high_threshold):
            continue
        if intensities[index] == 0: # from masking the center or scoping with angles, includes naturally zero points
            continue
        radius = np.sqrt((x[index] - centerX) ** 2 + (y[index] - centerY) ** 2) + 1 #plus one is to say the first pixel is of size one
        bin_radius = int(radius / bin_size + .99999) # decides to push partial pixels to a higher bin
        area = 1
        if scaled:
            area = area *(arcsec_per_pixel*1000)**2 #pixel^2 to mas^2
        area_array[bin_radius] = area_array[bin_radius] +area
        master_array[bin_radius] = master_array[bin_radius] + intensities[index]
        if intensities[index]  > 24:
            number_of_circle_pixels = number_of_circle_pixels + 1
        # master_array[bin_radius] = master_array[bin_radius] + 1 # for testing
    for index in range(len(area_array)):
        if area_array[index] > 0:
            master_array[index] = master_array[index]/area_array[index]

    if throw_zero_bin:
        master_array = master_array[1:]
        mas_array = mas_array[1:]
        pixel_array = pixel_array[1:]

    fwhm = 0
    peak_value = 0
    peak_index = 0
    a = -1

    try:
        if perform_gaussian_analysis:
            start = int(gaussian_start_pixel / bin_size) - 1
            end = int(gaussian_end_pixel / bin_size)
            sub_array_pixels = pixel_array[start:end]
            sub_array_intensities = master_array[start:end]

            mean = 0
            length = len(sub_array_pixels)
            for index in range(length):
                mean = mean + sub_array_pixels[index]/length
            my_sum = 0
            for index in range(length):
                my_sum = my_sum + (sub_array_pixels[index] - mean)**2 / length
            stddev = np.sqrt(my_sum)

            # plot the data
            # plt.plot(sub_array_pixels,sub_array_intensities,'b+:',label='data')

            # fit the data
            popt, pcov = curve_fit(gaus2, sub_array_pixels,sub_array_intensities, p0=[300, mean, stddev, 1])
            # print ("popt ", popt)

            # guassian statistic are contained in the popt array
            a = popt[0]
            x0 = popt[1]
            sigma = popt[2]
            c = popt[3]

            # approx solution to FWHM
            fwhm = 2.3548 * math.fabs(sigma)

            #increase resolution to plot the fit
            high_resolution_array = np.linspace(gaussian_start_pixel, gaussian_end_pixel, 100)

            #plot the higher resolution data
            # plt.plot(high_resolution_array, gaus2(high_resolution_array, *popt))

            #floor the gaussian fit
            popt[3] = 0

            #plot the higher resolution data as a floored function
            # plt.plot(high_resolution_array, gaus2(high_resolution_array, *popt))
            # plt.show()

            # move the guassian statistic to common variable names
            peak_value = gaus2(x0, *popt)
            peak_index = x0

    except:
        print ("no fit")

    results = {}
    results['fwhm'] = fwhm
    results['peak_value'] = peak_value
    results['peak_index'] = peak_index
    results['a'] = a

    return  results

def parseLine(line):
    return re.findall(r'\S+', line)

# finds the center of intensity in a heat map
def find_center_of_intensity_function(x, y, z, centerX, centerY, inner_percentage, max_size):
    x_sum = 0
    x_sub_sum = 0 # represents the mass (M) in computing the center of mass  1/M *Sum m*r
    y_sum = 0
    y_sub_sum = 0

    for index in range(len(z)):
        if ((x[index]-centerX)**2+(y[index]-centerY)**2)/max_size**2 < inner_percentage:
            x_sum = x_sum + x[index]*z[index]
            x_sub_sum = x_sub_sum +z[index]
            y_sum = y_sum + y[index]*z[index]
            y_sub_sum = y_sub_sum +z[index]
    x_dist = x_sum / x_sub_sum
    y_dist = y_sum / y_sub_sum

    print ("center of intensity pixel values")
    print ('x ' + str(x_dist))
    print ('y ' + str(y_dist))

main()
