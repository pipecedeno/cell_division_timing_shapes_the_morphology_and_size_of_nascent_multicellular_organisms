'''
Date: 20May2025
Version: 4

Modified version of measure_all_rois_by_group_v3.py with updated image processing steps.
Changes include:
- Replaced background subtraction with Biovoxxel pseudo flat field correction
- Using Yen thresholding
- Kept the dilate, close, and fill holes operations

Inputs:
In the window that appears when running the program you need to select the tif file that you want to analyze, this file 
should have in the same directory a txt file of the masks that needs to have the sufix "_cp_outlines.txt", if there's no
file with that sufix then the program won't be executed.

Output:
The output will be written in the same directory of the tif file and will have the suffix '_measurements_by_group.csv'.
From this file you need to delete or filter the rows that have an NA in the column 'group', as those rows belong to the ROIs
of the groups or to ROIs that were not assign to any group so they must be error of the detection step performed by cellpose.
Also it was added a jpg image of every file processed so that the user can check if the thresholding worked as expected.

WARNING:
When using this program always check that the steps used to threshold the image are good for the images that you have collected,
you can check if the process was run as expected by checking the jpg images in the directory
'''

from ij import IJ
from ij.plugin.frame import RoiManager
from ij.gui import PolygonRoi
from ij.gui import Roi
from java.awt import FileDialog
from ij.process import ImageConverter
from ij.measure import ResultsTable

import os  # library needed to check the files in a directory

'''
Function used to get the x and y coordinates of the center of the rectangle of a ROI

The ROIs are saved in imagej as a rectangle, so to get the center of the ROI you just need to get the center of the rectangle.
Formula used:
center coordinates = (x+w/2, y+h/2)
'''
def get_center_from_floatpolygon(floatpolygon):
    parts = floatpolygon.toString().split(' ')[1:]
    dict = {}
    for elem in parts:
        values = elem.split('=')  # separating the expressions 'x=n' to 2 separate values
        variable = values[0][0]
        number = values[1].replace(',', '').replace(']', '')  # the replace functions are to delete some missing 
        # characters that aren't deleted by the split
        dict[variable] = int(number)
    temp_center = [int(dict['x'] + dict['w'] / 2), int(dict['y'] + dict['h'] / 2)]
    # converting to int to avoid getting float values
    return(temp_center)


'''
This is the main function of the program, the function expects a string as input and will analyze 
that file if a txt file with the masks is found
the function return False if there's no txt file with the masks, and True if the file exists
'''
def analyze_image(file):
    
    # Parameters used:
    rolling_ball_radius = 10  # kept for backward compatibility but not used anymore
    min_area_group = 40  # um^2
    num_dilates = 5
    num_close = 3
    flag_substract_background = False  # set to False as we're using pseudo flat field correction instead
    
    # checking if the file exists, if it doesn't exist then it won't execute the rest of the function
    try:
        f = open(file.rstrip(".tif") + "_cp_outlines.txt")
        f.close()
    except IOError:
        return(False)

    # opening the image
    imp = IJ.openImage(file)
    imp.show()

    ##Processing the image
    # converting to 8-bit
    ImageConverter.setDoScaling(True)
    # IJ.run(imp, "8-bit", "")

    # IJ.run(imp, "Gaussian Blur...", "sigma=0.10 scaled")

    # Apply Biovoxxel pseudo flat field correction 
    IJ.run(imp, "Pseudo flat field correction", "blurring=10 hide")
    
    # find edges
    IJ.run(imp, "Find Edges", "")

    # thresholding with Yen
    IJ.setAutoThreshold(imp, "Yen dark")
    
    # applying threshold
    IJ.run(imp, "Convert to Mask", "")

    # dilate, close, fill holes operations
    IJ.run(imp, "EDM Binary Operations", "iterations=" + str(num_dilates) + " operation=dilate")
    IJ.run(imp, "EDM Binary Operations", "iterations=" + str(num_close) + " operation=close")
    IJ.run(imp, "Fill Holes", "")

    # analyze particle, the minimum size of the particle should also be tested
    IJ.run(imp, "Analyze Particles...", "size=" + str(min_area_group) + "-Infinity add")

    # count the number of groups found
    # this value is important to later iterate through the values of the masks of the groups
    n_groups = rm.getCount()

    # opening the specific txt file of the image and loading the rois as shown in imagej_roi_converter_load_all_files.py
    textfile = open(file.rstrip(".tif") + "_cp_outlines.txt", "r")

    cont = 1
    for line in textfile:
        xy = map(int, line.rstrip().split(","))
        X = xy[::2]
        Y = xy[1::2]
        imp.setRoi(PolygonRoi(X, Y, Roi.POLYGON))

        roi = imp.getRoi()
        rm.addRoi(roi)

    textfile.close()

    rm.runCommand(imp, "Measure")

    # getting total amount of ROIs, so n_rois include the rois of the groups and the rois of cellpose
    n_rois = rm.getCount()

    # iterating through the rois of cellpose
    for index_roi in range(n_groups, n_rois):
        # getting information of the ROI
        cell_roi = rm.getRoi(index_roi)

        # getting the center of the ROI of the cell
        center = get_center_from_floatpolygon(cell_roi)
        
        in_a_group = False  # this variable is important as a Roi may not
        # be contained in a group
        
        # iterating through the groups found
        for index_group in range(0, n_groups):
            group_roi = rm.getRoi(index_group)
            # checking if the center of the cell is in the rectangle of the group ROI, if
            # it is then there's no need to continue iterating, and the index were it stopped+1
            # is the group were the cell belongs
            if(group_roi.contains(center[0], center[1])):
                in_a_group = True
                break
        
        # if the roi of the cell was found within a roi of a group
        # then that is going to be added to the group column
        if(in_a_group):
            rt.setValue('group', index_roi, index_group + 1)
            rt.show('Results')
        else:
            # if it isn't found in any group then an NA is added and this will work as a filter for debris
            # that is miss labelled by cellpose
            rt.setValue('group', index_roi, 'NA')
            rt.show('Results')

    # this last for will be used to set the group ROIs as NA so that we
    # can differentiate them from the ROIs of cellpose,
    for i in range(0, n_groups):
        rt.setValue('group', i, 'NA')
        rt.show('Results')

    # Saving the results table
    IJ.saveAs("Results", file.rstrip(".tif") + '_measurements_by_group.csv')

    # saving a jpg image of all the process so that afterwards I can check if it worked correctly
    flat = imp.flatten()
    IJ.saveAs(flat, "Jpeg", file.rstrip(".tif") + '_masks.jpg')
    flat.close()

    # to close every window and clean the ROI manager
    IJ.run("Clear Results", "")
    imp.changes = False
    imp.close()  # close image
    rm.reset()  # reset ROI manager

    return(True)

# creating instances of the objects needed to execute the program, which are the ROI manager and the results table
# this variables are used as global variables so they are accessed in the functions
RM = RoiManager()
rm = RM.getRoiManager()
rt = ResultsTable.getResultsTable()

# getting the information of the directory where all the tif files that should be analyzed are located
fd = FileDialog(IJ.getInstance(), "Open", FileDialog.LOAD)
fd.show()

file_dir = fd.getDirectory()

# creating a list with all the tif files in a directory
list_files = []
for file in os.listdir(file_dir):
    if(".tif" in file):
        list_files.append(file_dir + file)

# this variables are just to output the progress of the program
total_files = len(list_files)
num_files_processed = 0

# this will count the files that didn't had a txt file and save the name of the file so that later they can be 
# printed in the terminal
cont_files_with_error = 0
list_files_error = []

# iterate though each tif file in the directory
for file_name in list_files:
    
    # execute the main function, and if a False is returned then the file will be counted as if it had an error
    if(not analyze_image(file_name)):
        cont_files_with_error += 1
        list_files_error.append(file_name)

    # This is just to print the progress of the program
    num_files_processed += 1
    print(str(num_files_processed) + '/' + str(total_files))

# print the final message of the program just to tell if all files were analyzed correctly or not
if(cont_files_with_error == 0):
    print('All files were analyzed correctly')
else:
    print('Warning!!!\nNot all files had their txt file.')
    print('The following files had the error:')
    for elem in list_files_error:
        print(elem)