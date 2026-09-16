
'''

Date:September 21 2023


This program is a modification of imagej_roi_converter_ignore_border_cells_v2.py. 
The idea of this modification is to allow the user to load the segmentation information from txt files or zip files if the manual
correction was done in imagej.

Inputs:
	-tif folder: folder were the .tif files are saved. All files that end with .tif are going to be analyzed

	-txt/zip folder: folder were either the txt files from cellpose or the zip files from imagej are saved
		-The txt file needs to end in "_cp_outlines.txt" to be processed
		-The zip file needs to end in "_rois.zip" to be processed
	Note: if the folder has both a txt and a zip file with the correct sufix, the zip file is going to be used.

	-output folder: folder where the following files are going to be saved:
		-csv files with the measurements
		-images_folder which will have jpg images of the bright field of the images with the ROIs in it
		-if the rois were loaded from a txt file, then a zip file will be saved with the ROIs from fiji


Warning:
	-The input images for which the script was designed had blue and bright field, to save the images in the field that you want change the
	variable to the number that you want, if images with only one field are used then an error may appear.


The user needs to define the variables:
	-dist_border_input: distance of the border that the cell needs to be in to be counted
	-flag_microns: if True then the units of the distance variable is considered as microns, else if False the units of the distance
	is pixels
	-flag_save_images: if the flag is set to True, then a folder is going to be created which will contain the bright field image
	-pixels_stroke: size of the stroke when saving the images of the masks
	-channel_for_image: which channel is going to be used in the images of the masks


Modification history:
Changed the measure parameters
Added to save the ROIs .zip file only if ROIs were extracted from txt files
Added a window that will ask for the input folders

'''

#libraries needed for the dialogs and input windows
from ij.gui import GenericDialog
import os
from java.awt import FileDialog, Frame
from ij import IJ

#libraries needed for the rest of the program
from ij.plugin.frame import RoiManager
from ij.gui import PolygonRoi
from ij.gui import Roi
from ij.io import DirectoryChooser

import time
import os
from ij.gui import Overlay
from java.awt import BasicStroke

from ij.measure import ResultsTable
from java.awt import Color


'''
#####
Input dialog functions
#####
'''

#Main function that creates the window for the different inputs
def get_user_inputs():
	# Create a new GenericDialog
	gd = GenericDialog("Select Folders and Input Values")

	gd.addMessage("Select tif files folder:")
	gd.addStringField("Tif folder:", "")
	gd.addButton("Select", lambda event: select_folder(gd, 1))  # Call select_folder() without arguments
	gd.addMessage("Select txt/zip files folder:")
	gd.addStringField("Txt/zip folder:", "")
	gd.addButton("Select", lambda event: select_folder(gd, 3))  # Call select_folder() without arguments
	gd.addMessage("Output directory:")
	gd.addStringField("Output dir:", "")
	gd.addButton("Select", lambda event: select_folder(gd, 5))  # Call select_folder() without arguments

	# Show the dialog
	gd.showDialog()

	# Check if the user clicked the continue button
	if gd.wasOKed():
		# save the folder paths
		tif_folder = gd.getNextString()
		txt_folder = gd.getNextString()
		output_dir = gd.getNextString()

		# Return the inputs as a dictionary
		return ({
			'tif_folder': tif_folder,
			'txt_folder': txt_folder,
			'output_dir': output_dir,
		})

	else:
		# Return None if the dialog was canceled
		return (None)


#function to select each individual folder, this function is called if the "Select" button is pressed
def select_folder(gd, index):

	dc = DirectoryChooser("Select Folder")
	folder_path = dc.getDirectory()

	# Update the selected folder path in the GenericDialog
	if folder_path:
		gd.getStringFields()[index//2].setText(os.path.abspath(folder_path))


#function to check that non of the inputs in the dictionary of results is empty
def check_inputs(user_inputs):
	# Iterate through all the keys in the user_inputs dictionary
	for key, value in user_inputs.items():
		# Check if the value is empty or None
		if not value:
			print("Input ",key," is missing or empty.")
			return False
	
	# All inputs are filled
	return(True)


'''
####
File processing functions
####
'''

def is_roi_close_to_border(X, Y, dist_border, width, height):
	for x, y in zip(X, Y):
		if x < dist_border or x > width - dist_border or y < dist_border or y > height - dist_border:
			return True
	return False

def save_masks_in_jpg(imp, rm, images_folder, file_name):
	# Change to the selected channel
	imp.setC(channel_for_image)
	
	# Create an Overlay object to store the masks
	overlay = Overlay()

	# Iterate through the ROIs and add them to the overlay
	for i in range(rm.getCount()):
		roi = rm.getRoi(i)
		roi.setStroke(BasicStroke(pixels_stroke))  # Set the stroke width to X pixels
		overlay.add(roi)

	# Set the overlay to the image and save it
	imp.setOverlay(overlay)
	imp.updateAndDraw()
	IJ.saveAs(imp, "Jpeg", os.path.join(images_folder, file_name.split("/")[-1].rstrip(".tif") + '_masks.jpg'))
	imp.setOverlay(None)  # Remove the overlay before processing the next image


def process_each_tif(tif_file, user_inputs, flag_microns, flag_save_images, images_folder, log_file):
	# Retrieve the input values from the dictionary
	tif_folder = user_inputs['tif_folder']
	txt_folder = user_inputs['txt_folder']
	output_dir = user_inputs['output_dir']

	# print(tif_file)
	imp = IJ.openImage(tif_file)



	#Code to change the channel used for calculation, important for pixel intensity
	# imp.setC(channel_for_image)

	width = imp.getWidth()
	height = imp.getHeight()
	calibration = imp.getCalibration()

	txt_file=os.path.join(txt_folder, tif_file.split('/')[-1].replace(".tif","_cp_outlines.txt"))
	zip_file=os.path.join(txt_folder, tif_file.split('/')[-1].replace(".tif","_rois.zip"))

	#checking if a zip file exists
	flag_zip_file_used=False #this flag is to know if later the zip file should be saved or not
	if(os.path.exists(zip_file)):
		#Open zip file to ROI manager
		rm.open(zip_file)
		flag_zip_file_used=True
	#if the zip file doesn't exist then use the txt file
	elif(os.path.exists(txt_file)):
		if flag_microns:
			dist_border = dist_border_input / calibration.pixelWidth

		textfile = open(txt_file, "r")
		for line in textfile:
			xy = map(int, line.rstrip().split(","))
			X = xy[::2]
			Y = xy[1::2]

			if is_roi_close_to_border(X, Y, dist_border, width, height):
				continue

			imp.setRoi(PolygonRoi(X, Y, Roi.POLYGON))
			roi = imp.getRoi()
			rm.addRoi(roi)

		textfile.close()

	rm.runCommand(imp, "Measure")

	if(flag_save_images):
		save_masks_in_jpg(imp, rm, images_folder, tif_file)

	output_file_name=os.path.join(output_dir, tif_file.split('/')[-1].replace('.tif', '_measurements_rois_border.csv'))
	IJ.saveAs("Results", output_file_name)
	IJ.run("Clear Results")

	if(flag_zip_file_used==False):
		zip_file_name=os.path.join(output_dir, tif_file.split('/')[-1].replace('.tif', '_rois.zip'))
		rm.save(zip_file_name);

	rm.runCommand("Reset")
	imp.close()



'''
####
Main part of the program
####
'''

# User input variables
dist_border_input = 0.5
flag_microns = True
flag_save_images = True
pixels_stroke=4 #size of the stroke in pixels
channel_for_image=2 


# Call the function to get user inputs
user_inputs = get_user_inputs()

# Check if user_inputs is not None (dialog was not canceled)
if (user_inputs is None) or (check_inputs(user_inputs)==False):
	print("Program execution stopped.")

else:

	# dc = DirectoryChooser("Select Folder")
	# folder_path = dc.getDirectory()
	# print(folder_path)

	# Retrieve the input values from the dictionary
	tif_folder = user_inputs['tif_folder']
	txt_folder = user_inputs['txt_folder']
	output_dir = user_inputs['output_dir']


	list_files = []
	for file in os.listdir(tif_folder):
		if(".tif" in file):
			list_files.append(os.path.join(tif_folder, file))

	# print(list_files)

	total_files = len(list_files)
	num_files_processed = 0

	RM = RoiManager()
	rm = RM.getRoiManager()

	# IJ.run("Set Measurements...", "area mean fit shape integrated redirect=None decimal=3")
	IJ.run("Set Measurements...", "area mean standard min centroid perimeter fit shape integrated median area_fraction stack redirect=None decimal=5");

	log_file_path = os.path.join(output_dir, "log_measure_all_rois_in_dir_ignore_border_cells_v2.1.txt")
	log_file = open(log_file_path, "w")
	log_file.write("dist_border_input: " + str(dist_border_input) + "\n")
	log_file.write("flag_microns: " + str(flag_microns) + "\n")
	log_file.write("tif_folder: " + tif_folder + "\n")
	log_file.write("txt/zip_folder: " + txt_folder + "\n")
	log_file.write("output_folder: " + output_dir + "\n\n")
	log_file.write("Files analyzed:\n")


	#create output directory if it doesn't exist
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	if(flag_save_images):
		images_folder = os.path.join(output_dir, "images_masks")
		if not os.path.exists(images_folder):
			os.makedirs(images_folder)

	for file_name in list_files:
		
		#process each file code
		# print(file_name)
		process_each_tif(file_name, user_inputs, flag_microns, flag_save_images, images_folder, log_file)

		num_files_processed += 1
		print(str(num_files_processed) + '/' + str(total_files))

	#if(not flag_save_images):
	IJ.run("Close")
	rm.close()

	log_file.write("\nCode finished running without errors\n")
	log_file.close()  # Close the log file after processing all files
	print('Process finished running without any error')

