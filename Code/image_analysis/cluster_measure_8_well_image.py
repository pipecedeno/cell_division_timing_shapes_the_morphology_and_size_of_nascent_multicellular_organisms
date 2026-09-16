
'''
Date: 27apr2023
This program was created to analyze and create the masks of the clusters for the images of a well of a 24-well plate. The images
are only bright field, using the 4x lense and the images are a 5x7 grid of 4x images.

Variables that the user may need to modify:
-min_circularity: variable used in Analize particles, it is the minimum circularity that the object need to have to be loaded 
	in the ROI manager
-max_circularity: variable used in Analize particles, it is the maximum circularity that the object need to have to be loaded 
	in the ROI manager
-min_size: variable used in Analize particles, it is the minimum size that the object need to have to be loaded 
	in the ROI manager, used to not process debri
-rolling_ball_radius: Variable used in the substract background
-blurring_radius: Variable used by the pseudo flat field correction of the Biovoxxel plugin

'''

from ij import IJ
from ij.plugin.frame import RoiManager
from ij.gui import PolygonRoi
from ij.gui import Roi
from java.awt import FileDialog
from ij.io import DirectoryChooser
from ij.process import ImageConverter
from ij import Prefs


import time
import os
from ij.gui import Overlay
from java.awt import BasicStroke


#This function saves the masks that were processed by the program on top of the bright field image, so it needs
#to open the image again and at the end it closes the image.
#Note: it uses the variable "roi" as a global variable
def save_masks_in_original(images_folder, file_name):
	imp = IJ.openImage(file_name)

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
	IJ.saveAs(imp, "Jpeg", images_folder + file_name.split("/")[-1].rstrip(".tif") + '_masks_original.jpg')
	imp.setOverlay(None)  # Remove the overlay before processing the next image
	imp.close()


#This function saves the masks of the ROI in the processed image, so in this case it is just to check more easily
#if the process is working as expected or not.
def save_masks_in_processed(imp, images_folder, file_name):
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
	IJ.saveAs(imp, "Jpeg", images_folder + file_name.split("/")[-1].rstrip(".tif") + '_masks_processed.jpg')
	imp.setOverlay(None)  # Remove the overlay before processing the next image


#Function to analyze each of the images
def analyze_image(file):
	
	imp = IJ.openImage(file_name)

	#converting to 8-bit
	ImageConverter.setDoScaling(True)
	IJ.run(imp, "8-bit", "")

	#Test gaussian blur
	#IJ.run(imp, "Gaussian Blur...", "sigma=0.5 scaled")

	#Pseudo flat field correction from the Biovoxxel plugin
	IJ.run(imp, "Pseudo flat field correction", "blurring="+str(blurring_radius)+" hide")

	#background substraction
	IJ.run(imp, "Subtract Background...", "rolling="+str(rolling_ball_radius)+" light")

	IJ.run(imp, "Find Edges", "")

	#Thresholding
	IJ.setAutoThreshold(imp, "Yen dark no-reset")
	Prefs.blackBackground = True
	IJ.run(imp, "Convert to Mask", "")

	#fill holes of the thresholding
	IJ.run(imp, "Fill Holes", "")

	#analyze particle using the global variables of size and circularity, it also excludes ROIs in the border of the image
	IJ.run(imp, "Analyze Particles...", "size="+str(min_size)+"-Infinity circularity="+str(min_circularity)+"-"+str(max_circularity)+" exclude add")

	rm.runCommand(imp, "Measure")

	#saving results
	IJ.saveAs("Results", file_name.replace('.tif', '_measurements_rois_border.csv'))

	#Creating image of the processed file with the masks used
	save_masks_in_processed(imp, images_folder, file_name)

	imp.close() #close image

	#save original bright field image with the ROIs
	save_masks_in_original(images_folder, file_name)

	IJ.run("Clear Results") #clear results window
	rm.runCommand("Reset") #reset ROI manager


# User input variables
pixels_stroke=4 #size of the stroke in pixels
min_circularity=0.2
max_circularity=1
min_size=50 #micron^2
rolling_ball_radius=50
blurring_radius=50


#Folder selection and search of the .tif files in the folder
dc = DirectoryChooser("Select Folder")
folder_path = dc.getDirectory()
print(folder_path)
list_files = []
for file in os.listdir(folder_path):
	if(".tif" in file):
		list_files.append(folder_path + file)

total_files = len(list_files)
num_files_processed = 0


#creating instances of the objects needed to execute the program, which are the ROI manager and the results table
#this variables are used as global variables so they are accesed in the functions
RM = RoiManager()
rm = RM.getRoiManager()

#this variables are just to output the progress of the program
total_files=len(list_files)
num_files_processed=0

IJ.run("Set Measurements...", "area mean standard min centroid perimeter fit shape integrated median area_fraction stack redirect=None decimal=5");

#Creating and writing the header of the log file
log_file_path = os.path.join(folder_path, "log_cluster_measure_8_well_plate.txt")
log_file = open(log_file_path, "w")
log_file.write("min_circularity: " + str(min_circularity) + "\n")
log_file.write("max_circularity: " + str(max_circularity) + "\n")
log_file.write("min_size: " + str(min_size) + "\n")
log_file.write("rolling_ball_radius: " + str(rolling_ball_radius) + "\n")
log_file.write("directory: " + folder_path + "\n\n")
log_file.write("Files analyzed:\n")

#creating folder with the images that will have the masks
images_folder = folder_path + "images_masks/"
if (not os.path.exists(images_folder)):
	os.makedirs(images_folder)

#iterate though each tif file in the directory
for file_name in list_files:
	
	#execute the main function, and if a False is returned then the file will be counted as if it had an error
	analyze_image(file_name)


	#This is just to print the progress of the program
	num_files_processed+=1
	print(str(num_files_processed)+'/'+str(total_files))

	log_file.write(file_name.split("/")[-1]+"\n")

#closing results window and roi manager
IJ.run("Close")
rm.close()


log_file.write("\nCode finished running without errors\n")
log_file.close()  # Close the log file after processing all files

#print the final message of the program just to tell if all files were analyzed correctly or not
print('All files were analyzed correctly')




