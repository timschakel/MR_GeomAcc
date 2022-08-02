#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 10:15:31 2022

@author: tschakel
"""
from wad_qc.modulelibs import wadwrapper_lib
import numpy as np
from skimage.color import rgb2gray
import pydicom
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Circle,Arc

from MR_GeomAcc_util import arc_patch

### Helper functions
def getValue(ds, label):
    """Return the value of a pydicom DataElement in Dataset identified by label.

    ds: pydicom Dataset
    label: dicom identifier, in either pydicom Tag object, string or tuple form.
    """
    if isinstance(label, str):
        try:
            # Assume form "0x0008,0x1030"
            tag = pydicom.tag.Tag(label.split(','))
        except ValueError:
            try:
                # Assume form "SeriesDescription"
                tag = ds.data_element(label).tag
            except (AttributeError, KeyError):
                # `label` string doesn't represent an element of the DataSet
                return None
    else:
        # Assume label is of form (0x0008,0x1030) or is a pydicom Tag object.
        tag = pydicom.tag.Tag(label)

    try:
        return str(ds[tag].value)
    except KeyError:
        # Tag doesn't exist in the DataSet
        return None

def isFiltered(ds, filters):
    """Return True if the Dataset `ds` complies to the `filters`,
    otherwise return False.
    """
    for tag, value in filters.items():
        if not str(getValue(ds, tag)) == str(value):
            # Convert both values to string before comparison. Reason is that
            # pydicom can return 'str', 'int' or 'dicom.valuerep' types of data.
            # Similarly, the user (or XML) supplied value can be of any type.
            return False
    return True

def applyFilters(series_filelist, filters):
    """Apply `filters` to the `series_filelist` and return the filtered list.

    First, convert `filters` from an ElementTree Element to a dictionary
    Next, create a new list in the same shape as `series_filelist`, but only
    include filenames for which isFiltered returns True.
    Only include sublists (i.e., series) which are non empty.
    """
    # Turn ElementTree element attributes and text into filters
    #filter_dict = {element.attrib["name"]: element.text for element in filters}
    filter_dict = filters

    filtered_series_filelist = []
    # For each series in the series_filelist (or, study):
    for instance_filelist in series_filelist:
        # Filter filenames within each series
        filtered_instance_filelist = [fn for fn in instance_filelist
                                      if isFiltered(
                pydicom.read_file(fn, stop_before_pixels=True), filter_dict)]
        # Only add the series which are not empty
        if filtered_instance_filelist:
            filtered_series_filelist.append(filtered_instance_filelist)

    return filtered_series_filelist

def acqdatetime(data, results, action):
    """
    Get the date and time of acquisition
    """
    params = action["params"]
    datetime_series = data.getSeriesByDescription(params["datetime_series_description"])
    dt = wadwrapper_lib.acqdatetime_series(datetime_series[0][0])
    results.addDateTime('AcquisitionDateTime', dt)
    
def GeomAcc_rgb(data, results, action):
    """
    Process Geometric Accuracy measurement
     - Read all the images
     - Stack slices properly
     - Find isolines
     - Create spheres from geometric center
         - Diameters: 20, 34, 45 cm
     - Check if spheres are within isolines
         - 20 cm --> 1 mm isoline
         - 34 cm --> 2 mm isoline
         - 45 cm --> 5 mm isoline
     - Output values:
         - 1 mm passed
         - 2 mm passed
         - 5 mm passed
     - Extra statistics: 
         - total area of 1/2/3/5 mm accuracy?
         - largest sphere possible within 1/2/5 mm isoline?
     - Output images:
         - Plot with 7 the slices with the sphere outlines
         - Extra ?
    """
    params = action["params"]
    filters = action["filters"]
    
    series_filter = {"SeriesDescription":filters.get(item)for item in ["series_description"]}
    data_series = applyFilters(data.series_filelist, series_filter)
    
    # The usual reader function will fail on this data. Probably missing RescaleSlope/Intercept
    # dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(data_series[0],headers_only=False)
    # check if dicom contains RGB:
    # data1 = pydicom.dcmread(data_series[0][3])
    # data1.PhotometricInterpretation (should be 'RGB')
    
    files = []
    for fname in data.series_filelist[0]:
        files.append(pydicom.dcmread(fname))
    
    # sort on ImagePositionPatient[2] (z-position)
    sortedfiles = sorted(files, key=lambda f: f.ImagePositionPatient[2])

    # set some constants:
    radiusDSV = [100,170,225] #radii of the spheres
    colorDSV = ['red','teal','green']
    
    breakpoint()
    for file in sortedfiles:
        image_data = file.pixel_array
        image_data_gray = rgb2gray(image_data)
        
        # red = 5 mm
        # yellow = 3 mm
        # teal = 2 mm
        # green = 1 mm
        # the cutoff values for the different colors are determined from the boxes in the bottom legend
        # could be done automatically, position of boxes is always the same (but pixelvalues too?)
        lineRed = image_data_gray == 0.2125
        lineYellow = image_data_gray == 0.9279000000000001
        lineTeal = image_data_gray == 0.7875000000000001
        lineGreen = image_data_gray == 0.1683294117647059
        
        # TO DO: lines to mask
        
        # find pixel of the geometric center
        x0 = np.round( (0 - file.ImagePositionPatient[0] + (file.PixelSpacing[0] / 2)) / file.PixelSpacing[0])
        y0 = np.round( (0 - file.ImagePositionPatient[1] + (file.PixelSpacing[1] / 2)) / file.PixelSpacing[1])
        
        # find pixel where the table starts
        tablePix = np.round((float(params['table_cutoff']) - file.ImagePositionPatient[1] + (file.PixelSpacing[1] / 2)) / file.PixelSpacing[1])
        
        # the 7 different slices are from 7 z-locations
        # they are not equidistant
        # the radius of the spheres at the intersections with the phantom slices are given by:
        # radius(z) = sqrt(R^2 - z^2), with R the radius of the current sphere
        # Note, the smaller spheres will not intersect with all slices
        curZ = file.ImagePositionPatient[2]
        
        for n in range(len(radiusDSV)):
            if radiusDSV[n]**2 - curZ**2 > 0: #check if current slice intersects with sphere
                rdsv = np.sqrt(radiusDSV[n]**2 - curZ**2) / file.PixelSpacing[0]
                
                #plot for testing
                fig,axs = plt.subplots()
                axs.imshow(image_data_gray,cmap='gray')
                
                # the larger spheres can extent into the table, where there is no phantom
                # avoid this region for plotting/analysis
                # this means cutting of a part of the circle (or drawing just a part of it)
                if (tablePix - y0) < rdsv:
                    arcstart = 90 + np.arccos((tablePix-y0)/rdsv) * 180/np.pi
                    arcend = 90 - np.arccos((tablePix-y0)/rdsv) * 180/np.pi
                    axs.add_patch(Arc([x0,y0],2*rdsv,2*rdsv,angle=0,theta1=arcstart,theta2=arcend,lw=2,ec=colorDSV[n]))
                    
                    # TO DO: connect the start and end of the Arc with a line     
                else:
                    axs.add_patch(Circle([x0,y0],rdsv,fc='none',lw=2,ec=colorDSV[n])) 
                
                axs.axis('off')
                
                # TO DO: sphere to mask
                
                # TO DO: Check if spheres are within isolines (use masks?)