#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 10:15:31 2022

@author: tschakel
"""
from wad_qc.modulelibs import wadwrapper_lib
from MR_GeomAcc_util import *
import numpy as np
from skimage.color import rgb2gray
import pydicom
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Circle,Arc

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
     - Read all the images (done)
     - Stack slices properly (done)
     - Find isolines (done)
     - Create spheres from geometric center
         - Diameters: 20, 34, 45 cm (done)
     - Check if spheres are within isolines
         - 20 cm --> 1 mm isoline (done)
         - 34 cm --> 2 mm isoline (done)
         - 45 cm --> 5 mm isoline (done)
     - Output values:
         - 1 mm passed (done)
         - 2 mm passed (done)
         - 5 mm passed (done)
     - Extra statistics: 
         - total area of 1/2/3/5 mm accuracy? (done)
         - largest sphere possible within 1/2/5 mm isoline?
     - Output images:
         - Plot with 7 the slices with the sphere outlines (done)
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
    for fname in data_series[0]:
        files.append(pydicom.dcmread(fname))
    
    # sort on ImagePositionPatient[2] (z-position)
    sortedfiles = sorted(files, key=lambda f: f.ImagePositionPatient[2])

    # set some constants:
    radiusDSV = [100,170,225] #radii of the spheres in mm
    colorDSV = [[0,1,0],[0,1,1],[1,0,0]]
    
    # output 20cm -> 1mm, 34cm -> 2mm, 45cm -> 5mm
    spheres_inside = [True, True, True]
    # over all 7 slices 1mm ... 5mm
    total_area = [0.0, 0.0, 0.0, 0.0]
    # output figure
    fig, axs = plt.subplots(2,4)
    idx_axs = 0
    
    for file in sortedfiles:
        image_data = file.pixel_array
        
        set_bright_green(image_data)
        
        # red = 5 mm
        # yellow = 3 mm
        # teal = 2 mm
        # green = 1 mm
        # the cutoff values for the different colors are determined from the boxes in the bottom legend
        
        # find the isolines and fill them
        b_green, b_teal, b_yellow, b_red = get_rgb_lines_slice(image_data)
        
        # in b_green/teal/yellow/red are the points for all the boundaries and the corresponding masks
        # we can used these masks to create and overall mask maybe with something like a volume 
        # largest volume > 1's check if smaller volume is inside larger volume -> and so on
        m_green, m_teal, m_yellow, m_red = create_masks(b_green, b_teal, b_yellow, b_red)
        
        
        # add area of masks to total
        area_pixel = file.PixelSpacing[0]*1e-3*file.PixelSpacing[1]*1e-3
        total_area[0] += m_green.num_voxels * area_pixel
        total_area[1] += m_teal.num_voxels * area_pixel
        total_area[2] += m_yellow.num_voxels * area_pixel
        total_area[3] += m_red.num_voxels * area_pixel
        
        # the 1mm, 2mm, 5mm masks for which we check if certain spheres fit inside
        mask_to_check = [m_green, m_teal, m_red]
          
        # find pixel of the geometric center
        x0 = np.round( (0 - file.ImagePositionPatient[0] + (file.PixelSpacing[0] / 2)) / file.PixelSpacing[0])
        y0 = np.round( (0 - file.ImagePositionPatient[1] + (file.PixelSpacing[1] / 2)) / file.PixelSpacing[1])
        
        # find pixel where the table starts
        tablePix = np.round((float(params['table_cutoff']) - file.ImagePositionPatient[1] + (file.PixelSpacing[1] / 2)) / file.PixelSpacing[1])
        # cutoff for the top part of the phantom (should only matter for slice 3,4,5)
        topPix = np.round((-float(params['phantom_top_cutoff']) - file.ImagePositionPatient[1] + (file.PixelSpacing[1] / 2)) / file.PixelSpacing[1])
        
        # the 7 different slices are from 7 z-locations
        # they are not equidistant
        # the radius of the spheres at the intersections with the phantom slices are given by:
        # radius(z) = sqrt(R^2 - z^2), with R the radius of the current sphere
        # Note, the smaller spheres will not intersect with all slices
        curZ = file.ImagePositionPatient[2]
        
        #add to output plot
        axs[divmod(idx_axs,4)].imshow(image_data)
        axs[divmod(idx_axs,4)].set_title("slice = " + str(idx_axs+1))
        axs[divmod(idx_axs,4)].axis('off')
        for n in range(len(radiusDSV)):
            if radiusDSV[n]**2 - curZ**2 > 0: #check if current slice intersects with sphere
               rdsv = np.sqrt(radiusDSV[n]**2 - curZ**2) / file.PixelSpacing[0]
               
               #calculation
               if spheres_inside[n]: # once the value is set to false we do not need to check the other slices.
                    # the larger spheres can extent into the table, where there is no phantom
                    # avoid this region for plotting/analysis
                    # this means cutting of a part of the circle (or drawing just a part of it)
                    mask_circle = np.zeros(image_data.shape[0:2])
                    for y in range(int(topPix),int(tablePix)): #we start from cutoff of phantom and go until we find the table
                        for x in range(mask_circle.shape[1]):
                            if np.sqrt((x0-x)**2 + (y0-y)**2) < rdsv :
                                mask_circle[y,x] = 1
                    
                    spheres_inside[n] = mask_to_check[n].inside(mask_circle)
                    
               #plotting
               if (y0 - topPix) < rdsv:
                   # part 1
                   arcstart = 90 + np.arccos((tablePix-y0)/rdsv) * 180/np.pi
                   arcend = 90 + np.arccos((topPix-y0)/rdsv) * 180/np.pi
                   axs[divmod(idx_axs,4)].add_patch(Arc([x0,y0],2*rdsv,2*rdsv,angle=0,theta1=arcstart,theta2=arcend,linestyle=(3,(6,2)),ec=colorDSV[n]))
                   
                   # part 2
                   arcend = 90 - np.arccos((tablePix-y0)/rdsv) * 180/np.pi
                   arcstart = 90 - np.arccos((topPix-y0)/rdsv) * 180/np.pi
                   axs[divmod(idx_axs,4)].add_patch(Arc([x0,y0],2*rdsv,2*rdsv,angle=0,theta1=arcstart,theta2=arcend,linestyle=(3,(6,2)),ec=colorDSV[n]))
               elif (tablePix - y0) < rdsv:
                   arcstart = 90 + np.arccos((tablePix-y0)/rdsv) * 180/np.pi
                   arcend = 90 - np.arccos((tablePix-y0)/rdsv) * 180/np.pi
                   axs[divmod(idx_axs,4)].add_patch(Arc([x0,y0],2*rdsv,2*rdsv,angle=0,theta1=arcstart,theta2=arcend,linestyle=(3,(6,2)),ec=colorDSV[n]))
                    
                    # TO DO: connect the start and end of the Arc with a line     
               else:
                   axs[divmod(idx_axs,4)].add_patch(Circle([x0,y0],rdsv,fc='none',linestyle=(3,(6,2)),ec=colorDSV[n])) 
    
        idx_axs += 1 

    # to remove the empty subplot
    fig.delaxes(axs[divmod(idx_axs, 4)])
    filename = 'Geom_acc.png'
    fig.savefig(filename,dpi=300)
    
    # after looping over files add results          
    results.addBool("20cm inside 1mm iso", bool(spheres_inside[0]))
    results.addBool("34cm inside 2mm iso", bool(spheres_inside[1]))
    results.addBool("45cm inside 5mm iso", bool(spheres_inside[2]))
    results.addFloat("Area 1mm iso over 7 slices (m^2)", total_area[0])
    results.addFloat("Area 2mm iso over 7 slices (m^2)", total_area[1])
    results.addFloat("Area 3mm iso over 7 slices (m^2)", total_area[2])
    results.addFloat("Area 5mm iso over 7 slices (m^2)", total_area[3])
    results.addObject("Isolines for 7 slices", filename)
