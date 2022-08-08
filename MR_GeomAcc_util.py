#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 13:53:06 2022

@author: pstijnm2
"""
import numpy as np
from skimage.segmentation import flood_fill
import skimage.morphology

class boundary:
    def __init__(self, points, image,x0,y0,tablePix,topPix,curZ,pixelSpacing):
        self.points = points
        self.tablePix = tablePix
        self.topPix = topPix
        self.curZ = curZ
        self.pixelSpacing = pixelSpacing
        self.mask = self.fill(image)
        self.num_voxels = self.set_num_voxels()
        self.min_rad = self.set_min_rad(x0,y0)
        
        
    def __lt__(self,other): # for sorting
        return self.num_voxels < other.num_voxels
        
    def fill(self, image):
        mask = np.ones(image.shape[0:2]) # create 2d mask
        for point in self.points:
            mask[point[0], point[1]] = 2
            
        mask = flood_fill(mask, (0,0), 0, selem=np.array([0,1,0,1,0,1,0,1,0]).reshape([3,3]))
        
        for point in self.points:
            mask[point[0], point[1]] = 1
        
        return mask
    
    def set_num_voxels(self):
        return np.sum(self.mask)
    
    # check if other mask is inside
    def inside(self,other):
        tmp = np.multiply(self.mask, other.mask)
        return np.sum(tmp) == other.num_voxels
    
    # for each point, add distance to origin (x0,y0)
    def set_min_rad(self,x0,y0):
        min_rad = 999999
        #if proper vector/array is created first, we can skip for loop? (speed gain?)
        for point in self.points:
            if point[0] > self.topPix and point[0] < self.tablePix:
                rad2d = np.sqrt((point[0] - x0)**2 + (point[1] - y0)**2) * self.pixelSpacing
                rad3d = np.sqrt(self.curZ**2 + rad2d**2)
                if rad3d < min_rad:
                    min_rad = rad3d
        
        return min_rad
        
class Mask:
    def __init__(self, mask):
        self.mask = mask
        self.num_voxels = self.set_num_voxels()
        
    def set_num_voxels(self):
        return np.sum(self.mask)

    def inside(self, mask):
        tmp = np.multiply(self.mask, mask)
        return np.sum(tmp) == np.sum(mask)
            
def create_mask(boundaries):
    if len(boundaries) == 1:
        return Mask(boundaries[0].mask)
    
    boundaries.sort()
    boundaries.reverse()
    mask = np.zeros(boundaries[0].mask.shape)
    mask += boundaries[0].mask
    
    for i in range(1,len(boundaries)):
        count = 0
        for j in range(i):
            count += boundaries[j].inside(boundaries[i])
        
        if count%2: # inside odd number of boundaries
            mask -= boundaries[i].mask
        else:
            mask += boundaries[i].mask
    
    return Mask(mask)
'''
create the overall masks given de all the filled in boundaries
'''
def create_masks(b_green, b_teal, b_yellow, b_red):
    m_green = create_mask(b_green)
    m_teal = create_mask(b_teal)
    m_yellow = create_mask(b_yellow)
    m_red = create_mask(b_red)
    return m_green, m_teal, m_yellow, m_red

def get_boundaries_from_points(points, image,x0,y0,tablePix,topPix,curZ,pixelSpacing):
    tmp_image = np.zeros(image.shape[0:2], dtype=np.bool)
    boundaries = list[boundary]()
    
    for p in points:
        tmp_image[p] = True
        
    labeled = skimage.morphology.label(tmp_image)
    coords = { i: (labeled == i).nonzero() for i in range(1,labeled.max()+1) }
    for key in coords:
        boundaries.append(boundary(list(zip(coords[key][0],coords[key][1])), image, x0, y0,tablePix,topPix,curZ,pixelSpacing))
        
    return boundaries
'''
given an image find yellow red teal and green pixels
returns the boundaries for each color
'''
def get_rgb_lines_slice(image,x0,y0,tablePix,topPix,curZ,pixelSpacing):
    red = [255,0,0]
    yellow = [255,255,0]
    teal = [0,255,255]
    green = [0,255,0]
    
    # 5mm
    indices_red = np.where(np.all(image == red, axis=-1))
    indices_legend = np.where(indices_red[0] > image.shape[0]-100) #skip the legend at bottom
    line_red = list(zip(np.delete(indices_red[0],indices_legend),np.delete(indices_red[1],indices_legend)))
    boundaries_red = get_boundaries_from_points(line_red, image,x0,y0,tablePix,topPix,curZ,pixelSpacing)
    
    # 3mm
    indices_yellow = np.where(np.all(image == yellow, axis=-1))
    indices_legend = np.where(indices_yellow[0] > image.shape[0]-100)
    line_yellow = list(zip(np.delete(indices_yellow[0],indices_legend),np.delete(indices_yellow[1],indices_legend)))
    boundaries_yellow = get_boundaries_from_points(line_yellow, image,x0,y0,tablePix,topPix,curZ,pixelSpacing)
    
    # 2mm 
    indices_teal = np.where(np.all(image == teal, axis=-1))
    indices_legend = np.where(indices_teal[0] > image.shape[0]-100)
    line_teal = list(zip(np.delete(indices_teal[0],indices_legend),np.delete(indices_teal[1],indices_legend)))
    boundaries_teal = get_boundaries_from_points(line_teal, image,x0,y0,tablePix,topPix,curZ,pixelSpacing)
    
    # 1mm
    indices_green = np.where(np.all(image == green, axis=-1))
    indices_legend = np.where(indices_green[0] > image.shape[0]-100)
    line_green= list(zip(np.delete(indices_green[0],indices_legend),np.delete(indices_green[1],indices_legend)))
    boundaries_green = get_boundaries_from_points(line_green, image,x0,y0,tablePix,topPix,curZ,pixelSpacing)
    
    return boundaries_green, boundaries_teal, boundaries_yellow, boundaries_red
    
def set_bright_green(image):
    green = [0,60,0]
    bright_green = [0,255,0]
    
    indices_green = np.where(np.all(image == green, axis=-1))
    image[indices_green] = bright_green
                
def calc_min_rad_from_origin(b_green, b_teal, b_yellow, b_red):
    # This does not take the cutoff for top/table into account
    # Returns the smallest distance from the origin to the isolines
    # Will give counter intuitive results: largest possible spheres are smaller than the passed spheres....
    b_to_check = [b_green, b_teal, b_yellow, b_red]
    min_rad_from_ori = [999999.0, 999999.0, 999999.0, 999999.0]
    for n in range(len(b_to_check)):
        for b in range(len(b_to_check[n])):
            min_rad = b_to_check[n][b].min_rad
            if min_rad < min_rad_from_ori[n]:
                min_rad_from_ori[n] = min_rad
        
    return min_rad_from_ori