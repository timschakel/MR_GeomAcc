#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 13:53:06 2022

@author: pstijnm2
"""
import numpy as np
from skimage.segmentation import flood_fill
class boundary:
    def __init__(self, points, image):
        self.points = points
        self.mask = self.fill(image)
        self.num_voxels = self.set_num_voxels()
        
        
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
    
    

'''
given a list of points
start with the first point
find all the neighbours
once the list of connected points stops growing and all the points in the 
list are checked we return the boundary points
'''
def get_connected_points(points,xs,ys):
    boundary_points = list()
    current_point = points[0]
    boundary_points.append(current_point)
    idx = 0
    while (True):
        for y in ys:
            for x in xs:
                check_point = (current_point[0]-x, current_point[1]-y)
                if (check_point in points and check_point not in boundary_points):
                    boundary_points.append(check_point)
    
        idx += 1
        if (idx >= len(boundary_points)):
            break
        current_point = boundary_points[idx]
        
    return boundary_points

"given the points of a specific rgb color return the connected points in sets of boundaries"
def get_boundaries_from_points(points, image):
    xs = [-1,0,1]
    ys = [-1,0,1]
    boundaries = list()
    
    while (len(points) > 0):
        connected = get_connected_points(points,xs,ys)
        boundaries.append(boundary(connected, image))
        points = list(set(points) ^ set(connected)) # returns list of points not in connected
        
    return boundaries

                
'''
if a -> [r,g,b] == b -> [r,g,b] return 1
'''
def compare_rgb(a,b):
    return a[0] == b[0] and a[1] == b[1] and a[2] == b[2]

'''
go over the image and found all points with a specific rgb value
'''
def get_points_rgb_color(image, rgb):
    points = list()
    
    for y in range(image.shape[1]): 
        for x in range(image.shape[0]-100):# skip part with the legend
            if compare_rgb(image[x,y,:], rgb):
                points.append((x,y))
    
    return points

'''
given an image find yellow red teal and green pixels
returns the boundaries for each color
'''
def get_rgb_lines_slice(image):
    red = [255,0,0]
    yellow = [255,255,0]
    teal = [0,255,255]
    green = [0,255,0]
    
    # 5mm
    line_red = get_points_rgb_color(image, red)
    boundaries_red = get_boundaries_from_points(line_red, image)
    # 3mm
    line_yellow = get_points_rgb_color(image, yellow)
    boundaries_yellow = get_boundaries_from_points(line_yellow, image)
    # 2mm 
    line_teal = get_points_rgb_color(image, teal)
    boundaries_teal = get_boundaries_from_points(line_teal, image)
    # 1mm
    line_green = get_points_rgb_color(image, green)
    boundaries_green = get_boundaries_from_points(line_green, image)
    
    return boundaries_green, boundaries_teal, boundaries_yellow, boundaries_red
    

def set_bright_green(image):
    green = [0,60,0]
    bright_green = [0,255,0]
    
    for y in range(image.shape[1]):
        for x in range(image.shape[0]):
            if compare_rgb(image[x,y], green):
                image[x,y] = bright_green