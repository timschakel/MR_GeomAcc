#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 13:53:06 2022

@author: pstijnm2
"""
import numpy as np
from skimage.segmentation import flood_fill
class boundary:
    def __init__(self, points, image,x0,y0):
        self.points = points
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
            rad = np.sqrt((point[0] - x0)**2 + (point[1] - y0)**2)
            if rad < min_rad:
                min_rad = rad
        
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
def get_boundaries_from_points(points, image,x0,y0):
    xs = [-1,0,1]
    ys = [-1,0,1]
    boundaries = list()
    
    while (len(points) > 0):
        connected = get_connected_points(points,xs,ys)
        boundaries.append(boundary(connected, image,x0,y0))
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
def get_rgb_lines_slice(image,x0,y0):
    red = [255,0,0]
    yellow = [255,255,0]
    teal = [0,255,255]
    green = [0,255,0]
    
    # 5mm
    #line_red = get_points_rgb_color(image, red)
    indices_red = np.where(np.all(image == red, axis=-1))
    indices_legend = np.where(indices_red[0] > image.shape[0]-100) #skip the legend at bottom
    line_red = list(zip(np.delete(indices_red[0],indices_legend),np.delete(indices_red[1],indices_legend)))
    boundaries_red = get_boundaries_from_points(line_red, image,x0,y0)

    # 3mm
    #line_yellow = get_points_rgb_color(image, yellow)
    indices_yellow = np.where(np.all(image == yellow, axis=-1))
    indices_legend = np.where(indices_yellow[0] > image.shape[0]-100)
    line_yellow = list(zip(np.delete(indices_yellow[0],indices_legend),np.delete(indices_yellow[1],indices_legend)))
    boundaries_yellow = get_boundaries_from_points(line_yellow, image,x0,y0)
    
    # 2mm 
    #line_teal = get_points_rgb_color(image, teal)
    indices_teal = np.where(np.all(image == teal, axis=-1))
    indices_legend = np.where(indices_teal[0] > image.shape[0]-100)
    line_teal = list(zip(np.delete(indices_teal[0],indices_legend),np.delete(indices_teal[1],indices_legend)))
    boundaries_teal = get_boundaries_from_points(line_teal, image,x0,y0)
    
    # 1mm
    #line_green = get_points_rgb_color(image, green)
    indices_green = np.where(np.all(image == green, axis=-1))
    indices_legend = np.where(indices_green[0] > image.shape[0]-100)
    line_green= list(zip(np.delete(indices_green[0],indices_legend),np.delete(indices_green[1],indices_legend)))
    boundaries_green = get_boundaries_from_points(line_green, image,x0,y0)
    
    return boundaries_green, boundaries_teal, boundaries_yellow, boundaries_red
    
def set_bright_green(image):
    green = [0,60,0]
    bright_green = [0,255,0]
    
    for y in range(image.shape[1]):
        for x in range(image.shape[0]):
            if compare_rgb(image[x,y], green):
                image[x,y] = bright_green
                
def calc_min_rad_from_origin(b_green, b_teal, b_yellow, b_red,curZ,pixSpacing):
    # This does not take the cutoff for top/table into account
    # Returns the smallest distance from the origin to the isolines
    # Will give counter intuitive results: largest possible spheres are smaller than the passed spheres....
    b_to_check = [b_green, b_teal, b_yellow, b_red]
    min_rad_from_ori = [999999.0, 999999.0, 999999.0, 999999.0]
    for n in range(len(b_to_check)):
        for b in range(len(b_to_check[n])):
            min_rad = b_to_check[n][b].min_rad * pixSpacing #convert to mm, curZ is also mm
            min_rad_ori = np.sqrt(min_rad**2 + curZ**2)
            if min_rad_ori < min_rad_from_ori[n]:
                min_rad_from_ori[n] = min_rad_ori
        
    return min_rad_from_ori