#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 13:53:06 2022

@author: pstijnm2
"""
import numpy as np
from skimage.segmentation import flood_fill
class boundary:
    def __init__(self,points):
        self.points = points

    def fill(self, image):
        self.mask = np.ones(image.shape[0:2]) # create 2d mask
        for point in self.points:
            self.mask[point[0], point[1]] = 2
            
        self.mask = flood_fill(self.mask, (0,0), 0, selem=np.array([0,1,0,1,0,1,0,1,0]).reshape([3,3]))
        
        for point in self.points:
            self.mask[point[0], point[1]] = 1
        
        

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
def get_boundaries_from_points(points):
    xs = [-1,0,1]
    ys = [-1,0,1]
    boundaries = list()
    
    while (len(points) > 0):
        connected = get_connected_points(points,xs,ys)
        boundaries.append(boundary(connected))
        points = list(set(points) ^ set(connected)) # returns list of points not in connected
        
    return boundaries

                
'''
if a -> [r,g,b] == b -> [r,g,b] return 1
'''
def compare_rgb(a,b):
    return np.prod(a == b)

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
    green = [0,60,0]
    
    # 5mm
    line_red = get_points_rgb_color(image, red)
    boundaries_red = get_boundaries_from_points(line_red)
    # 3mm
    line_yellow = get_points_rgb_color(image, yellow)
    boundaries_yellow = get_boundaries_from_points(line_yellow)
    # 2mm 
    line_teal = get_points_rgb_color(image, teal)
    boundaries_teal = get_boundaries_from_points(line_teal)
    # 1mm
    line_green = get_points_rgb_color(image, green)
    boundaries_green = get_boundaries_from_points(line_green)
    
    return boundaries_green, boundaries_teal, boundaries_yellow, boundaries_red
    
def fill_boundaries(image, boundaries_green, boundaries_teal, boundaries_yellow, boundaries_red):
    for b in boundaries_red:
        b.fill(image)
    for b in boundaries_yellow:
        b.fill(image)    
    for b in boundaries_teal:
        b.fill(image)    
    for b in boundaries_red:
        b.fill(image)    
        
    return