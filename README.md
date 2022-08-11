# MR_GeomAcc

## Summary
This module performs analysis of results of measurements with the Philips single-slab Geometric Fidelity Phantom. 

## Status
Initial version released (20220810)

## Dependencies
pip:
- numpy
- pydicom
- matplotlib
- skimage

## Acquisition protocol
Use the recommended Philips protocol (Geometric Fidelity)

## Selector
A selector for this module should run on dcm_series with SeriesDescription: Geometric Distortion_Rgb

### Input data
dcm_series

### Config
A simple config file is provided:
- mr_philips_geometric_fidelity_config.json

It provides offsets for the top of the phantom and the start of the table.

### Meta
- mr_philips_geometric_fidelity_meta.json

Currently, there are 3 limits defined: passing the 20, 34, 45 cm DSV constraint (True/False).

### Rules
StationName, SeriesDescription

## Analysis
Process Geometric Accuracy measurement
- Find isolines
	- Detect colors in RGB dicoms
- Create spheres from geometric center
	- Diameters: 20, 34, 45 cm
- Check if spheres are within isolines
	- 20 cm --> 1 mm isoline
	- 34 cm --> 2 mm isoline
	- 45 cm --> 5 mm isoline
- Extra statistics: 
	- total area of 1/2/3/5 mm accuracy
	- largest sphere (3D) possible within 1/2/5 mm isoline
- Output images:
	- Plot with 7 the slices with the sphere outlines

## Results
- Figure with the 7 slice positions and the sphere outlines
- Isolines passed
- Total area (m^2) within the isolines
- Largest sphere (3D) possible within 1/2/3/5 mm isolines