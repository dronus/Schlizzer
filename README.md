Schlizzer
=========

Early experimental Slic3r remake using C++

The reason for this remake is the low performance of Slic3r on embedded platforms like the RaspberryPi. 
I like to have a cheap embedded platform capable of quick slicing attached to the printer. 

Important features missing in respect to Slic3r:

- No discrimination of  'solid' and 'fill' areas, thus they are handled equally with fill_density=1 always.
- Only fill pattern 'rectlinear' available
- No contour correction of perimeter and infill, thus the object exceeds the specified .stl
- Bad route planning leading to bad movement order with useless many and long travels
- Even worse planning for 'non manifold' objects (most of thingyverse i guess..)
- No brim, skirt, cooling etc. 
- No automatic placement and z leveling
- No Gcode init header, layer, and footer, so generated Gcode need to be extended before printing

Current performance gain in respect to Slic3r: 
Note this is no fair comparison at the moment with so many features missing
- Object sliced:  http://www.thingiverse.com/download:59484 - a nice extruder body
- Performance gain on Linux / AMD64 platform: 0.3s vs 21.4s, about 70 times
- Performance gain on Linux / RaspberryPi ARMv6 platform:  5.9s vs 1270s, about 200 times