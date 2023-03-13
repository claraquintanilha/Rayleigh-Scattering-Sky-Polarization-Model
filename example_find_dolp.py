#!/usr/bin/env python3
'''
Given the azimuth and elevation of the light source (az0, el0), and
the azimuth and elevation of the viewing direction (az, el), along 
with the mesh refinement parameter npts, estimate the percent DoLP
at the viewing location based on the tables in Natraj et al. (2009).
'''

# Import stuff from natraj_table_utils.py
import natraj_table_utils as nt

# Input tau (my_tau). Has to be one of [0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1]
# Input albedo (my_alb). Has to be one of 0, 0.25, 0.8
my_tau, my_alb = 0.15, 0.25

# Azimuth and elevation of light source in degrees
az0, el0 = 20, 10

# Azimuth and elevation of viewing direction in degrees
az, el = 200, 80


dolp = nt.find_dolp (az0, el0, az, el, my_tau, my_alb, 99)
print('\nFor (az0,el0)=(%.1f, %.1f) (az,el)=(%.1f, %.1f) : DoLP_pct= %.1f\n' \
        % (az0, el0, az, el, dolp) )

