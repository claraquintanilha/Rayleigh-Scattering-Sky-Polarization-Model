#!/usr/bin/env python3

# Import stuff from natraj_table_utils.py
import natraj_table_utils as nt

# Input tau (my_tau) has to be one of [0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1]
my_tau = 0.15

# Input albedo (my_alb) has to be one of 0, 0.25, 0.8
my_alb = 0.25

# Azimuth and elevation in degrees
my_az, my_el = 40, 10

nt.make_skymap_dolp(my_el, my_az, my_tau, my_alb, 299)
