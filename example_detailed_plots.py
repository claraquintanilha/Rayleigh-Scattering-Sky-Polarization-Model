#!/usr/bin/env python3

import numpy as np

import natraj_table_utils as nt  # Import stuff from natraj_table_utils.py

# Input tau (my_tau) has to be one of [0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1]
my_tau = 0.1

# Input albedo (my_alb) has to be one of 0, 0.25, 0.8
my_alb = 0.25

# Input mu0 (my_mu0) has to be within 0.1 <= my_mu0 <= 1 
my_mu0 = np.cos (np.deg2rad(46))

# Regrid the mu and phi grids so that each have nNew grid points
nNew = 499
muNew = np.linspace( nt.mu[0],  nt.mu[-1], nNew)
phiNew = np.linspace( nt.phi[0],  nt.phi[-1], nNew)

# Generate the list of points where we'd like to get I/Q/U/DoLP etc.
# The process to generate the points list is based on the following
# https://stackoverflow.com/questions/57659322/evaluate-a-regulargridinterpolator-on-an-another-regular-grid
pts_grid = np.meshgrid([my_mu0], muNew, phiNew, indexing='ij')
pts_list = np.reshape(pts_grid, (3, -1), order='C').T


# Get the list of I/Q/U values at the locations on the list
Ivals, Qvals, Uvals, DoLPvals = nt.ReadTableLoadRefinedGrid (my_tau, my_alb, 
                                                          nNew, pts_list)

# Turn the I/Q/U/DoLP values into 2D arrays for visualizations
I = Ivals.reshape(nNew, nNew)
Q = Qvals.reshape(nNew, nNew)
U = Uvals.reshape(nNew, nNew)
DoLP = DoLPvals.reshape(nNew, nNew)


# Visualize results
my_filetype = "png"
my_dpi = 150
nt.make_detailed_plots(my_tau, my_alb, my_mu0, phiNew, muNew, I, Q, U,
                    DoLP, my_filetype, my_dpi)

