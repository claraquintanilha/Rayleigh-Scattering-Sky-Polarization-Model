#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions to read tables in Natraj et al. (2009, ApJ, 691, 1909) 
tables with I/Q/U and interpolates them to obtain I/Q/U/DoLP values at 
arbitrary viewing directions in the sky, for light source in the 
mu0 range of [0.1:1], i.e., the zenith distance of the light source in the 
range of [84.26:0] degrees.

The observing direction can have zenith distance range of [88.85:0] degrees.

2022-Jan-30: First version (DM)
"""

#######################################################################
# Stuff to import
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from scipy.interpolate import RegularGridInterpolator as rgi

import sys
from os.path import exists
from urllib import request
from urllib.parse import urljoin
#######################################################################

#######################################################################
# The grid of model parameters that are tabulated in the Natraj et al. paper

# Tabulated tau values
tau = np.array([0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1]) 

# Ground albedo should be one of 0, 0.25, 0.8

# Tabulated mu0 values
mu0 = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 0.92, 1]) 

# Tabulated mu and phi values
mu = np.array([0.02, 0.06, 0.1, 0.16, 0.2, 0.28, 0.32, 0.4, 0.52, \
        0.64, 0.72, 0.84, 0.92, 0.96, 0.98, 1.])
phi = np.array([0,30,60,90,120,150,180])
#######################################################################


#######################################################################
# Getting tables from web repository
#
# Probably not the safest way to read the tables from the web  ...
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

# The data table files are here
tbl_url = 'https://web.gps.caltech.edu/~vijay/Rayleigh_Scattering_Tables/CDS/'
#######################################################################


#######################################################################
# Given tau, this function checks which data table files we need, and
# downloads them if they are not in the current working directory.
def get_files(tau):
    Ifile = "I_DN_TAU_%s" % tau
    Qfile = "Q_DN_TAU_%s" % tau
    Ufile = "U_DN_TAU_%s" % tau
    flist = [Ifile, Qfile, Ufile]
    print('Status of data files to read:')
    
    for file_name in flist:
        if exists(file_name):
            print(' ', file_name, 'is already downloaded')
        else:
            file_url = urljoin(tbl_url, file_name)
            request.urlretrieve(file_url, file_name)
            print('getting',file_name,'from', file_url)
            if not exists(file_name):
                print('could not get file from', file_url)
                sys.exit(0)

    return flist
#######################################################################


#######################################################################
# Given albedo, this function returns which parts of the 
# datafiles to read.
# Albedo should be one of 0, 0.25, 0.8, and reading the data for these 
# albedo require skipping 7, 122, 237 lines in the table data files.
def get_lines2skip (albedo):
    if albedo == 0:
        a_skip = 7

    elif albedo == 0.25:
        a_skip = 122

    elif albedo == 0.8:
        a_skip = 237

    else:
        print('Albedo has to be one of 0, 0.25, or 0.8')
        sys.exit(0)

    return a_skip
#######################################################################


#######################################################################
# Given albedo, tau, and a list of mu0, mu, phi values, this function 
# returns interpolated values of I/Q/U/DoLP at the listed (mu0,mu,phi) points.
def ReadTableLoadRefinedGrid (my_tau, my_alb, nRefine, points_list):
    '''

    Parameters
    ----------
    my_tau : optical depth of atmosphere.

    my_alb : ground albedo.
        
    nRefine : refine the mu and phi table grids, using cubic spline 
            interpolation so that the refined grid for each mu and phi has 
            nRefine grid points.
            
    points_list: list of mu0, mu, phi values where interpolated values of
            I/Q/U are to be output.
        
        
    Returns
    -------
    Lists of interpolated values of I/Q/U/DoLP at the locations defined 
    in points_list.
    '''
    
    # Refined mu-phi grid at each tabulated mu0, and set the 2D 
    # interpolation type to be cubic
    muRefined = np.linspace( mu[0],  mu[-1], nRefine)
    phiRefined = np.linspace(phi[0], phi[-1], nRefine)
    TwoDintType='cubic'

    nmu0, nmu, nphi = len(mu0), len(mu), len(phi)
    nrows = nmu0 * nmu

    # Given tau, get the right data files
    (Ifile, Qfile, Ufile) = get_files(my_tau)

    # Given albedo, read the right portion of the data files
    skp = get_lines2skip (my_alb)

    # From the I/Q/U datafiles read the relevant data
    I_table = np.loadtxt(Ifile, skiprows=skp, max_rows=nrows, 
                   usecols = (2,3,4,5,6,7,8))
    Q_table = np.loadtxt(Qfile, skiprows=skp, max_rows=nrows, 
                   usecols = (2,3,4,5,6,7,8))
    U_table = np.loadtxt(Ufile, skiprows=skp, max_rows=nrows, 
                   usecols = (2,3,4,5,6,7,8))

    # Make 3-dimensional grid of table data
    I3d = I_table.reshape(nmu0, nmu, nphi)
    Q3d = Q_table.reshape(nmu0, nmu, nphi)
    U3d = U_table.reshape(nmu0, nmu, nphi)

    # Create 3D grids or arrays to store the refined table values
    I3Refined = np.zeros((nmu0, nRefine, nRefine))
    Q3Refined = np.zeros((nmu0, nRefine, nRefine))
    U3Refined = np.zeros((nmu0, nRefine, nRefine))
 
    # Refine the mu-phi grid for each tabulated mu0
    for ii in range(nmu0):
        
        # Take slice in mu0
        I, Q, U  = I3d[ii], Q3d[ii], U3d[ii]
        
        # Set up the interpolation functions for that slice
        fI = interp2d(phi, mu, I, kind=TwoDintType)
        fQ = interp2d(phi, mu, Q, kind=TwoDintType)
        fU = interp2d(phi, mu, U, kind=TwoDintType)

        # Get the interpolated I/Q/U values, and the DoLPs,
        # at the new grid points, and stack different mu0 vertically
        I3Refined[ii] = fI(phiRefined, muRefined)
        Q3Refined[ii] = fQ(phiRefined, muRefined)
        U3Refined[ii] = fU(phiRefined, muRefined)
        
    # Set up the 3D interpolator, which unfortunately only accepts
    # 'linear' or 'nearest' as methods. Cubic spline would be nice!
    fI = rgi ((mu0, muRefined, phiRefined), I3Refined, method='linear')
    fQ = rgi ((mu0, muRefined, phiRefined), Q3Refined, method='linear')
    fU = rgi ((mu0, muRefined, phiRefined), U3Refined, method='linear')
        
    # Get the list of I/Q/U/DoLP values at the requested points
    Ilist = fI(points_list)
    Qlist = fQ(points_list)
    Ulist = fU(points_list)
    
    # Compute the DoLPs at the requested points too
    DoLPs = 100 * np.sqrt(Qlist*Qlist + Ulist*Ulist) / Ilist
    
    # Return the interpolated values of I/Q/U/DoLP at the list of points
    return Ilist, Qlist, Ulist, DoLPs
#######################################################################


#######################################################################
def locate_neutral_points (mu0, muArr, DoLP):
    '''
    Search for neutral points (NP) in the sky.
    NOTE: Erroneous results for ZD~<25 deg, i.e. mu>0.9? Return that
    no neutral points found if mu>0.9 !!!

    Parameters
    ----------
    mu0 : TYPE: float.
        DESCRIPTION: cos(zenith distance) to Sun/Moon.
    muArr : TYPE: 1-D float array. 
        DESCRIPTION: mu values of observing directions.
    DoLP : TYPE: 2-D float array.
        DESCRIPTION: DOLP values on the sky for given mu0.

    Returns
    -------
    A 1-D array with 3 elements. The first element is the zenith distance 
    to the Arago NP, the Babinet NP is the second element of the 
    array, and the Brewster NP the third element of the array.
    The units are in degrees.
    Zenith distances in the antisolar meridian (opposite side of the
    meridian as the Sun) have negative values
    Missing NPs are assigned a physically implausible zenith 
    distance value of -999.
    '''
    
    zd_np = np.full(3,-999, dtype=float)   # Array to store NPs.
    
    if mu0 > 0.90:       # See NOTE above
        print('Source too close to zenith. No good neutral points.')
        return zd_np
    
    zdsun = np.rad2deg( np.arccos(mu0) )
    zdArr = np.rad2deg( np.arccos(muArr) )
    meridian_solar = DoLP[:,0]
    meridian_antisolar = DoLP[:,-1]

    # Search for a NP in the antisolar meridian direction
    minLoc = np.argmin(meridian_antisolar)

    if minLoc > 0 and minLoc < len(meridian_antisolar) - 1:
        zd_min_as = -zdArr[minLoc]
        if zdsun - zd_min_as > 90:
            ptName = "Arago"
            degAboveAntisun = 180 - zdsun + zd_min_as
            print('Found',ptName,'point',degAboveAntisun,
                    'deg above antisun.')
            zd_np[0] = zd_min_as
        else:
            ptName = "Babinet"
            degAboveSun = zdsun - zd_min_as
            print('Found',ptName,'point',degAboveSun,
                    'deg above Sun.')
            zd_np[1] = zd_min_as

    # Search for Babinet and Brewster points in the solar meridian direction
    idBelowSun = np.searchsorted(muArr, mu0)

    belowSun = meridian_solar[:idBelowSun]
    zdBelowSunArr = zdArr[:idBelowSun]
    minLoc = np.argmin(belowSun)
    if minLoc > 0 and minLoc < len(belowSun) - 1:
        zdMin = zdBelowSunArr[minLoc]
        ptName = "Brewster"
        degBelowSun = zdMin - zdsun
        print('Found',ptName,'point',degBelowSun,
                'deg below Sun.')
        zd_np[2] = zdMin

    aboveSun = meridian_solar[idBelowSun:]
    zdAboveSunArr = zdArr[idBelowSun:]
    minLoc = np.argmin(aboveSun)
    if minLoc > 0 and minLoc < len(aboveSun) - 1:
        zdMin = zdAboveSunArr[minLoc]
        ptName = "Babinet"
        degAboveSun = zdsun - zdMin
        print('Found',ptName,'point',degAboveSun,
                'deg above Sun.')
        zd_np[1] = zdMin

    return zd_np
#######################################################################


#######################################################################
def find_dolp (az0, el0, az, el, tau, albedo, npts):
    '''
    Given the azimuth and elevation of the light source (az0, el0), and
    the azimuth and elevation of the viewing direction (az, el), along 
    with the mesh refinement parameter npts, estimate the percent DoLP
    at the viewing location based on the tables in Natraj et al. (2009).
    '''

    # Convert the elevations into mu-s
    mu0 = np.cos( np.deg2rad(90-el0) )
    mu  = np.cos( np.deg2rad(90-el ) )

    # Compute relative azimuth of viewing direction w.r.t. light source
    rel_az = np.absolute(az - az0)
    if rel_az > 180:
        rel_az = 360 - rel_az

    pt_list = np.array([mu0, mu, rel_az])
    
    # Get the I/Q/U/DoLP value at the location in pt_list
    I, Q, U, DoLP = ReadTableLoadRefinedGrid (tau, albedo, npts, pt_list)

    return DoLP[0]

#######################################################################


#######################################################################
def make_skymap_dolp (el, az, tau, albedo, npts, opfilename="None"):
    '''
    Given the elevation and azimuth of the Sun/Moon, along with the mesh
    refinement parameter npts, creates a map of the DoLP distributions on the
    sky.

    Parameters
    ----------
    el : elevation of Sun/Moon in degrees
    az : azimuth of Sun/Moon in degrees
    tau: optical depth of the atmosphere. Tau has to be one of the following
         [0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1]
    albedo : ground albedo. Has to be one of [0, 0.25, 0.8]
    npts : Refine the grid such that each the azimuth grid and the elevation
           grid has npts points.
           
    opfilename : if supplied, then use this as the output filename.

    Returns
    -------
    A png image of the sky with DoLP values color coded.

    '''
    mu0 = np.cos( np.deg2rad(90-el) )
    
    # Check to make sure that mu0 >= 0.1 (smallest mu0 in tables)
    if mu0 < 0.1 :
        print('Too low elevation of light source.')
        sys.exit(0)
        
        
    # Regrid the mu and phi grids so that each have npts grid points
    muNew = np.linspace( mu[0],  mu[-1], npts)
    phiNew = np.linspace( phi[0],  phi[-1], npts)

    # Generate the list of points where we'd like to get I/Q/U/DoLP etc.
    # The process to generate the points list is based on the following
    # https://stackoverflow.com/questions/57659322/evaluate-a-regulargridinterpolator-on-an-another-regular-grid
    pts_grid = np.meshgrid([mu0], muNew, phiNew, indexing='ij')
    pts_list = np.reshape(pts_grid, (3, -1), order='C').T


    # Get the list of I/Q/U values at the locations on the list
    Ivals, Qvals, Uvals, DoLPvals = ReadTableLoadRefinedGrid (tau, 
                                            albedo, npts, pts_list)

    # Turn the DoLP values into 2D arrays for visualizations
    DoLP = DoLPvals.reshape(npts, npts)

    # Visualize results
    my_filetype = "png"
    my_dpi = 150
    if opfilename == "None" :
        skymap_dolp(tau, albedo, az, el, phiNew, muNew,
                    DoLP, my_filetype, my_dpi)
    else :
        skymap_dolp(tau, albedo, az, el, phiNew, muNew,
                    DoLP, my_filetype, my_dpi, opfilename=opfilename)
#######################################################################



#######################################################################
# Visualizations

# Detailed plots of I, Q, U, DoLP skymap etc
def make_detailed_plots (tau, albedo, mu0, phiArr, muArr, I, Q, U, DoLP,
                         opfiletype, opdpi):
    '''
    Produces a detailed visualization showing I/Q/U/DoLP etc.
    
    Parameters
    ----------
    tau : optical depth of atmosphere
        
    albedo : ground albedo
        
    mu0 : cos(zenith-distance to light source)
        
    phiArr : viewing direction azimuth grid points
        
    muArr : viewing direction zenith-distance grid points
        
    I/Q/U/DoLP: 2D arrays with I/Q/U/DoLP values
    
    opfiletype : output file type (png/pdf/jpg)
    
    opdpi : DPI of output file

    Returns
    -------
    PDF file with visualizations.
    '''
    
    # Create the name of the output graphics file
    op = 'tau_%.2f_A_%.2f_mu0_%.2f.%s' % (tau, albedo, mu0, opfiletype)
    print('Output graphics file:',op)

    X, Y = np.meshgrid(phiArr, muArr)

    zdsun = np.rad2deg( np.arccos(mu0) )
    zdArr = np.rad2deg( np.arccos(muArr) )
    meridian_solar = DoLP[:,0]
    meridian_antisolar = DoLP[:,-1]
    
    # Get locations of neutral points
    neutral_pts = locate_neutral_points (mu0, muArr, DoLP)
    zd_nps = neutral_pts[neutral_pts > -999] # Select only detected NPs
    az_nps_rad = np.where(zd_nps>0,0, np.pi) # Azimuths of detected NPs (rad)
    az_nps_deg = np.where(zd_nps>0,0, 180)   # Azimuths of detected NPs (deg)
    zd_nps = np.absolute(zd_nps)             # Zenith distances of NPs (deg)
    mu_nps = np.cos(np.deg2rad(zd_nps))


    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(15, 9))
    fig.suptitle(r'Interpolations for $(\mu_0, \tau, A)=$ (%.2f, %.2f, %.2f).'\
                 ' The red/white dots show the location of the Sun/neutral points.' % (mu0, tau, albedo) , fontsize=14)

    ax = axs[0, 0]
    c = ax.pcolor(X, Y, I, cmap='plasma')
    ax.plot([0],[mu0], 'ro', markersize=10)
    ax.plot(az_nps_deg, mu_nps, 'wo', markersize=5)
    ax.set_title('I')
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\mu$')
    fig.colorbar(c, ax=ax)

    ax = axs[0, 1]
    c = ax.pcolor(X, Y, Q, cmap='plasma')
    ax.plot([0],[mu0], 'ro', markersize=10)
    ax.plot(az_nps_deg, mu_nps, 'wo', markersize=5)
    ax.set_title('Q')
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\mu$')
    fig.colorbar(c, ax=ax)

    ax = axs[0, 2]
    c = ax.pcolor(X, Y, U, cmap='plasma')
    ax.plot([0],[mu0], 'ro', markersize=10)
    ax.plot(az_nps_deg, mu_nps, 'wo', markersize=5)
    ax.set_title('U')
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\mu$')
    fig.colorbar(c, ax=ax)

    ax = axs[1, 0]
    ax.set_xlabel("Zenith distance in solar/antisolar direction in green/blue. \
        \nSun's location shown by vertical red line.")
    ax.set_ylabel('DoLP')
    ax.plot(zdArr, meridian_solar, 'g-')
    ax.plot(-zdArr, meridian_antisolar, 'b-')
    ax.axvline(zdsun, color ='r')
    
    #ax.set(xlim =(-10, zdsun+10), ylim=(-2, 10)) # Uncomment to zoom near Sun
    ymin, ymax = ax.get_ylim()
    tpos = ymax - 0.3*(ymax-ymin)

    for ii in range(3):
        if neutral_pts[ii] !=-999:
            ax.axvline(neutral_pts[ii], color='k', linestyle='dotted')

    if neutral_pts[0] !=-999:
        ax.text(neutral_pts[0], tpos, 'Arago',    rotation=90,
                rotation_mode='anchor')
    if neutral_pts[1] !=-999:
        ax.text(neutral_pts[1], tpos, 'Babinet',  rotation=90, 
                rotation_mode='anchor')
    if neutral_pts[2] !=-999:
        ax.text(neutral_pts[2], tpos, 'Brewster', rotation=90,
                rotation_mode='anchor')
    ax.grid()

    ax = axs[1, 1]
    c = ax.pcolor(X, Y, DoLP)
    ax.plot([0],[mu0], 'ro', markersize=10)
    ax.plot(az_nps_deg, mu_nps, 'wo', markersize=5)
    ax.set_title('DoLP')
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\mu$')
    fig.colorbar(c, ax=ax)

    # The final plot will show distribution of DoLP on the
    # sky, using a polar plot. So we first remove the linear
    # matplotlib axis at row=2, col=3, and instead add
    # a new axis there with polar projection
    axs[1,2].remove()
    ax = fig.add_subplot(2, 3, 6, projection='polar')
    ax.grid(False)
    azm = np.deg2rad(phiArr)
    rad = np.rad2deg( np.arccos(muArr) )
    th, r = np.meshgrid(azm, rad)
    ax.pcolormesh(th, r, DoLP)
    ax.pcolormesh(-th, r, DoLP)
    ax.plot([0],[zdsun], 'ro', markersize=10)       # Plot Sun
    ax.plot(az_nps_rad, zd_nps, 'wo', markersize=5) # Plot neutral points
    ax.grid()

    fig.tight_layout(pad=1.0)
    plt.savefig(op, dpi=opdpi)
    plt.close()
    
    return 0

# Skymap of DoLPs given azimuth and altitude of light source
def skymap_dolp (tau, albedo, az, el, phiArr, muArr, DoLP,
                         opfiletype, opdpi, opfilename="None"):
    '''
    Produces a detailed visualization showing I/Q/U/DoLP etc.
    
    Parameters
    ----------
    tau : optical depth of atmosphere
        
    albedo : ground albedo
    
    az : azimuth of light source in degrees
        
    el : elevation of light source in degrees
        
    phiArr : viewing direction azimuth grid points
        
    muArr : viewing direction zenith-distance grid points
        
    DoLP: 2D array with DoLP values
    
    opfiletype : output file type (png/pdf/jpg)
    
    opdpi : DPI of output file
    
    opfilename : if supplied, then use this as the output filename.

    Returns
    -------
    Graphics file with visualizations.
    '''
    
    # Create the name of the output graphics file
    if opfilename == "None" :
        op = 'tau_%.2f_A_%.2f_az_%06.2f_el_%05.2f.%s' \
            % (tau, albedo, az, el, opfiletype)
    else :
        op = '%s.%s' % (opfilename, opfiletype)
        
    print('Output graphics file:',op)


    plt_title = r'DoLP map for (elevation, azimuth, $\tau$, albedo) = (%.1f, %.1f, %.2f, %.2f)' \
                 % (el, az, tau, albedo)
    plt_title = plt_title + \
        '\n Light source as red point. Neutral points in white.'
    fig = plt.figure(figsize=(9, 9))
    fig.suptitle(plt_title , fontsize=12)

    ax = fig.add_subplot(projection='polar')
    
    zdsun = 90 - el
    mu0 = np.cos(np.deg2rad(zdsun))
    zdArr = np.rad2deg( np.arccos(muArr) )
    
    # Get locations of neutral points
    neutral_pts = locate_neutral_points (mu0, muArr, DoLP)
    zd_nps = neutral_pts[neutral_pts > -999] # Select only detected NPs
    az_nps_rad = np.where(zd_nps>0,0, np.pi) # Azimuths of detected NPs (rad)
    zd_nps = np.absolute(zd_nps)             # Zenith distances of NPs (deg)

    
    azm = np.deg2rad(phiArr)
    th, r = np.meshgrid(azm, zdArr)
    to_az = np.pi/2 + np.deg2rad(az)   # Rotate to appropriate azimuth
    ax.plot([to_az],[zdsun], 'ro', markersize=10)
    ax.plot(to_az + az_nps_rad, zd_nps, 'wo', markersize=5)
    ax.pcolormesh(to_az - th, r, DoLP, vmin=0, vmax=100)
    c = ax.pcolormesh(to_az + th, r, DoLP, vmin=0, vmax=100)

    # Customize radial and azimuthal ticks    
    ax.set_rticks([20, 40, 60, 80])
    tticks = np.arange(0,360,45)
    ax.set_thetagrids(tticks, ('W', 'NW', 'N', 'NE', 'E', 'SE', 'S', 'SW'))

    ax.grid(True)
    
    # Add colorbar
    cbar = fig.colorbar(c, ax=ax, orientation="vertical", fraction=0.04, \
                        pad=0.12)
    cbar.set_label('DoLP (%)', rotation=90, labelpad=-65)

    fig.tight_layout(pad=1.0)
    plt.savefig(op, dpi=opdpi)
    plt.close()
    
    return 0
#######################################################################

