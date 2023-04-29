'''
Sky Model of Linear Polarization due to Rayleigh Scattering
Written by Clara Quintanilha // Wheaton College '23

This program takes
    Input: Tau, Albedo, and Mu_0 values
to interpolate between Stokes parameter values I/Q/U as tabulated by Natraj et al. (2009, ApJ, 691, 1909) 
(the files of which can be downloaded from https://web.gps.caltech.edu/~vijay/Rayleigh_Scattering_Tables/ )
to calculate values of DoLP and AoLP for each point in the sky, ultimately producing 
    Output: Maps displaying DoLP and AoLP values for all locations in the sky above horizon

This program only models the light source (sun/moon) along one straight line, however code posted to our GitHub written by 
Professor Dipankar Maitra improves upon this limitation and can model the light source's movement with more degrees of freedom.
'''

## Import Libraries

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

def main():
    
    ## Input
    # Possible options:
        # tau = 0.02, 0.05, 0.1, 0.15, 0.25, 1
        # a = 0, 0.25, 0.8
        # mn = 0.1, 0.2, 0.4, 0.6, 0.8, 0.92, 1
        
    tau = 0.25
    a = 0
    mn = 0.8
    
    ## Define Functions
    
    def readfile(var,tau):
        filename = f"{var}_{tau}.txt"
        fileopen = open(filename,"r")
        fileread = fileopen.readlines()
        fileopen.close()
        return fileread
    
    def convert(var):
        var_list = []
        for i in var:
            j = " ".join(i.split())
            k = list(j.split(" "))
            y = []
            for x in k:
                z = float(x)
                y.append(z)
            var_list.append(y)
        newvar = np.array(var_list)
        newvar = newvar[:,2:]
        return newvar
    
    ## Selecting files according to Tau
    
    U = readfile("U",tau)
    Q = readfile("Q",tau)
    I = readfile("I",tau)
    
    ## Select Lines according to Albedo
    
    if a == 0:
        min_range = 7
        max_range = 119
    elif a == 0.25:
        min_range = 122
        max_range = 234
    elif a == 0.80:
        min_range = 237
        max_range = 349
    
    U = U[min_range:max_range]
    Q = Q[min_range:max_range]
    I = I[min_range:max_range]
    
    ## Convert Elements to Floats, & Lists to Arrays
    
    U = convert(U)
    Q = convert(Q)
    I = convert(I)
    
    ## Select Lines according to Mu0
    
    if mn == 0.1:
        min_range2 = 0
        max_range2 = 16
    elif mn == 0.2:
        min_range2 = 16
        max_range2 = 32
    elif mn == 0.4:
        min_range2 = 32
        max_range2 = 48
    elif mn == 0.6:
        min_range2 = 48
        max_range2 = 64
    elif mn == 0.8:
        min_range2 = 64
        max_range2 = 80
    elif mn == 0.92:
        min_range2 = 80
        max_range2 = 96
    elif mn == 1:
        min_range2 = 96
        max_range2 = 112
    
    U = U[min_range2:max_range2]
    Q = Q[min_range2:max_range2]
    I = I[min_range2:max_range2]
    
    ## Interpolate
    
    intType = 'cubic'
    phi = np.array([0,30,60,90,120,150,180])
    mu = np.array([0.02, 0.06, 0.1, 0.16, 0.2, 0.28, 0.32, 0.4, 0.52, 0.64, 0.72, 0.84, 0.92, 0.96, 0.98, 1])
    
    fI = interp2d(phi, mu, I, kind=intType)
    fQ = interp2d(phi, mu, Q, kind=intType)
    fU = interp2d(phi, mu, U, kind=intType)
    
    xnew = np.linspace(phi[0], phi[-1], 99)
    ynew = np.linspace( mu[0],  mu[-1], 99)
    
    Inew = fI(xnew, ynew)
    Qnew = fQ(xnew, ynew)
    Unew = fU(xnew, ynew)
    
    DoLP = 100*np.sqrt(Qnew**2 + Unew**2)/Inew
    
    # np.arctan2 calculates x1/x2 while choosing the correct quadrant
    AoLP = np.rad2deg((np.arctan2(Unew,Qnew))/2)
    
    # Plotting
    
    X, Y = np.meshgrid(xnew, ynew)
    
    fig, (ax1, ax2) = plt.subplots(1,2, subplot_kw=dict(projection='polar'), figsize = (10,5))
    fig.suptitle(f"DoLP and AoLP maps of the sky for values $\mu_0$ = {mn}, $\u03c4$ = {tau}, Albedo = {a}.\nThe white dot indicates the location of the Sun, and the white stars represent neutral points.", fontsize=14)
    tticks = np.arange(0,360,45)
    ax1.set_rticks([20, 40, 60, 80])
    ax1.set_thetagrids(tticks, ('W', 'NW', 'N', 'NE', 'E', 'SE', 'S', 'SW'))
    ax2.set_rticks([20, 40, 60, 80])
    ax2.set_thetagrids(tticks, ('W', 'NW', 'N', 'NE', 'E', 'SE', 'S', 'SW'))
    
    zdsun = np.rad2deg(np.arccos(mn))
    azm = np.deg2rad(xnew)                
    rad = np.rad2deg(np.arccos(ynew))   
    th, r = np.meshgrid(azm, rad)
    
    sun_marker = 'wo'
    sun_markersize = 10
    ax1.grid(False)
    bar1 = ax1.pcolor(th, r, DoLP)
    ax1.pcolor(-th, r, DoLP)
    ax1.plot([0],[zdsun], sun_marker, markersize=sun_markersize)
    ax1.set_title('DoLP')
    ax1.grid()
    cbar1 = fig.colorbar(bar1, ax=ax1, orientation="vertical", fraction=0.04, pad=0.2)
    cbar1.set_label('DoLP (%)', rotation=90, labelpad=-65)
    
    ax2.grid(False)
    c = 'hsv'
    bar2 = ax2.pcolor(th, r, AoLP,cmap=c)
    ax2.pcolor(-th, r, -AoLP, cmap=c)
    ax2.plot([0],[zdsun], sun_marker, markersize=sun_markersize)
    ax2.set_title('AoLP')
    ax2.grid()
    cbar2 = fig.colorbar(bar2, ax=ax2, orientation="vertical", fraction=0.04, pad=0.2)
    cbar2.set_label('AoLP (degrees)', rotation=90, labelpad=-65) 
    
    #Loosely approximate location of neutral points
    zdnew = np.rad2deg(np.arccos(ynew))
    meridian_solar = DoLP[:,0]
    meridian_antisolar = DoLP[:,-1]
    
    zdnew = zdnew.tolist() # convert vars to lists for indexing
    meridian_solar = meridian_solar.tolist()
    meridian_antisolar = meridian_antisolar.tolist()
    
    sunloc = zdnew.index(zdsun) # find index locations of sun & neutral points
    zero_babinet = min(meridian_solar[:sunloc:]) # neutral points are not quite zero, so look for mins
    zero_brewster = min(meridian_solar[sunloc::])
    zero_arago = min(meridian_antisolar)
    
    zdloc_babinet = meridian_solar.index(zero_babinet) # find zd index for neutral points
    zdloc_brewster = meridian_solar.index(zero_brewster)
    zdloc_arago = meridian_antisolar.index(zero_arago)
    
    babinet = zdnew[zdloc_babinet] # find zd coordinate for neutral points
    brewster = zdnew[zdloc_brewster]
    arago = zdnew[zdloc_arago]
    
    np_marker = 'w*'
    np_size = 10
    ax1.plot([0],[brewster], np_marker, markersize = np_size)
    ax2.plot([0],[brewster], np_marker, markersize = np_size)
    if babinet < max(rad) and zero_babinet < meridian_solar[sunloc]:
        ax1.plot([0],[babinet], np_marker, markersize = np_size)
        ax2.plot([0],[babinet], np_marker, markersize = np_size)
    if arago < max(rad) and arago > 0:
        ax1.plot([np.pi],[arago], np_marker, markersize = np_size)
        ax2.plot([np.pi],[arago], np_marker, markersize = np_size)

    fig.tight_layout(pad=1.0)

if __name__ == "__main__":
    main()
