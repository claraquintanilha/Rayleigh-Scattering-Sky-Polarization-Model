#!/usr/bin/env python3
'''
This sample code shows how to get the azimuths and
elevations of the Sun from any given location on any
given day, and then create visualizations showing how 
the DoLP distribution on the sky changes throughout the
day. As an example, we show the the DoLP variation 
throughout the day on June 21 and Dec 21, as seen from
Wheaton College Observatory.

First the user selects a day and creates an array of times
when to calculate the Sun's positions on this day.
Then for each of these times the code computes the azimuth 
and elevation of the Sun. We only keep those instances when
the zenith distance is less than 84.26 degrees because of
the limitations of the tables in Natraj et al. (2009).

Finally for each of the good instances (when Sun's zenith
distance is less than 84.26 degrees) the code produces a
PNG image showing the distribution of DoLP in the sky.

Once the code has completed, the user can use a tool of their
choice, e.g. imagemagick's "convert" tool, to create a movie
out of the PNGs produced by the code.

2023-Feb-12: First usable version. (DM)
'''

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, get_sun, AltAz


# Import stuff from natraj_table_utils.py
import natraj_table_utils as nt

# Input tau (my_tau) has to be one of [0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1]
# Input albedo (my_alb) has to be one of 0, 0.25, 0.8
my_tau, my_alb = 0.15, 0.25


# Define the observer at Wheaton College Observatory (WCO)
wco = EarthLocation(lat=41.9656*u.deg, lon=-71.1836*u.deg, height=40*u.m)

# Select the day.
myday = Time('2023-03-21 00:00:00')
utcoffset = -4 * u.hour  # EDT=UTC-4 during summer


# Set up the array of times for computing solar positions and 
# output file names based on local time
t_utc = myday - utcoffset
hmin, hmax, hinc = 4, 21, 1/40    # Stepsize = 1/40 hour = 1.5 minutes 
npts = 1 + int( (hmax-hmin)/hinc )
dt = np.linspace(hmin, hmax, npts) * u.hour
my_utc_times = t_utc + dt
my_loc_times = my_utc_times + utcoffset
opfilenames = my_loc_times.strftime('%Y-%m-%d-%H-%M-%S')


# Set up the AltAz coordinate frame for the above time array
AltAzFrame = AltAz(obstime = my_utc_times, location = wco)

# Compute AltAz of Sun at these times
sunAltAz = get_sun(my_utc_times).transform_to(AltAzFrame)


# Filter arrays to keep quantities only when the zenith distance (ZD) to 
# the Sun is less than 84.26 degrees [so that cos(ZD) > 0.1].
altmin = 90 - 84.26

good_alt = sunAltAz.alt.deg [sunAltAz.alt.deg > altmin]
good_az  = sunAltAz.az.deg  [sunAltAz.alt.deg > altmin]
good_tim = my_loc_times     [sunAltAz.alt.deg > altmin]

ngood = len(good_alt)
nmesh = 599            # Refine the mu-phi grid 

# Create visualizations for each time
for ii in range(ngood):
    pct_complete = 100*ii/ngood
    print('*** Creating DoLP map for',good_tim[ii], \
            'Percent of job done =',pct_complete)
    nt.make_skymap_dolp(good_alt[ii], good_az[ii], my_tau, my_alb, nmesh,
                        opfilename=opfilenames[ii])

