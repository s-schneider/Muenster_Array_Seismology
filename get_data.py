# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 15:12:53 2015

@author: Simon
"""

from __future__ import print_function
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import Stream
import obspy
import numpy as np

import Muenster_Array_Seismology as MAS
from Muenster_Array_Seismology import plot_gcp

def gcp_extrapol(lat1, lon1, lat2, lon2):
    """
    calc new coord on extrapolated gcp of lat1|lon1 to lat2|lon2
    """
    #convert degree to radian
    phi1=lat1 * np.pi/180
    rho1=lon1 * np.pi/180
    
    phi2=lat2 * np.pi/180
    rho2=lon2 * np.pi/180
    
    rho12 = np.abs(rho2 - rho1)

    #central angle between array and piercepoint
    sigma12 = np.arccos(np.sin(phi1) * np.sin(phi2) + np.cos(phi1) * np.cos(phi2) * np.cos( rho12 ))

    #calc azimuth of array
    azi1 = np.arctan( np.sin(rho12) / ( np.cos(phi1) * np.tan(phi2) - np.sin(phi1)*np.cos(rho12)) )

    #calc azimuth at meridian 0
    azi0 = np.arcsin( np.sin(azi1)*np.cos(phi1) )

    #central angle between array and meridian gcp-extrapolation
    sigma01=np.arctan( np.tan(phi1)/np.cos(azi1) )

    #calc meridian 0 longitude
    rho0= rho1 - ( np.arctan(np.sin(azi0) * np.sin(sigma01)/np.cos(sigma01)) )
    print(rho0 * 180/np.pi)
    print(rho1 * 180/np.pi)
    print(np.arctan(np.sin(azi0) * np.sin(sigma01)/np.cos(sigma01)) * 180/np.pi)
    #calc new central angle, double distance between array and piercepoint
    """
    Hier if abfrage wegen lat/lon1
    """
    sigmanew=2*sigma12+sigma01
    phinew = np.arcsin(np.cos(azi0) * np.sin(sigmanew));
    rhonew = np.arctan(abs(np.sin(azi0) * np.sin(sigmanew)/np.cos(sigmanew)) + rho0 )

    qlat = phinew * 180/np.pi
    qlon = rhonew * 180/np.pi
    
    return(qlat,qlon)


    
"""
Ask User to enter desired coordinates
"""

#bcoords = input("Please enter min. latitude, max. latitude, max. longitude and min. longitude")
#blatmin = bcoords.split()[0]
#blatmax = bcoords.split()[1]
#blonmin = bcoords.split()[2]
#blonmax = bcoords.split()[3]


"""
desired array coordinates and pierce-point
"""
slat = 10
slon = 10

plat = [0]
plon = [30]

def APsearch(slat,slon,plat,plon):
    """
    Input coordinates of an arrayregion and a region of interest of pierce- or
    reflectionpoints. It will search for matching arrays, times and events
    """    
    
    #################
    array_lat = []
    array_lon = []
    
    event_lat = []
    event_lon = []
    
    pierce_lat = []
    pierce_lon = []
    for i in range(len(plat)):    
        pierce_lat.append(plat[i])
        pierce_lon.append(plon[i])
        array_lat.append(slat)
        array_lon.append(slon)   

        qlat, qlon = gcp_extrapol(lat1=slat, lon1=slon, lat2=plat[i], lon2=plon[i])
        event_lat.append(qlat)
        event_lon.append(qlon)
    #################
    plot_gcp(slat=array_lat, slon=array_lon, qlat=event_lat, qlon=event_lon,
          plat=pierce_lat, plon=pierce_lon)

APsearch(slat,slon,plat,plon)

#print("Nice to meet you " + str(name) + "!")
#age = input("Your age? ")
#print("So, you are are already " + str(age) + " years old, " + name + "!")

"""
calculate scoord and qcoord couple
Use scoord to look for network and timespan
"""
"""
Use the calculated qlat and qlon and time to get data
"""
#client = Client("IRIS")
#t = UTCDateTime("2011-03-11T05:46:23")  # Tohoku
#catalog = client.get_events(starttime=t - 100, endtime=t + 24 * 3600,
#                            minmagnitude=7)
##catalog = client.get_events(starttime=None, endtime=None, minlatitude=None,
##                   maxlatitude=None, minlongitude=None, maxlongitude=None,
##                   latitude=None, longitude=None, minradius=None,
##                   maxradius=None, mindepth=None, maxdepth=None,
##                   minmagnitude=None)
#print(catalog)


