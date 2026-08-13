"""Rotinas básicas para aplicações GNSS em Geodésia."""
import numpy as np
from .coordinates import ecef_to_enu
from .heights import ellipsoidal_to_orthometric_height
def geometric_range(receiver_xyz,satellite_xyz): return np.linalg.norm(np.asarray(satellite_xyz,dtype=float)-np.asarray(receiver_xyz,dtype=float),axis=0)
def pseudorange_model(receiver_xyz,satellite_xyz,receiver_clock_bias=0.0,satellite_clock_bias=0.0,c=299792458.0): return geometric_range(receiver_xyz,satellite_xyz)+c*(receiver_clock_bias-satellite_clock_bias)
def satellite_elevation_azimuth(receiver_lat,receiver_lon,receiver_h,sat_X,sat_Y,sat_Z,degrees=True):
    E,N,U=ecef_to_enu(sat_X,sat_Y,sat_Z,receiver_lat,receiver_lon,h0=receiver_h,degrees=degrees); hor=np.sqrt(E*E+N*N); el=np.arctan2(U,hor); az=np.arctan2(E,N)%(2*np.pi); return (np.rad2deg(el),np.rad2deg(az)) if degrees else (el,az)
def simple_tropospheric_delay(elevation_deg,zenith_delay=2.3): return zenith_delay/np.maximum(np.sin(np.deg2rad(elevation_deg)),1e-3)
def gnss_height_to_orthometric(h_ellipsoidal,geoid_height): return ellipsoidal_to_orthometric_height(h_ellipsoidal,geoid_height)
