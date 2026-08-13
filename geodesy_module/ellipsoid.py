"""Geometria do elipsoide de revolução."""
import numpy as np
from .constants import WGS84, DEG2RAD, RAD2DEG
def semi_minor_axis(a=WGS84.a,f=WGS84.f): return a*(1.0-f)
def first_eccentricity_squared(f=WGS84.f): return 2.0*f-f**2
def first_eccentricity(f=WGS84.f): return np.sqrt(first_eccentricity_squared(f))
def second_eccentricity_squared(f=WGS84.f):
    e2=first_eccentricity_squared(f); return e2/(1.0-e2)
def second_eccentricity(f=WGS84.f): return np.sqrt(second_eccentricity_squared(f))
def linear_eccentricity(a=WGS84.a,f=WGS84.f):
    b=semi_minor_axis(a,f); return np.sqrt(a*a-b*b)
def prime_vertical_radius(lat,a=WGS84.a,f=WGS84.f,degrees=True):
    lat=np.asarray(lat,dtype=float); lat=lat*DEG2RAD if degrees else lat
    e2=first_eccentricity_squared(f); return a/np.sqrt(1.0-e2*np.sin(lat)**2)
def meridian_radius(lat,a=WGS84.a,f=WGS84.f,degrees=True):
    lat=np.asarray(lat,dtype=float); lat=lat*DEG2RAD if degrees else lat
    e2=first_eccentricity_squared(f); return a*(1.0-e2)/(1.0-e2*np.sin(lat)**2)**1.5
def mean_radius_of_curvature(lat,a=WGS84.a,f=WGS84.f,degrees=True):
    return np.sqrt(meridian_radius(lat,a,f,degrees)*prime_vertical_radius(lat,a,f,degrees))
def geocentric_radius(lat,a=WGS84.a,f=WGS84.f,degrees=True):
    lat=np.asarray(lat,dtype=float); lat=lat*DEG2RAD if degrees else lat
    b=semi_minor_axis(a,f)
    num=(a*a*np.cos(lat))**2+(b*b*np.sin(lat))**2
    den=(a*np.cos(lat))**2+(b*np.sin(lat))**2
    return np.sqrt(num/den)
def geocentric_latitude(geodetic_lat,f=WGS84.f,degrees=True):
    phi=np.asarray(geodetic_lat,dtype=float); phi=phi*DEG2RAD if degrees else phi
    e2=first_eccentricity_squared(f); th=np.arctan((1.0-e2)*np.tan(phi))
    return th*RAD2DEG if degrees else th
def reduced_latitude(geodetic_lat,f=WGS84.f,degrees=True):
    phi=np.asarray(geodetic_lat,dtype=float); phi=phi*DEG2RAD if degrees else phi
    beta=np.arctan((1.0-f)*np.tan(phi)); return beta*RAD2DEG if degrees else beta
def authalic_radius(a=WGS84.a,f=WGS84.f):
    e=first_eccentricity(f); area=2*np.pi*a*a*(1+(1-e*e)/e*np.arctanh(e)); return np.sqrt(area/(4*np.pi))
def ellipsoid_summary(a=WGS84.a,f=WGS84.f):
    b=semi_minor_axis(a,f); e2=first_eccentricity_squared(f); ep2=second_eccentricity_squared(f)
    return {'a':a,'b':b,'f':f,'e2':e2,'e':np.sqrt(e2),'ep2':ep2,'ep':np.sqrt(ep2),'linear_eccentricity':linear_eccentricity(a,f)}
