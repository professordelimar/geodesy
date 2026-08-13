"""Campo normal da gravidade do elipsoide de referência."""
import numpy as np
from .constants import WGS84,GAMMA_E_WGS84,K_SOMIGLIANA_WGS84,DEG2RAD,SI_TO_MGAL
from .ellipsoid import first_eccentricity_squared
def normal_gravity_somigliana(lat,a=WGS84.a,f=WGS84.f,gamma_e=GAMMA_E_WGS84,k=K_SOMIGLIANA_WGS84,degrees=True):
    lat=np.asarray(lat,dtype=float); lat=lat*DEG2RAD if degrees else lat; e2=first_eccentricity_squared(f); s2=np.sin(lat)**2
    return gamma_e*(1+k*s2)/np.sqrt(1-e2*s2)
def free_air_gradient_si(): return 3.086e-6
def free_air_gradient_mgal_per_m(): return free_air_gradient_si()*SI_TO_MGAL
def normal_gravity_height(lat,h,degrees=True): return normal_gravity_somigliana(lat,degrees=degrees)-free_air_gradient_si()*np.asarray(h,dtype=float)
def centrifugal_acceleration(lat,h=0.0,a=WGS84.a,omega=WGS84.omega,degrees=True):
    lat=np.asarray(lat,dtype=float); lat=lat*DEG2RAD if degrees else lat; return omega**2*(a+np.asarray(h,dtype=float))*np.cos(lat)
def centrifugal_potential(lat,h=0.0,a=WGS84.a,omega=WGS84.omega,degrees=True):
    lat=np.asarray(lat,dtype=float); lat=lat*DEG2RAD if degrees else lat; rp=(a+np.asarray(h,dtype=float))*np.cos(lat); return 0.5*omega**2*rp**2
def normal_potential_spherical(r,GM=WGS84.GM): return GM/np.asarray(r,dtype=float)
