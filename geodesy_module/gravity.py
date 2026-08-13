"""Funções básicas para potencial e atração gravitacional."""
import numpy as np
from .constants import G_CONSTANT,WGS84
def newtonian_potential_point_mass(x,y,z,mass,x0=0,y0=0,z0=0,G=G_CONSTANT):
    dx=np.asarray(x,dtype=float)-x0; dy=np.asarray(y,dtype=float)-y0; dz=np.asarray(z,dtype=float)-z0; r=np.sqrt(dx*dx+dy*dy+dz*dz); return G*mass/r
def newtonian_attraction_point_mass(x,y,z,mass,x0=0,y0=0,z0=0,G=G_CONSTANT):
    dx=np.asarray(x,dtype=float)-x0; dy=np.asarray(y,dtype=float)-y0; dz=np.asarray(z,dtype=float)-z0; r=np.sqrt(dx*dx+dy*dy+dz*dz); fac=-G*mass/r**3; return fac*dx,fac*dy,fac*dz
def gravity_spherical_earth(r,GM=WGS84.GM): return GM/np.asarray(r,dtype=float)**2
def gravity_potential_spherical(r,GM=WGS84.GM): return GM/np.asarray(r,dtype=float)
def disturbing_potential_from_geoid(N,gamma): return np.asarray(N,dtype=float)*np.asarray(gamma,dtype=float)
def radial_gravity_gradient_spherical(r,GM=WGS84.GM): return -2*GM/np.asarray(r,dtype=float)**3
