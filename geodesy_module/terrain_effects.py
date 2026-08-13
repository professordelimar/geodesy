"""Correções e efeitos topográficos simples."""
import numpy as np
from .constants import G_CONSTANT,SI_TO_MGAL
def simple_bouguer_plate(height,density=2670.0): return 2*np.pi*G_CONSTANT*density*np.asarray(height,dtype=float)*SI_TO_MGAL
def residual_terrain_model(topography,reference_topography): return np.asarray(topography,dtype=float)-np.asarray(reference_topography,dtype=float)
def rtm_gravity_effect_simple(topography,reference_topography,density=2670.0): return simple_bouguer_plate(residual_terrain_model(topography,reference_topography),density)
def point_mass_topographic_effect(x_obs,y_obs,z_obs,x_mass,y_mass,z_mass,mass):
    dx=np.asarray(x_obs,dtype=float)-x_mass; dy=np.asarray(y_obs,dtype=float)-y_mass; dz=np.asarray(z_obs,dtype=float)-z_mass; r=np.sqrt(dx*dx+dy*dy+dz*dz); return (-G_CONSTANT*mass*dz/r**3)*SI_TO_MGAL
def prism_mass(dx,dy,dz,density=2670.0): return dx*dy*dz*density
def grid_topographic_point_mass_effect(X,Y,H,density=2670.0,cell_area=1.0,observation_height=0.0):
    H=np.asarray(H,dtype=float); mass=H*cell_area*density; dz=observation_height-H/2; r=np.sqrt(cell_area/np.pi+dz**2); return (-G_CONSTANT*mass*dz/r**3)*SI_TO_MGAL
