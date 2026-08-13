"""Funções básicas para geoide e quasigeoide."""
import numpy as np
def bruns_formula(T,gamma): return np.asarray(T,dtype=float)/np.asarray(gamma,dtype=float)
def disturbing_potential_from_geoid_height(N,gamma): return np.asarray(N,dtype=float)*np.asarray(gamma,dtype=float)
def geoid_height_from_T(T,gamma): return bruns_formula(T,gamma)
def height_anomaly_from_T(T,gamma): return bruns_formula(T,gamma)
def geoid_to_quasigeoid(N,separation): return np.asarray(N,dtype=float)-np.asarray(separation,dtype=float)
def quasigeoid_to_geoid(zeta,separation): return np.asarray(zeta,dtype=float)+np.asarray(separation,dtype=float)
def ellipsoidal_to_orthometric(h,N): return np.asarray(h,dtype=float)-np.asarray(N,dtype=float)
def orthometric_to_ellipsoidal(H,N): return np.asarray(H,dtype=float)+np.asarray(N,dtype=float)
def geoid_residual(N_observed,N_model): return np.asarray(N_observed,dtype=float)-np.asarray(N_model,dtype=float)
def mean_dynamic_topography(sea_surface_height,geoid_height): return np.asarray(sea_surface_height,dtype=float)-np.asarray(geoid_height,dtype=float)
