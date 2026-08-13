"""Sistemas de alturas em Geodésia Física."""
import numpy as np
def geopotential_number(W0,WP): return np.asarray(W0,dtype=float)-np.asarray(WP,dtype=float)
def dynamic_height(C,gamma0=9.80665): return np.asarray(C,dtype=float)/gamma0
def orthometric_height_approx(C,mean_gravity=9.80665): return np.asarray(C,dtype=float)/mean_gravity
def normal_height(C,mean_normal_gravity=9.8062): return np.asarray(C,dtype=float)/mean_normal_gravity
def ellipsoidal_to_orthometric_height(h,N): return np.asarray(h,dtype=float)-np.asarray(N,dtype=float)
def orthometric_to_ellipsoidal_height(H,N): return np.asarray(H,dtype=float)+np.asarray(N,dtype=float)
def ellipsoidal_to_normal_height(h,zeta): return np.asarray(h,dtype=float)-np.asarray(zeta,dtype=float)
def normal_to_ellipsoidal_height(H_normal,zeta): return np.asarray(H_normal,dtype=float)+np.asarray(zeta,dtype=float)
def geoid_quasigeoid_separation(N,zeta): return np.asarray(N,dtype=float)-np.asarray(zeta,dtype=float)
