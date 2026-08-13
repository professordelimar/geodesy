"""Transformações simples entre referenciais geodésicos."""
import numpy as np
def rotation_matrix_x(angle_rad):
    c,s=np.cos(angle_rad),np.sin(angle_rad); return np.array([[1,0,0],[0,c,-s],[0,s,c]],dtype=float)
def rotation_matrix_y(angle_rad):
    c,s=np.cos(angle_rad),np.sin(angle_rad); return np.array([[c,0,s],[0,1,0],[-s,0,c]],dtype=float)
def rotation_matrix_z(angle_rad):
    c,s=np.cos(angle_rad),np.sin(angle_rad); return np.array([[c,-s,0],[s,c,0],[0,0,1]],dtype=float)
def small_angle_rotation_matrix(rx,ry,rz): return np.array([[1,-rz,ry],[rz,1,-rx],[-ry,rx,1]],dtype=float)
def helmert_transform_7params(X,Y,Z,tx=0,ty=0,tz=0,rx=0,ry=0,rz=0,scale_ppm=0):
    P=np.vstack([np.asarray(X,dtype=float),np.asarray(Y,dtype=float),np.asarray(Z,dtype=float)]); R=small_angle_rotation_matrix(rx,ry,rz); s=1+scale_ppm*1e-6; T=np.array([[tx],[ty],[tz]],dtype=float); Q=T+s*(R@P); return Q[0],Q[1],Q[2]
def inverse_helmert_transform_7params(X,Y,Z,tx=0,ty=0,tz=0,rx=0,ry=0,rz=0,scale_ppm=0): return helmert_transform_7params(X,Y,Z,-tx,-ty,-tz,-rx,-ry,-rz,-scale_ppm)
def apply_rotation(X,Y,Z,R):
    Q=np.asarray(R)@np.vstack([np.asarray(X,dtype=float),np.asarray(Y,dtype=float),np.asarray(Z,dtype=float)]); return Q[0],Q[1],Q[2]
def arcsec_to_rad(value_arcsec): return value_arcsec*np.pi/(180.0*3600.0)
