"""Ajustamento de redes geodésicas por mínimos quadrados."""
import numpy as np
def normal_equations(A,l,P=None):
    A=np.asarray(A,dtype=float); l=np.asarray(l,dtype=float); P=np.eye(A.shape[0]) if P is None else np.asarray(P,dtype=float); return A.T@P@A, A.T@P@l
def least_squares_adjustment(A,l,P=None):
    N,u=normal_equations(A,l,P); x=np.linalg.solve(N,u); v=np.asarray(A,dtype=float)@x-np.asarray(l,dtype=float); return x,v,N
def covariance_of_parameters(N,sigma0_squared=1.0): return sigma0_squared*np.linalg.inv(N)
def posterior_variance_factor(v,P=None,dof=None):
    v=np.asarray(v,dtype=float); P=np.eye(len(v)) if P is None else np.asarray(P,dtype=float); dof=len(v) if dof is None else dof; return (v.T@P@v)/dof
def error_ellipse_2d(Cxx):
    vals,vecs=np.linalg.eigh(np.asarray(Cxx,dtype=float)[:2,:2]); order=np.argsort(vals)[::-1]; vals=vals[order]; vecs=vecs[:,order]; return np.sqrt(vals[0]),np.sqrt(vals[1]),np.arctan2(vecs[1,0],vecs[0,0])
