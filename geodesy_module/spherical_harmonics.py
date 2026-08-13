"""Ferramentas básicas para harmônicos esféricos."""
import numpy as np
from math import factorial
from .constants import WGS84
def legendre_polynomial(n,x):
    x=np.asarray(x,dtype=float)
    if n==0: return np.ones_like(x)
    if n==1: return x
    P0=np.ones_like(x); P1=x
    for k in range(2,n+1): P0,P1=P1,((2*k-1)*x*P1-(k-1)*P0)/k
    return P1
def associated_legendre(l,m,x):
    if m<0 or m>l: raise ValueError('Use 0 <= m <= l.')
    x=np.asarray(x,dtype=float); pmm=np.ones_like(x)
    if m>0:
        somx2=np.sqrt(np.maximum(0,1-x*x)); fact=1.0
        for _ in range(1,m+1): pmm*= -fact*somx2; fact+=2
    if l==m: return pmm
    pmmp1=x*(2*m+1)*pmm
    if l==m+1: return pmmp1
    p2,p1=pmm,pmmp1
    for ll in range(m+2,l+1): pll=((2*ll-1)*x*p1-(ll+m-1)*p2)/(ll-m); p2,p1=p1,pll
    return pll
def fully_normalized_legendre(l,m,x):
    delta=1.0 if m==0 else 0.0; norm=np.sqrt((2-delta)*(2*l+1)*factorial(l-m)/factorial(l+m)); return norm*associated_legendre(l,m,x)
def spherical_harmonic(l,m,lat,lon,C=1.0,S=0.0,degrees=True):
    lat=np.asarray(lat,dtype=float); lon=np.asarray(lon,dtype=float)
    if degrees: lat=np.deg2rad(lat); lon=np.deg2rad(lon)
    return fully_normalized_legendre(l,m,np.sin(lat))*(C*np.cos(m*lon)+S*np.sin(m*lon))
def potential_from_coefficients(lat,lon,r,C,S=None,GM=WGS84.GM,a=WGS84.a,lmax=None,degrees=True):
    C=np.asarray(C,dtype=float); S=np.zeros_like(C) if S is None else np.asarray(S,dtype=float); lmax=C.shape[0]-1 if lmax is None else lmax
    summ=np.zeros_like(np.asarray(lat,dtype=float),dtype=float)
    for l in range(lmax+1):
        inner=np.zeros_like(summ,dtype=float)
        for m in range(l+1): inner+=spherical_harmonic(l,m,lat,lon,C[l,m],S[l,m],degrees)
        inner=np.asarray(inner,dtype=float); summ+=(a/np.asarray(r,dtype=float))**l*inner
    return GM/np.asarray(r,dtype=float)*summ
def degree_variance(C,S=None,lmax=None):
    C=np.asarray(C,dtype=float); S=np.zeros_like(C) if S is None else np.asarray(S,dtype=float); lmax=C.shape[0]-1 if lmax is None else lmax
    return np.array([np.sum(C[l,:l+1]**2+S[l,:l+1]**2) for l in range(lmax+1)])
def truncate_coefficients(C,S=None,lmax=10):
    Ct=np.asarray(C,dtype=float)[:lmax+1,:lmax+1].copy();
    if S is None: return Ct,None
    return Ct,np.asarray(S,dtype=float)[:lmax+1,:lmax+1].copy()
