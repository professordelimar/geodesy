"""Transformações de coordenadas geodésicas, ECEF e ENU."""
import numpy as np
from .constants import WGS84, DEG2RAD, RAD2DEG
from .ellipsoid import first_eccentricity_squared
def geodetic_to_ecef(lat,lon,h,a=WGS84.a,f=WGS84.f,degrees=True):
    lat=np.asarray(lat,dtype=float); lon=np.asarray(lon,dtype=float); h=np.asarray(h,dtype=float)
    if degrees: lat=lat*DEG2RAD; lon=lon*DEG2RAD
    e2=first_eccentricity_squared(f); N=a/np.sqrt(1.0-e2*np.sin(lat)**2)
    X=(N+h)*np.cos(lat)*np.cos(lon); Y=(N+h)*np.cos(lat)*np.sin(lon); Z=(N*(1.0-e2)+h)*np.sin(lat)
    return X,Y,Z
def ecef_to_geodetic(X,Y,Z,a=WGS84.a,f=WGS84.f,degrees=True,tol=1e-12,max_iter=20):
    X=np.asarray(X,dtype=float); Y=np.asarray(Y,dtype=float); Z=np.asarray(Z,dtype=float)
    e2=first_eccentricity_squared(f); lon=np.arctan2(Y,X); p=np.sqrt(X**2+Y**2); lat=np.arctan2(Z,p*(1.0-e2)); h=np.zeros_like(lat,dtype=float)
    for _ in range(max_iter):
        old=lat.copy(); N=a/np.sqrt(1.0-e2*np.sin(lat)**2); h=p/np.cos(lat)-N; lat=np.arctan2(Z,p*(1.0-e2*N/(N+h)))
        if np.nanmax(np.abs(lat-old))<tol: break
    N=a/np.sqrt(1.0-e2*np.sin(lat)**2); h=p/np.cos(lat)-N
    return (lat*RAD2DEG,lon*RAD2DEG,h) if degrees else (lat,lon,h)
def rotation_ecef_to_enu(lat0,lon0,degrees=True):
    if degrees: lat0=lat0*DEG2RAD; lon0=lon0*DEG2RAD
    slat,clat=np.sin(lat0),np.cos(lat0); slon,clon=np.sin(lon0),np.cos(lon0)
    return np.array([[-slon,clon,0.0],[-slat*clon,-slat*slon,clat],[clat*clon,clat*slon,slat]])
def ecef_to_enu(X,Y,Z,lat0,lon0,h0=0.0,a=WGS84.a,f=WGS84.f,degrees=True):
    X0,Y0,Z0=geodetic_to_ecef(lat0,lon0,h0,a,f,degrees); R=rotation_ecef_to_enu(lat0,lon0,degrees)
    d=np.vstack([np.asarray(X,dtype=float)-X0,np.asarray(Y,dtype=float)-Y0,np.asarray(Z,dtype=float)-Z0]); enu=R@d
    return enu[0],enu[1],enu[2]
def enu_to_ecef(E,N,U,lat0,lon0,h0=0.0,a=WGS84.a,f=WGS84.f,degrees=True):
    X0,Y0,Z0=geodetic_to_ecef(lat0,lon0,h0,a,f,degrees); R=rotation_ecef_to_enu(lat0,lon0,degrees)
    dxyz=R.T@np.vstack([np.asarray(E,dtype=float),np.asarray(N,dtype=float),np.asarray(U,dtype=float)])
    return X0+dxyz[0],Y0+dxyz[1],Z0+dxyz[2]
def geodetic_to_enu(lat,lon,h,lat0,lon0,h0=0.0,a=WGS84.a,f=WGS84.f,degrees=True):
    X,Y,Z=geodetic_to_ecef(lat,lon,h,a,f,degrees); return ecef_to_enu(X,Y,Z,lat0,lon0,h0,a,f,degrees)
def enu_to_geodetic(E,N,U,lat0,lon0,h0=0.0,a=WGS84.a,f=WGS84.f,degrees=True):
    X,Y,Z=enu_to_ecef(E,N,U,lat0,lon0,h0,a,f,degrees); return ecef_to_geodetic(X,Y,Z,a,f,degrees)
def geodetic_to_spherical(lat,lon,h,a=WGS84.a,f=WGS84.f,degrees=True):
    X,Y,Z=geodetic_to_ecef(lat,lon,h,a,f,degrees); r=np.sqrt(X**2+Y**2+Z**2); la=np.arcsin(Z/r); lo=np.arctan2(Y,X)
    return (r,la*RAD2DEG,lo*RAD2DEG) if degrees else (r,la,lo)
def degrees_to_dms(deg):
    sign=np.sign(deg); da=np.abs(deg); d=np.floor(da); mf=(da-d)*60; m=np.floor(mf); s=(mf-m)*60; return sign*d,m,s
def dms_to_degrees(d,m,s):
    sign=np.sign(d) if d!=0 else 1.0; return sign*(abs(d)+m/60.0+s/3600.0)
def haversine_distance(lat1,lon1,lat2,lon2,radius=6371008.8,degrees=True):
    if degrees: lat1,lon1,lat2,lon2=np.array([lat1,lon1,lat2,lon2])*DEG2RAD
    dlat=lat2-lat1; dlon=lon2-lon1; a=np.sin(dlat/2)**2+np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    return radius*2*np.arctan2(np.sqrt(a),np.sqrt(1-a))
