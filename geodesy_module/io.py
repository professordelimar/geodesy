"""Leitura, escrita e organização de dados geodésicos em arquivos texto."""
import numpy as np
def read_xyz_file(filename,skiprows=0,delimiter=None):
    data=np.loadtxt(filename,skiprows=skiprows,delimiter=delimiter)
    if data.shape[1]<3: raise ValueError('O arquivo deve conter pelo menos três colunas.')
    return data[:,0],data[:,1],data[:,2]
def write_xyz_file(filename,x,y,z,header='x y z',fmt='%.10f'):
    np.savetxt(filename,np.column_stack([x,y,z]),fmt=fmt,header=header,comments='')
def xyz_to_grid(x,y,z):
    x=np.asarray(x,dtype=float); y=np.asarray(y,dtype=float); z=np.asarray(z,dtype=float); xu=np.unique(x); yu=np.unique(y); X,Y=np.meshgrid(xu,yu); Z=np.full_like(X,np.nan,dtype=float); xi={v:i for i,v in enumerate(xu)}; yi={v:i for i,v in enumerate(yu)}
    for xx,yy,zz in zip(x,y,z): Z[yi[yy],xi[xx]]=zz
    return X,Y,Z
def grid_to_xyz(X,Y,Z): return X.ravel(),Y.ravel(),Z.ravel()
def read_lon_lat_value(filename,skiprows=0,delimiter=None): return read_xyz_file(filename,skiprows,delimiter)
def save_grid_xyz(filename,X,Y,Z,header='lon lat value',fmt='%.10f'):
    x,y,z=grid_to_xyz(X,Y,Z); write_xyz_file(filename,x,y,z,header,fmt)
def check_regular_grid(x,y):
    x=np.asarray(x); y=np.asarray(y); nx=len(np.unique(x)); ny=len(np.unique(y)); return len(x)==nx*ny,nx,ny
