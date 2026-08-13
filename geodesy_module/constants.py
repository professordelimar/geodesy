"""Constantes geodésicas e parâmetros de elipsoides de referência."""
from dataclasses import dataclass
import numpy as np
G_CONSTANT = 6.67430e-11
SPEED_OF_LIGHT = 299792458.0
MEAN_EARTH_RADIUS = 6371008.8
EARTH_MASS = 5.9722e24
WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_GM = 3.986004418e14
WGS84_OMEGA = 7.292115e-5
GRS80_A = 6378137.0
GRS80_F = 1.0 / 298.257222101
GRS80_GM = 3.986005e14
GRS80_OMEGA = 7.292115e-5
GAMMA_E_WGS84 = 9.7803253359
GAMMA_P_WGS84 = 9.8321849378
K_SOMIGLIANA_WGS84 = 0.00193185265241
DEG2RAD = np.pi/180.0
RAD2DEG = 180.0/np.pi
MGAL_TO_SI = 1e-5
SI_TO_MGAL = 1e5
@dataclass
class Ellipsoid:
    name: str; a: float; f: float; GM: float=WGS84_GM; omega: float=WGS84_OMEGA
    @property
    def b(self): return self.a*(1.0-self.f)
    @property
    def e2(self): return 2.0*self.f-self.f**2
    @property
    def ep2(self): return self.e2/(1.0-self.e2)
WGS84=Ellipsoid('WGS84',WGS84_A,WGS84_F,WGS84_GM,WGS84_OMEGA)
GRS80=Ellipsoid('GRS80',GRS80_A,GRS80_F,GRS80_GM,GRS80_OMEGA)
def get_ellipsoid(name='WGS84'):
    name=name.upper()
    if name=='WGS84': return WGS84
    if name=='GRS80': return GRS80
    raise ValueError("Use 'WGS84' ou 'GRS80'.")
