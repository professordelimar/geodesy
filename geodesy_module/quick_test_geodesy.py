from geodesy.coordinates import geodetic_to_ecef, ecef_to_geodetic
from geodesy.normal_gravity import normal_gravity_somigliana
from geodesy.heights import ellipsoidal_to_orthometric_height
from geodesy.geoid import bruns_formula
lat, lon, h = -3.2, -52.2, 100.0
X,Y,Z = geodetic_to_ecef(lat, lon, h)
lat2,lon2,h2 = ecef_to_geodetic(X,Y,Z)
print('Entrada:', lat, lon, h)
print('ECEF:', X, Y, Z)
print('Volta:', lat2, lon2, h2)
gamma = normal_gravity_somigliana(lat)
N = 12.5
print('Gravidade normal:', gamma)
print('Altura ortométrica:', ellipsoidal_to_orthometric_height(h,N))
print('Geoide por Bruns:', bruns_formula(N*gamma,gamma))
