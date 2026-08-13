import numpy as np
from .constants import MEAN_EARTH_RADIUS, GAMMA_E_WGS84, MGAL_TO_SI

def spherical_distance(lat1, lon1, lat2, lon2, degrees=True):
    lat1 = np.asarray(lat1, dtype=float)
    lon1 = np.asarray(lon1, dtype=float)
    lat2 = np.asarray(lat2, dtype=float)
    lon2 = np.asarray(lon2, dtype=float)
    if degrees:
        lat1 = np.deg2rad(lat1); lon1 = np.deg2rad(lon1)
        lat2 = np.deg2rad(lat2); lon2 = np.deg2rad(lon2)
    dlon = lon2 - lon1
    cospsi = np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(dlon)
    return np.arccos(np.clip(cospsi, -1.0, 1.0))

def stokes_kernel(psi):
    psi = np.asarray(psi, dtype=float)
    psi_safe = np.maximum(psi, 1e-12)
    s = np.sin(psi_safe/2.0)
    return 1.0/s - 6.0*s + 1.0 - 5.0*np.cos(psi_safe) - 3.0*np.cos(psi_safe)*np.log(s + s**2)

def stokes_integral_point(lat0, lon0, lat, lon, gravity_anomaly_mgal,
                          cell_area_spherical, R=MEAN_EARTH_RADIUS,
                          gamma=GAMMA_E_WGS84, degrees=True,
                          exclude_singularity=True, min_psi=1e-10):
    psi = spherical_distance(lat0, lon0, lat, lon, degrees=degrees)
    S = stokes_kernel(psi)
    dg_si = np.asarray(gravity_anomaly_mgal, dtype=float) * MGAL_TO_SI
    dOmega = np.asarray(cell_area_spherical, dtype=float)
    if exclude_singularity:
        mask = (psi > min_psi) & np.isfinite(S) & np.isfinite(dg_si) & np.isfinite(dOmega)
    else:
        mask = np.isfinite(S) & np.isfinite(dg_si) & np.isfinite(dOmega)
    return R/(4.0*np.pi*gamma) * np.nansum(dg_si[mask] * S[mask] * dOmega[mask])

def stokes_integral_grid(lat_grid, lon_grid, gravity_anomaly_mgal,
                         dlat_deg=None, dlon_deg=None,
                         R=MEAN_EARTH_RADIUS, gamma=GAMMA_E_WGS84,
                         exclude_singularity=True, min_psi=1e-10):
    lat_grid = np.asarray(lat_grid, dtype=float)
    lon_grid = np.asarray(lon_grid, dtype=float)
    dg = np.asarray(gravity_anomaly_mgal, dtype=float)
    if lat_grid.shape != lon_grid.shape or lat_grid.shape != dg.shape:
        raise ValueError("lat_grid, lon_grid e gravity_anomaly_mgal devem ter o mesmo shape.")
    if dlat_deg is None:
        dlat_deg = np.nanmedian(np.abs(np.diff(np.unique(lat_grid[:, 0]))))
    if dlon_deg is None:
        dlon_deg = np.nanmedian(np.abs(np.diff(np.unique(lon_grid[0, :]))))
    dOmega = np.cos(np.deg2rad(lat_grid)) * np.deg2rad(dlat_deg) * np.deg2rad(dlon_deg)
    N = np.zeros_like(dg, dtype=float)
    flat_lat = lat_grid.ravel(); flat_lon = lon_grid.ravel()
    flat_dg = dg.ravel(); flat_area = dOmega.ravel(); flat_N = N.ravel()
    for idx, (lt0, ln0) in enumerate(zip(flat_lat, flat_lon)):
        flat_N[idx] = stokes_integral_point(
            lt0, ln0, flat_lat, flat_lon, flat_dg, flat_area,
            R=R, gamma=gamma, degrees=True,
            exclude_singularity=exclude_singularity, min_psi=min_psi
        )
    return N

def planar_stokes_fft_approx(gravity_anomaly_mgal, dx, dy, gamma=GAMMA_E_WGS84):
    dg_si = np.asarray(gravity_anomaly_mgal, dtype=float) * MGAL_TO_SI
    if dg_si.ndim != 2:
        raise ValueError("gravity_anomaly_mgal deve ser uma matriz 2D.")
    ny, nx = dg_si.shape
    kx = 2.0*np.pi*np.fft.fftfreq(nx, d=dx)
    ky = 2.0*np.pi*np.fft.fftfreq(ny, d=dy)
    KX, KY = np.meshgrid(kx, ky)
    k = np.sqrt(KX**2 + KY**2)
    k[0, 0] = np.inf
    N = np.real(np.fft.ifft2(np.fft.fft2(dg_si)/(gamma*k)))
    N -= np.nanmean(N)
    return N
