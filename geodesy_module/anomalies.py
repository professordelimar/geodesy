"""
anomalies.py
Grandezas anômalas do campo de gravidade.

As funções trabalham preferencialmente em mGal para anomalias,
exceto quando indicado.

Inclui:
- gravidade teórica/normal por latitude;
- correção de latitude;
- gravidade corrigida pela latitude;
- distúrbio da gravidade;
- correção de ar livre;
- anomalia de ar livre;
- correção de Bouguer simples;
- anomalia de Bouguer;
- anomalia de altura;
- ondulação geoidal;
- deflexão da vertical.
"""

import numpy as np

from .constants import SI_TO_MGAL, MGAL_TO_SI
from .normal_gravity import (
    normal_gravity_somigliana,
    normal_gravity_height,
    free_air_gradient_mgal_per_m
)


def theoretical_gravity(latitude, degrees=True):
    """
    Gravidade teórica/normal em função da latitude.

    Parâmetros
    ----------
    latitude : float ou array
        Latitude geodésica.
    degrees : bool
        Se True, latitude em graus. Se False, em radianos.

    Retorna
    -------
    gamma_mgal : float ou array
        Gravidade normal em mGal.

    Observação
    ----------
    Esta função é a forma explícita da gravidade teórica usada
    na correção de latitude das anomalias gravimétricas.
    """
    return normal_gravity_somigliana(latitude, degrees=degrees) * SI_TO_MGAL


def latitude_correction(latitude, reference_lat=0.0, degrees=True):
    """
    Correção de latitude em mGal.

    Calcula a diferença entre a gravidade normal na latitude informada
    e a gravidade normal em uma latitude de referência.

    CL = gamma(latitude) - gamma(reference_lat)

    Por padrão, usa o Equador como referência.

    Parâmetros
    ----------
    latitude : float ou array
        Latitude geodésica dos pontos.
    reference_lat : float
        Latitude de referência.
    degrees : bool
        Se True, latitude em graus. Se False, em radianos.

    Retorna
    -------
    CL : float ou array
        Correção de latitude em mGal.
    """
    gamma_lat = theoretical_gravity(latitude, degrees=degrees)
    gamma_ref = theoretical_gravity(reference_lat, degrees=degrees)

    return gamma_lat - gamma_ref


def latitude_corrected_gravity(g_observed_mgal, latitude, reference_lat=0.0, degrees=True):
    """
    Gravidade observada corrigida da variação normal com a latitude.

    Uma forma simples de remover o efeito latitudinal é subtrair a
    diferença entre a gravidade normal local e a gravidade normal em
    uma latitude de referência:

    g_corr = g_obs - [gamma(latitude) - gamma(reference_lat)]

    Parâmetros
    ----------
    g_observed_mgal : float ou array
        Gravidade observada em mGal.
    latitude : float ou array
        Latitude geodésica dos pontos.
    reference_lat : float
        Latitude de referência.
    degrees : bool
        Se True, latitude em graus.

    Retorna
    -------
    g_corr : float ou array
        Gravidade corrigida da latitude em mGal.
    """
    return np.asarray(g_observed_mgal, dtype=float) - latitude_correction(
        latitude,
        reference_lat=reference_lat,
        degrees=degrees
    )


def gravity_disturbance(g_observed, gamma_at_point, input_si=False, output_mgal=True):
    """
    Distúrbio da gravidade: delta_g = g_P - gamma_P.

    Parâmetros
    ----------
    g_observed : float ou array
        Gravidade observada.
    gamma_at_point : float ou array
        Gravidade normal no ponto.
    input_si : bool
        Se True, entradas em m/s². Caso contrário, assume mGal.
    output_mgal : bool
        Se input_si=True e output_mgal=True, converte saída para mGal.

    Retorna
    -------
    delta_g : float ou array
        Distúrbio da gravidade.
    """
    dg = np.asarray(g_observed, dtype=float) - np.asarray(gamma_at_point, dtype=float)

    if input_si and output_mgal:
        return dg * SI_TO_MGAL

    return dg


def free_air_correction(height_m):
    """
    Correção de ar livre aproximada em mGal.

    FAC = 0.3086 h

    Parâmetros
    ----------
    height_m : float ou array
        Altura em metros.

    Retorna
    -------
    FAC : float ou array
        Correção de ar livre em mGal.
    """
    return free_air_gradient_mgal_per_m() * np.asarray(height_m, dtype=float)


def free_air_anomaly(g_observed_mgal, lat, height_m, gamma0_mgal=None, degrees=True):
    """
    Anomalia de ar livre aproximada:

    Δg_FA = g_obs + FAC - gamma0

    onde:
        FAC = 0.3086 h
        gamma0 = gravidade normal na latitude

    Parâmetros
    ----------
    g_observed_mgal : float ou array
        Gravidade observada absoluta em mGal.
    lat : float ou array
        Latitude geodésica.
    height_m : float ou array
        Altura em metros.
    gamma0_mgal : float ou array ou None
        Gravidade normal em mGal. Se None, calcula por Somigliana.
    degrees : bool
        Se True, latitude em graus.

    Retorna
    -------
    FA : float ou array
        Anomalia de ar livre em mGal.

    Atenção
    -------
    g_observed_mgal deve representar gravidade observada absoluta
    em mGal, geralmente da ordem de 978000 mGal. Não passe aqui
    uma anomalia já reduzida, como distúrbio ou anomalia local.
    """
    if gamma0_mgal is None:
        gamma0_mgal = theoretical_gravity(lat, degrees=degrees)

    return np.asarray(g_observed_mgal, dtype=float) + free_air_correction(height_m) - gamma0_mgal


def simple_bouguer_correction(height_m, density=2670.0, G=6.67430e-11):
    """
    Correção de Bouguer simples para placa infinita:

    BC = 2*pi*G*rho*h

    A saída é dada em mGal.

    Parâmetros
    ----------
    height_m : float ou array
        Altura/topografia em metros.
    density : float
        Densidade em kg/m³.
    G : float
        Constante gravitacional universal.

    Retorna
    -------
    BC : float ou array
        Correção de Bouguer simples em mGal.
    """
    return 2.0 * np.pi * G * density * np.asarray(height_m, dtype=float) * SI_TO_MGAL


def bouguer_anomaly(g_observed_mgal, lat, height_m, density=2670.0, gamma0_mgal=None, degrees=True):
    """
    Anomalia de Bouguer simples:

    Δg_B = g_obs + FAC - BC - gamma0

    ou:

    Δg_B = Δg_FA - BC

    Parâmetros
    ----------
    g_observed_mgal : float ou array
        Gravidade observada absoluta em mGal.
    lat : float ou array
        Latitude geodésica.
    height_m : float ou array
        Altura/topografia em metros.
    density : float
        Densidade em kg/m³.
    gamma0_mgal : float ou array ou None
        Gravidade normal em mGal. Se None, calcula por Somigliana.
    degrees : bool
        Se True, latitude em graus.

    Retorna
    -------
    BA : float ou array
        Anomalia de Bouguer em mGal.
    """
    return (
        free_air_anomaly(
            g_observed_mgal,
            lat,
            height_m,
            gamma0_mgal=gamma0_mgal,
            degrees=degrees
        )
        - simple_bouguer_correction(height_m, density=density)
    )


def bouguer_anomaly_from_free_air(free_air_anomaly_mgal, height_m, density=2670.0):
    """
    Calcula a anomalia de Bouguer a partir da anomalia de ar livre:

    Δg_B = Δg_FA - BC

    Parâmetros
    ----------
    free_air_anomaly_mgal : float ou array
        Anomalia de ar livre em mGal.
    height_m : float ou array
        Altura/topografia em metros.
    density : float
        Densidade em kg/m³.

    Retorna
    -------
    BA : float ou array
        Anomalia de Bouguer em mGal.
    """
    return np.asarray(free_air_anomaly_mgal, dtype=float) - simple_bouguer_correction(
        height_m,
        density=density
    )


def height_anomaly_from_disturbing_potential(T, gamma):
    """
    Anomalia de altura:

    zeta = T/gamma

    Parâmetros
    ----------
    T : float ou array
        Potencial perturbador em m²/s².
    gamma : float ou array
        Gravidade normal em m/s².

    Retorna
    -------
    zeta : float ou array
        Anomalia de altura em metros.
    """
    return np.asarray(T, dtype=float) / np.asarray(gamma, dtype=float)


def geoid_height_from_disturbing_potential(T, gamma):
    """
    Ondulação geoidal aproximada:

    N = T/gamma

    Parâmetros
    ----------
    T : float ou array
        Potencial perturbador em m²/s².
    gamma : float ou array
        Gravidade normal em m/s².

    Retorna
    -------
    N : float ou array
        Ondulação geoidal em metros.
    """
    return np.asarray(T, dtype=float) / np.asarray(gamma, dtype=float)


def deflection_of_vertical_components(dT_dnorth, dT_deast, gamma):
    """
    Componentes aproximadas da deflexão da vertical:

    xi  ≈ -1/gamma * dT/dnorth
    eta ≈ -1/gamma * dT/deast

    Parâmetros
    ----------
    dT_dnorth : float ou array
        Derivada do potencial perturbador na direção norte.
    dT_deast : float ou array
        Derivada do potencial perturbador na direção leste.
    gamma : float ou array
        Gravidade normal em m/s².

    Retorna
    -------
    xi, eta : float ou array
        Componentes da deflexão da vertical em radianos.
    """
    gamma = np.asarray(gamma, dtype=float)

    xi = -np.asarray(dT_dnorth, dtype=float) / gamma
    eta = -np.asarray(dT_deast, dtype=float) / gamma

    return xi, eta
