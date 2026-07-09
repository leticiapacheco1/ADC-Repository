import numpy as np
import spiceypy as spice
import ppigrf
from datetime import datetime


def load_spice_kernels():
    spice.furnsh("naif0012.tls")
    spice.furnsh("de440s.bsp")
    spice.furnsh("pck00010.tpc")


def normalize(vector):
    norm = np.linalg.norm(vector)

    if norm < 1e-6:
        return np.zeros(3)

    return vector / norm


def enu_to_ecef_vector(
    east,
    north,
    up,
    latitude_deg,
    longitude_deg
):
    lat = np.radians(latitude_deg)
    lon = np.radians(longitude_deg)

    east_unit = np.array([
        -np.sin(lon),
        np.cos(lon),
        0.0
    ])

    north_unit = np.array([
        -np.sin(lat) * np.cos(lon),
        -np.sin(lat) * np.sin(lon),
        np.cos(lat)
    ])

    up_unit = np.array([
        np.cos(lat) * np.cos(lon),
        np.cos(lat) * np.sin(lon),
        np.sin(lat)
    ])

    return (
        east * east_unit
        + north * north_unit
        + up * up_unit
    )


def ecef_to_j2000_vector(
    ecef_vector,
    utc_time
):
    et = spice.utc2et(utc_time)

    rotation_matrix = spice.pxform(
        "IAU_EARTH",
        "J2000",
        et
    )

    return rotation_matrix @ ecef_vector


def compute_reference_magnetic_field_vector(
    utc_time,
    latitude_deg,
    longitude_deg,
    altitude_m
):
    date = datetime.fromisoformat(utc_time)
    altitude_km = altitude_m / 1000.0

    Be_nT, Bn_nT, Bu_nT = ppigrf.igrf(
        longitude_deg,
        latitude_deg,
        altitude_km,
        date
    )

    Be_uT = Be_nT / 1000.0
    Bn_uT = Bn_nT / 1000.0
    Bu_uT = Bu_nT / 1000.0

    magnetic_ecef_uT = enu_to_ecef_vector(
        Be_uT,
        Bn_uT,
        Bu_uT,
        latitude_deg,
        longitude_deg
    )

    magnetic_j2000_uT = ecef_to_j2000_vector(
        magnetic_ecef_uT,
        utc_time
    )

    return normalize(magnetic_j2000_uT)
