import numpy as np
import spiceypy as spice


def normalize(vector):
    norm = np.linalg.norm(vector)

    if norm < 1e-6:
        return np.zeros(3)

    return vector / norm


def load_spice_kernels():
    spice.furnsh("naif0012.tls")
    spice.furnsh("de440s.bsp")
    spice.furnsh("pck00010.tpc")


def lla_to_ecef(
    latitude_deg,
    longitude_deg,
    altitude_m
):
    lat = np.radians(latitude_deg)
    lon = np.radians(longitude_deg)

    a = 6378137.0
    e2 = 6.69437999014e-3

    N = a / np.sqrt(
        1 - e2 * np.sin(lat)**2
    )

    x = (N + altitude_m) * np.cos(lat) * np.cos(lon)
    y = (N + altitude_m) * np.cos(lat) * np.sin(lon)
    z = ((1 - e2) * N + altitude_m) * np.sin(lat)

    return np.array([x, y, z]) / 1000.0


def ecef_to_j2000(
    ecef_position_km,
    utc_time
):
    et = spice.utc2et(utc_time)

    rotation_matrix = spice.pxform(
    "IAU_EARTH",
    "J2000",
    et
)

    return rotation_matrix @ ecef_position_km


def compute_reference_sun_vector(
    utc_time,
    latitude_deg,
    longitude_deg,
    altitude_m
):
    spacecraft_ecef_km = lla_to_ecef(
        latitude_deg,
        longitude_deg,
        altitude_m
    )

    spacecraft_j2000_km = ecef_to_j2000(
        spacecraft_ecef_km,
        utc_time
    )

    et = spice.utc2et(utc_time)

    sun_position_j2000_km, light_time = spice.spkpos(
        "SUN",
        et,
        "J2000",
        "LT+S",
        "EARTH"
    )

    spacecraft_to_sun_vector = (
        np.array(sun_position_j2000_km)
        - spacecraft_j2000_km
    )

    return normalize(spacecraft_to_sun_vector)
