import numpy as np


def normalize(vector):
    norm = np.linalg.norm(vector)

    if norm < 1e-6:
        return np.zeros(3)

    return vector / norm


def compute_sun_vector_body(
    deployed_3u_panel_1_photodiode,
    deployed_3u_panel_2_photodiode,
    stationary_3u_panel_photodiode,
    stationary_1_5u_panel_photodiode,
    deployed_3u_panel_1_normal,
    deployed_3u_panel_2_normal,
    stationary_3u_panel_normal,
    stationary_1_5u_panel_normal
):
    sun_vector_body = (
        deployed_3u_panel_1_photodiode * normalize(deployed_3u_panel_1_normal)
        + deployed_3u_panel_2_photodiode * normalize(deployed_3u_panel_2_normal)
        + stationary_3u_panel_photodiode * normalize(stationary_3u_panel_normal)
        + stationary_1_5u_panel_photodiode * normalize(stationary_1_5u_panel_normal)
    )

    return normalize(sun_vector_body)
