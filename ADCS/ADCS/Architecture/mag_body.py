import numpy as np


def normalize(vector):
    norm = np.linalg.norm(vector)

    if norm < 1e-6:
        return np.zeros(3)

    return vector / norm


def compute_measured_magnetic_field_vector(
    magnetic_strength_uT_x,
    magnetic_strength_uT_y,
    magnetic_strength_uT_z
):
    magnetic_vector_body = np.array([
        magnetic_strength_uT_x,
        magnetic_strength_uT_y,
        magnetic_strength_uT_z
    ])

    return normalize(magnetic_vector_body)
