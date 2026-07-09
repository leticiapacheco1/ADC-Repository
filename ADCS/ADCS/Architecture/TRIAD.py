import numpy as np


def normalize(vector):
    norm = np.linalg.norm(vector)

    if norm < 1e-6:
        return np.zeros(3)

    return vector / norm


def compute_body_triad_matrix(sun_body, mag_body):
    t1_body = normalize(sun_body)
    mag_body = normalize(mag_body)

    t2_body = normalize(
        np.cross(t1_body, mag_body)
    )

    t3_body = np.cross(
        t1_body,
        t2_body
    )

    return np.column_stack((
        t1_body,
        t2_body,
        t3_body
    ))


def compute_inertial_triad_matrix(sun_ref, mag_ref):
    t1_ref = normalize(sun_ref)
    mag_ref = normalize(mag_ref)

    t2_ref = normalize(
        np.cross(t1_ref, mag_ref)
    )

    t3_ref = np.cross(
        t1_ref,
        t2_ref
    )

    return np.column_stack((
        t1_ref,
        t2_ref,
        t3_ref
    ))


def compute_euler_angles(C_bi):
    roll = np.arctan2(C_bi[2, 1], C_bi[2, 2])
    pitch = np.arcsin(-C_bi[2, 0])
    yaw = np.arctan2(C_bi[1, 0], C_bi[0, 0])

    return (
        np.degrees(roll),
        np.degrees(pitch),
        np.degrees(yaw)
    )


def compute_attitude_euler_angles(
    sun_body,
    mag_body,
    sun_ref,
    mag_ref
):
    Tb = compute_body_triad_matrix(sun_body, mag_body)
    Ti = compute_inertial_triad_matrix(sun_ref, mag_ref)

    C_bi = Tb @ Ti.T

    return compute_euler_angles(C_bi)

