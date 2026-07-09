import numpy as np
from TRIAD import compute_attitude_euler_angles

sun_body = np.array([1, 0, 0])
mag_body = np.array([0, 1, 0])

sun_ref = np.array([1, 0, 0])
mag_ref = np.array([0, 1, 0])

roll_deg, pitch_deg, yaw_deg = compute_attitude_euler_angles(
    sun_body,
    mag_body,
    sun_ref,
    mag_ref
)

print("Roll:", roll_deg)
print("Pitch:", pitch_deg)
print("Yaw:", yaw_deg)