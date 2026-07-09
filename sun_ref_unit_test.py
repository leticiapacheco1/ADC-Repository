from sun_ref import load_spice_kernels
from sun_ref import compute_reference_sun_vector

print("1")
load_spice_kernels()

print("2")
sun_ref = compute_reference_sun_vector(
    "2026-05-29T12:00:00",
    38.6270,
    -90.1994,
    400000
)

print("3")
print(sun_ref)