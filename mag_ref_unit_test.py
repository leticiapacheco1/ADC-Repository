from mag_ref import load_spice_kernels
from mag_ref import compute_reference_magnetic_field_vector

load_spice_kernels()

mag_ref = compute_reference_magnetic_field_vector(
    "2026-05-29T12:00:00",
    38.6270,
    -90.1994,
    400000
)

print(mag_ref)