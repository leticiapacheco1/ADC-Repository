from sun_body import compute_sun_vector_body

sun_body = compute_sun_vector_body(
    1.0,
    0.0,
    0.0,
    0.0,
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0]
)

print(sun_body)