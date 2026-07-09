# ADC-Repository
CubeSat ADCS subsystem developed at SLU SSRL for NASA CSLI and UNP missions. Implements TRIAD-based attitude determination, B-dot and dipole control laws, and magnetorquer hardware testing on embedded systems. Includes full software-hardware integration, IMU-based estimation, and experimental validation of torque generation.

# Structure

**Attitude Determination** contains the TRIAD attitude solution pipeline, split into the three pieces that feed it. Sun Vectors builds the measured and reference sun vectors. Magnetic Field Vectors builds the measured and reference magnetic field vectors, with the reference side using an IGRF model. TRIAD combines both vector pairs into a single attitude solution and includes the SPICE kernels (`naif0012.tls`, `pck00010.tpc`, `de440s.bsp`) needed for the frame transformations and ephemeris lookups behind the reference vectors.

**Bdot Controller** contains the detumbling control law and its hardware interface, including the IMU driver, magnetometer calibration, and the Raspberry Pi to H-bridge to magnetorquer chain. It also contains the STK and MATLAB simulation used to validate the B-dot controller against a modeled DARLA-02 spin state before it ran on hardware.

**MTQ Testing** holds the experimental validation work, including the suspended test article campaign, a test that expects the box to spin along the string that holds the box, and a compass needle test, that expects the compass needle to deflect as soon as the MTQ is powered. Both tests are used to prove torque generation.

## Requirements traceability

Every script in this repository maps back to a specific requirement in `Docs/ADC_Requirements_Verification_Matrix.xlsx`. See that file for the full requirement text, verification method, and current status for each.

## Note on SPICE kernels

`naif0012.tls` and `pck00010.tpc` are small and included directly. `de440s.bsp` is a standard public JPL planetary ephemeris kernel, not something specific to this mission, and is included here at roughly 31 MB for convenience. If a future kernel update pushes the file size higher, it is better linked from NAIF's public kernel archive than committed directly.
