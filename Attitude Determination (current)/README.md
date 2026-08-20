The extended Kalman filter (EKF) is the current attitude determination strategy for this ADCS subsystem, replacing the earlier TRIAD-based solution documented in attitude_determination_old. Where TRIAD solves for attitude as an instantaneous, memoryless geometric problem, the EKF treats attitude determination as a continuous estimation problem, propagating a state forward in time and correcting it with measurements as they become available rather than requiring a complete measurement set at every instant.

**Why an EKF**

TRIAD requires two valid, non-parallel vector measurements at the same instant to produce a solution. During eclipse, the sun sensor has nothing to measure, and TRIAD simply stops producing an output for that portion of every orbit. An EKF does not have this requirement. It carries a continuously propagated state forward using the spacecraft's rotational dynamics as a process model, and it only incorporates a measurement update, from the magnetometer, the sun sensor, or both, whenever that particular measurement happens to be available. When the sun sensor drops out during eclipse, the filter does not stall, it continues propagating from the gyro-driven dynamics and whatever the magnetometer, which has no eclipse dependency, is still providing.

**State**

The filter state is attitude, represented as either a quaternion or Modified Rodrigues Parameters depending on the implementation stage, augmented with gyro bias. Carrying gyro bias in the state rather than treating the gyro as a perfect rate measurement is what keeps the propagated attitude from drifting unbounded between measurement updates, since real gyros have a slowly varying bias that would otherwise integrate into significant attitude error over an orbit.

**Process model**

Between measurement updates, the state is propagated forward using the spacecraft's kinematic equations of motion driven by the gyro rate measurement, corrected for the current bias estimate. The covariance is propagated alongside the state using the linearized dynamics (the "extended" part of EKF, linearizing about the current state estimate at each step rather than assuming linear dynamics globally).

**Measurement model**

Two measurement types feed the update step, each compared against its corresponding reference vector the same way the TRIAD implementation did:

Magnetometer measurement, compared against the IGRF-derived reference magnetic field vector for the spacecraft's current position and time. Available continuously, with no eclipse dependency.
Sun sensor measurement (photodiode-derived body vector), compared against the SPICE-derived reference sun vector. Available only when the spacecraft is in sunlight.

When a measurement is unavailable, its update step is simply skipped for that cycle rather than substituted with a default or assumed value. The covariance grows somewhat faster during the gap, which correctly reflects the estimate's true uncertainty increasing without new information, but the filter continues producing a usable attitude estimate throughout.

**Relationship to the TRIAD implementation**

The EKF's measurement model still depends on the same body-frame and reference-frame vector construction documented in attitude_determination_old, the sun and magnetic field vectors, the SPICE and IGRF reference pipeline, and the underlying geometry are unchanged. TRIAD is what validated that pipeline in isolation, with a method simple enough that a bad output pointed directly at a sensor or geometry problem rather than hiding inside a filter's tuning. The EKF builds on top of that validated pipeline rather than replacing it.

**Known limitations, current development focus**

Two weak points from the earlier flight EKF implementation are the current focus of the redesign:

1. Propagation between measurement updates when no external torque model is included, i.e. the assumption that the only forces acting on the state during propagation are captured correctly by the process model.
2. The eclipse observability gap itself, ensuring the filter remains well-conditioned with an intermittent sun vector rather than implicitly assuming it is always present.
Status

Active development. This is the attitude determination strategy currently being flown forward, as opposed to TRIAD, which is retained as validated history rather than active flight code.
