TRIAD (Tri-Axial Attitude Determination) is used as a deterministic attitude estimation method that reconstructs spacecraft orientation by aligning measured body-frame vectors with known reference-frame vectors. In the ADCS subsystem, TRIAD is implemented as an initial coarse attitude solution to establish spacecraft orientation using magnetic field and sun vector measurements.

In this implementation, TRIAD is used as a first-principles estimator because it does not require temporal filtering or tuning parameters, making it suitable for initial attitude acquisition and debugging hardware behavior during ground testing. It provides a direct geometric relationship between measured vectors and known environmental references, which is useful when validating sensor consistency before applying dynamic filters such as Kalman-based estimation.

The implementation uses two non-parallel reference vectors: the Earth magnetic field vector obtained from magnetometer measurements and a secondary reference vector derived from either sun sensor photodiode readings or assumed inertial alignment during testing. These vectors are normalized and used to construct two orthogonal body-frame bases, which are then mapped to their corresponding reference-frame bases to compute the direction cosine matrix representing spacecraft attitude.

TRIAD serves as the initialization layer of the attitude determination pipeline. Once an initial orientation is computed, higher-frequency gyro integration and filtering methods are used to maintain temporal continuity. This layered approach allows TRIAD to provide geometric grounding while dynamic estimators handle time evolution of attitude.

In the context of this ADCS subsystem, TRIAD is primarily used for sensor validation, initial attitude acquisition, and debugging of magnetometer-sun sensor consistency rather than continuous flight estimation. It provides a reference framework for verifying that measured vectors behave consistently under controlled test conditions before transitioning to full control law execution.

**Status**

This module is no longer the active attitude determination strategy for flight. It is kept here as the earlier implementation in this subsystem's development, since it is what validated the sensor and reference vector pipeline that the current EKF implementation still depends on. See the Attitude Determination Current folder for the strategy actually in use now.

**Why This Was Moved Off Flight Status**

TRIAD's core limitation is structural, not a bug: it needs two valid, non-parallel vector measurements at the same instant to produce any attitude solution at all. In sunlight that pairing is the sun vector and the magnetic field vector, and the algorithm performs exactly as intended.

The problem is eclipse. During the portion of the orbit spent in Earth's shadow, the sun sensor has no signal to provide, and TRIAD has no second vector to construct a solution from. It does not degrade gracefully or return a lower confidence estimate; it simply produces nothing. For a CubeSat that spends a real, recurring fraction of every orbit in shadow, this is not an edge case, it is a coverage gap that shows up on a fixed schedule every single orbit.

A second, related limitation is that TRIAD is memoryless. Each solution is computed independently from whatever vectors are available at that instant, with no way to carry information forward from the previous good estimate. Combined with the eclipse gap, this means the subsystem would have no attitude estimate at all for a predictable, recurring portion of every orbit, with no mechanism to bridge that gap using prior knowledge of the spacecraft's motion.

Both limitations pointed toward the same fix: an estimator that propagates continuously using the spacecraft's own dynamics and only pulls in a measurement update when a given sensor happens to be available, rather than requiring a complete vector pair at every instant. That is the extended Kalman filter implementation in Attitude Determination Current.
