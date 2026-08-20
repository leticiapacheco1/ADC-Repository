import numpy as np
from scipy.linalg import expm

def adcs_ekf(xhatIn, Pin, deltaT, sunline_i, magField_i, sensors, z):
    '''
    INPUTS:
    xhatIn : 1x7 vector representing the latest quaternion and angular velocity, angular velocities are in body coordinates
             [q_0, q_1, q_2, q_3, omega_x, omega_y, omega_z]
    Pin : 7x7 covariance estimate 
    deltaT : change in time from last calculated attitude
    sunline_i : the vector from the S/C to the sun (can be found with SPICE based on Earth position)
    magField_i : the magnetic field vector (can be found with IGRF model based on spacecraft position)
    sensors: 1x3 vector of 1s and 0s holding information for which sensors should be used in the update
             [useSolarPanelCurrents, useMagnetometer, useGyroscope]
    z : input of all the sensor readings, right now it is solar panel currents, magnetometer, and gyro
            I_x+ corresponds to the current value for the X+ solar panel, etc. for other panels
            magfield_x corresponds to the X component of the magnetic field vector, etc. for other values
            gyro_x corresponds to the roll rate around the X body axis in rad/s, etc. for other values
            [I_x+, I_x-, I_y+, I_y-, I_z-, magfield_x, magfield_y, magfield_z, gyro_x, gyro_y, gyro_z,]

    Sidenote: This algorithm uses the methodology of Babcock 2011 from the Small Satellite Conference, with some assistance from Lefferts 1982. 
    If something is confusing and unexplained in this code, then reading the corresponding section of Babcock 2011 may help.
    Babcock 2011 - https://digitalcommons.usu.edu/smallsat/2011/all2011/56/
    Lefferts 1982 - http://malcolmdshuster.com/Pub_1982c_J_LMS_scan.pdf
    '''

    '''
    CONSTANTS:
    Properties of the satellite that are necessary for the algorithm but will not change
    '''
    # Mass moment of inertia for the spacecraft, which is hard-coded in since it should not change on orbit
    # If someone in software thinks these values would be better placed elsewhere, be my guest
    # These values were also originally in units of lbm * mm^2 (vomit), so they have been converted to kg * m^2
    J = np.array([[15518.823,0,0],
         [0,15632.682,0],
         [0,0,4183.826]]) * (4.5359237 * 10**-7)
    Jinv = np.array([[0.0000644379,0,0],
            [0,0.0000639685,0],
            [0,0,0.000239016]]) / (4.5359237 * 10**-7)
    
    # Magnetic field vector for the spacecraft bar magnets in body coordinates, in units of A*m^2
    m_bar = np.array([0,0,0.817])

    # Maximum current for a single solar cell of the NanoAvionics solar panels
    # This number isn't actually important since we're turning the sun vector into a unit vector which removes this dependency,
    # but it doesn't hurt to keep it here right now
    I_max = 0.455 # amps
    # Minimum current for a single solar cell before it can be considered to be not in sunlight
    I_min = 0.010 # amps

    '''
    PROPAGATION:
    Updating the quaternion and angular rates based on the dynamics of the satellite
    '''
    # Process noise matrix
    Q = np.block([
        [np.eye(4)*0.1**2, np.zeros((4,3))],
        [np.zeros((3,4)), np.eye(3)*0.1**2]
    ])
    qHat = xhatIn[0:4]
    omegaHat = xhatIn[4:7]
    P = Pin

    # This is just breaking up our deltaT in order to get more accurate attitude propagation
    numChunks = 30
    dTinternal = deltaT / numChunks

    # Propagating a new attitude for each of our dTinternal
    for _ in range(numChunks):
        # Propagating q_hat
        qHat = qdot_spice(qHat, omegaHat, dTinternal)
        
        # tau = tau_build(qHat, m_bar, magField_i)
        tau = np.array([0,0,0])

        # This is a Runge-Kutta propagation which should be slightly more accurate (and prevent omega from blowing up)
        k1 = euler_dot(omegaHat, J, Jinv, tau)
        k2 = euler_dot(omegaHat + 0.5 * dTinternal * k1, J, Jinv, tau)
        k3 = euler_dot(omegaHat + 0.5 * dTinternal * k2, J, Jinv, tau)
        k4 = euler_dot(omegaHat + dTinternal * k3, J, Jinv, tau)

        omegaHat += (dTinternal / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

        # Updating the covariance
        P = P_update(P, qHat, omegaHat, dTinternal, J, Q)


    if np.max(sensors) > 0:
        '''
        QUATERNION TRUNCATION:
        A degree of freedom from the quaternion will be eliminated in order to ensure that all of the state variables remain independent and
        ensure convergence later on. 
        '''
        # This is the component of qHat that will be considered dependent and truncated
        # It is chosen as the element of q with the largest magnitude to improve the odds that the sign of the truncated component remains the same
        c = np.argmax(np.abs(qHat))
        sign = np.sign(qHat[c])

        # Truncating P by removing the appropriate row and column
        c_P = np.delete(P, (c), axis = 0)
        c_P = np.delete(c_P, (c), axis = 1)

        # Creating our truncated state vector
        c_xhat = np.concatenate((qHat[:c],qHat[c+1:],omegaHat))

        '''
        SENSOR MODELING
        Converting our propagated state and knowledge of sun and magnetic field vectors into an estimation of the observation vector, as well 
        as the H matrix. 
        '''
        # Direction Cosine Matrix to convert vectors in ECI to body coordinates
        DCM = np.array([[1-2*(qHat[2]**2+qHat[3]**2), 2*(qHat[1]*qHat[2]-qHat[0]*qHat[3]), 2*(qHat[1]*qHat[3]+qHat[0]*qHat[2])],
                [2*(qHat[1]*qHat[2]+qHat[0]*qHat[3]), 1-2*(qHat[1]**2+qHat[3]**2), 2*(qHat[2]*qHat[3]-qHat[0]*qHat[1])], 
                [2*(qHat[1]*qHat[3]-qHat[0]*qHat[2]), 2*(qHat[2]*qHat[3]+qHat[0]*qHat[1]), 1-2*(qHat[1]**2+qHat[2]**2)]])
        
        # DCM derivatives with respect to each quaternion
        dDCM_dq1 = np.array([
            [qHat[0], -qHat[3], qHat[2]], 
            [qHat[3], qHat[0], -qHat[1]],
            [-qHat[2], qHat[1], qHat[0]]         
            ])
        dDCM_dq2 = np.array([
            [qHat[1], qHat[2], qHat[3]], 
            [qHat[2], -qHat[1], -qHat[0]],
            [qHat[3], qHat[0], -qHat[1]]         
            ])    
        dDCM_dq3 = np.array([
            [-qHat[2], qHat[1], qHat[0]], 
            [qHat[1], qHat[2], qHat[3]],
            [-qHat[0], qHat[3], -qHat[2]]         
            ])    
        dDCM_dq4 = np.array([
            [-qHat[3], -qHat[0], qHat[1]], 
            [qHat[0], -qHat[3], qHat[2]],
            [qHat[1], qHat[2], qHat[3]]         
            ])
        
        # Initializing the observation matrix as 0x7 so that rows can be appended for each measurement
        H = np.zeros((0,7))

        # Initializing the measurement noise
        R = np.zeros((0,0))

        # Initializing the estimated measurement vector as 0x0 so that estimated measurements can be appended
        z_est = np.array([])

        # Initilizing the actual measurement vector as 0x0 so that the appropriate actual measurements can be appended
        z_act = np.array([])

        ## SUN VECTOR MODELING
        # Checks if we want to use the solar panel currents
        if sensors[0] == 1:
            # Adding to measurement noise matrix
            sun_var = 0.05**2
            R = np.block([
                [R, np.zeros((R.shape[0],3))],
                [np.zeros((3,R.shape[1])), np.eye(3)*sun_var]
            ])

            # Normal vectors for the faces with solar panels on them
            n = np.array([[1, 0, 0],
                          [-1, 0, 0],
                          [0, 1, 0],
                          [0, -1, 0],
                          [0, 0, -1]
                          ])

            # Getting the actual solar panel currents
            I_act = z[0:5]

            # Using the current measurements to construct a sun vector
            z_act_sun = np.zeros((3,))
            for i in range(n.shape[0]):
                z_act_sun = z_act_sun + (I_act[i] / I_max) * n[i,:]

            # Normalizing the estimated sun vector
            if np.linalg.norm(sunline_i) != 0:
                sunline_i = sunline_i / np.linalg.norm(sunline_i)

            # Converting the estimated sun vector to body coordinates
            z_est_sun = DCM @ sunline_i

            # Constructing the sun portion of the H matrix 
            H_sun = np.zeros((3,7))
            H_sun[:,0] = dDCM_dq1 @ sunline_i
            H_sun[:,1] = dDCM_dq2 @ sunline_i
            H_sun[:,2] = dDCM_dq3 @ sunline_i
            H_sun[:,3] = dDCM_dq4 @ sunline_i
            
            # This is removing the body z-axis from the estimated sun vector in case the z panel isn't illuminated
            if I_act[4] < I_min:
                z_act_sun = z_act_sun[0:2]
                z_est_sun = z_est_sun[0:2]
                H_sun = H_sun[0:2,:] / np.linalg.norm(z_est_sun)                
                R = R[:-1,:-1]
            
            # Normalizing our sun vectors
            if np.linalg.norm(z_act_sun) != 0:
                z_act_sun = z_act_sun / np.linalg.norm(z_act_sun)
            if np.linalg.norm(z_est_sun) != 0:
                z_est_sun = z_est_sun / np.linalg.norm(z_est_sun)

            # Appending the actual and estimated sun vectors to their respective measurement vectors
            z_act = np.append(z_act, z_act_sun)
            z_est = np.append(z_est, z_est_sun)

            # Appending the appropriate rows to the H matrix
            H = np.vstack((H, H_sun))

        ## MAGNETIC FIELD VECTOR MODELING
        # Checks if we want to use the magnetometer reading
        if sensors[1] == 1:
            # Adding to measurement noise matrix
            mag_var = 0.05**2
            R = np.block([
                [R, np.zeros((R.shape[0],3))],
                [np.zeros((3,R.shape[1])), np.eye(3)*mag_var]
            ])

            # Adding the magnetometer reading to the actual measurement vector
            z_act_mag = z[5:8]
            if (np.linalg.norm(z_act_mag) != 0):
                z_act_mag = z_act_mag / np.linalg.norm(z_act_mag)
            z_act = np.append(z_act, z_act_mag)

            # Normalizing the magnetic field vector
            if(np.linalg.norm(magField_i) != 0):
                magField_i = magField_i / np.linalg.norm(magField_i)

            # Converting from ECI to body coordinates
            z_est_mag = DCM @ magField_i

            # Appending the estimated magnetometer measurements as the estimated magnetic field vector in body coordinates
            z_est = np.append(z_est, z_est_mag)

            # Appending the appropriate rows to the H matrix       
            H_mag = np.zeros((3,7))
            H_mag[:,0] = dDCM_dq1 @ magField_i
            H_mag[:,1] = dDCM_dq2 @ magField_i
            H_mag[:,2] = dDCM_dq3 @ magField_i
            H_mag[:,3] = dDCM_dq4 @ magField_i

            H = np.vstack((H, H_mag))

        ## GYROSCOPE MODELING
        # Checks if we want to use the gyroscope reading
        if sensors[2] == 1:
            # Adding to measurement noise matrix
            gyro_var = 0.01**2
            R = np.block([
                [R, np.zeros((R.shape[0],3))],
                [np.zeros((3,R.shape[1])), np.eye(3)*gyro_var]
            ])

            # Adding the gyroscope reading to the actual measurement vector
            z_act_gyro = z[8:11]
            z_act = np.append(z_act, z_act_gyro)

            # Appending the estimated gyro measurements as the propagated angular rates
            z_est_gyro = omegaHat
            z_est = np.append(z_est, z_est_gyro)

            # Appending the appropriate rows to the H matrix
            H_gyro = np.array([
                [0,0,0,0,1,0,0],
                [0,0,0,0,0,1,0],
                [0,0,0,0,0,0,1]
                ])
            H = np.vstack((H, H_gyro))

        ## TRUNCATING H MATRIX
        # We computed the full H matrix with 7 columns, but we need a version with only 6
        c_H = H
        c_H = np.delete(c_H, c, axis=1) # Deleting the column corresponding to the redundant quaternion
        
        for i in range(H.shape[0]):
            for j in range(3):
                c_H[i,j] = c_H[i,j] - (qHat[eta(j,c)] / qHat[c])*H[i,c]

        '''
        STATE CORRECTION
        Updating the state vector based on the measurements
        '''
        c_xhat, c_P = ekf_update(c_xhat, c_P, z_act, z_est, c_H, R)

        '''
        UNTRUNCATING
        Recover the full state from the truncated state
        '''
        if np.sum(c_xhat[0:3]**2) < 1:
            q_c = np.sqrt(1 - np.sum(c_xhat[0:3]**2)) * sign
        else:
            q_c = np.abs(qHat[c]) * sign
        
        xhatOut = np.concatenate((c_xhat[:c],[q_c],c_xhat[c:]))

        if np.linalg.norm(xhatOut[0:4]) != 0:
            xhatOut[0:4] = xhatOut[0:4] / np.linalg.norm(xhatOut[0:4])
        qHat = xhatOut[0:4]
        Pout = P_recover(c_P,c,qHat)
    else:
        xhatOut = np.concatenate((qHat,omegaHat))
        Pout = P

    return xhatOut, Pout

'''
OTHER FUNCTIONS
Additional code that was outsourced to functions in order to make the main algorithm flow more smoothly
'''
# Propagates the quaternions in the propagation loop
def qdot_spice(q, I_omega_B, time):
    # This matrix is the conversion between q and q_dot, as in q_dot = 0.5 * omega_matrix * q
    omega_matrix = np.array([
        [0, -I_omega_B[0], -I_omega_B[1], -I_omega_B[2]],
        [I_omega_B[0], 0, I_omega_B[2], -I_omega_B[1]],
        [I_omega_B[1], -I_omega_B[2], 0, I_omega_B[0]],
        [I_omega_B[2], I_omega_B[1], -I_omega_B[0], 0]
    ])

    # Propagate q
    out = expm(0.5 * omega_matrix * time) @ q

    if np.linalg.norm(out) != 0:
        out /= np.linalg.norm(out)
    
    return out

# Used in propagating the angular velocities
def euler_dot(omega, J, Jinv, tau):
    return Jinv @ (tau - np.cross(omega, J @ omega))

# Build the environmental torque vector (tau)
# This isn't currently being used in the EKF since there was not a great deal of confidence in the torque modeling, specifically with the knowledge of the spacecraft's magnetic field
def tau_build(qHat, m_bar, magField_i):
    m_bar = np.array(m_bar, dtype=np.float64)
    magField_i = np.array(magField_i, dtype=np.float64)
    
    # Direction Cosine Matrix to convert vectors in ECI to body coordinates
    DCM = np.array([[1-2*(qHat[2]**2+qHat[3]**2), 2*(qHat[1]*qHat[2]-qHat[0]*qHat[3]), 2*(qHat[1]*qHat[3]+qHat[0]*qHat[2])],
                    [2*(qHat[1]*qHat[2]+qHat[0]*qHat[3]), 1-2*(qHat[1]**2+qHat[3]**2), 2*(qHat[2]*qHat[3]-qHat[0]*qHat[1])], 
                    [2*(qHat[1]*qHat[3]-qHat[2]*qHat[3]), 2*(qHat[2]*qHat[3]+qHat[0]*qHat[1]), 1-2*(qHat[1]**2+qHat[2]**2)]])

    # Torque due to the bar magnet
    magField_body = DCM @ magField_i * 10**-9

    tau_B = np.cross(m_bar, magField_body)

    tau = tau_B

    return tau

# Propagates the truncated covariance matrix
def P_update(P, qHat, omegaHat, deltaT, J, Q):
    omega_matrix = np.array([
        [0, -omegaHat[0], -omegaHat[1], -omegaHat[2]],
        [omegaHat[0], 0, omegaHat[2], -omegaHat[1]],
        [omegaHat[1], -omegaHat[2], 0, omegaHat[0]],
        [omegaHat[2], omegaHat[1], -omegaHat[0], 0]
    ])
    
    skew_matrix = np.array([
        [-qHat[1], -qHat[2], -qHat[3]],
        [qHat[0], -qHat[3], qHat[2]],
        [qHat[3], qHat[0], -qHat[1]],
        [-qHat[2], qHat[1], qHat[0]]
    ])
    
    pi_matrix = np.array([
        [0, omegaHat[2]*(J[1,1]-J[2,2])/J[0,0], omegaHat[1]*(J[1,1]-J[2,2])/J[0,0]],
        [-omegaHat[2]*(J[0,0]-J[2,2])/J[1,1], 0, -omegaHat[0]*(J[0,0]-J[2,2])/J[1,1]],
        [omegaHat[1]*(J[0,0]-J[1,1])/J[2,2], omegaHat[0]*(J[0,0]-J[1,1])/J[2,2], 0]
    ])
    
    F = np.block([
        [np.eye(4) + 0.5*deltaT*omega_matrix, 0.5*deltaT*skew_matrix],
        [np.zeros((3, 4)), np.eye(3) + deltaT*pi_matrix]
    ])

    # Predicted covariance of the estimate
    P = (F @ P @ F.T) + Q

    return P

# Indexing function necessary for the truncation
# Essentially just converts between indices of the truncated state representation and indices of the full state representation based on
# which quaternion component c was removed
def eta(i,c):
    if i < c:
        return i
    else:
        return i+1

# Recovers the full state representation covariance matrix from the truncated state representation covariance matrix
def P_recover(c_P,c,qHat):
    c_qHat = np.delete(qHat,c,axis=0)

    P = np.insert(c_P,(c),0,axis = 0)
    P = np.insert(P,(c),0,axis = 1)

    for j in range(7):
        elem = 0
        if j == c:
            for n in range(3):
                for m in range(3):
                    elem += (c_qHat[n]*c_P[n,m]*c_qHat[m])
            elem /= (qHat[c]**2)
            P[j,j] = elem
        elif j < c:
            for n in range(3):
                elem += (c_qHat[n]*c_P[n,j])
            elem /= (-qHat[c])
            P[c,j] = elem
            P[j,c] = elem
        elif j > c:
            for n in range(3):
                elem += (c_qHat[n]*c_P[n,j-1])
            elem /= (-qHat[c])
            P[c,j] = elem
            P[j,c] = elem
    return P

def ekf_update(xhat, P, z_act, z_est, H, R):
    n = len(xhat)
    
    # Covariance of residual (measurement error)
    S = (H @ P @ H.T) + R
    
    # Near-optimal Kalman gain
    # Note that this was changed from inv to pinv in case S is non-invertible
    K = P @ H.T @ np.linalg.pinv(S)

    # Update the state estimate
    xhat = xhat + K @ (z_act - z_est)
    
    # Updated covariance of the estimate
    P = (np.eye(n) - K @ H) @ P
    
    return xhat, P

# Function for multiplying quaternions (because all python packages related to quaternions hate me and me specifically)
def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    
    w = w1*w2 - x1*x2 - y1*y2 - z1*z2
    x = w1*x2 + x1*w2 + y1*z2 - z1*y2
    y = w1*y2 - x1*z2 + y1*w2 + z1*x2
    z = w1*z2 + x1*y2 - y1*x2 + z1*w2
    
    return np.array([w, x, y, z])

# Function for finding quaternion inverse (because all python packages related to quaternions hate me and me specifically)
def q_inv(q):
    # Extract the scalar (w) and vector part (x, y, z)
    w, x, y, z = q
    
    # Compute the norm (magnitude) of the quaternion
    norm_squared = w**2 + x**2 + y**2 + z**2
    
    # Inverse is conjugate of quaternion divided by norm squared
    return np.array([w, -x, -y, -z]) / norm_squared