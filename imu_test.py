# Imports the appropriate libraries
import time
import board
import busio
import adafruit_lsm6ds.lsm6dsox
import adafruit_lis3mdl

# Sets up the i2c connection
i2c = board.I2C()

# Sets up the accelerometer and gyroscope
accel_gyro = adafruit_lsm6ds.lsm6dsox.LSM6DSOX(i2c)

# Sets up the magnetometer
# This try statement is necessary as the magnetometer has been occasionally flipping between 0x1c and 0x1e during testing
try:
    mag = adafruit_lis3mdl.LIS3MDL(i2c,0x1e)
except:
    mag = adafruit_lis3mdl.LIS3MDL(i2c)

# Outputs acceleration, gryo, and magnetic data once per second
while True:

    acceleration = accel_gyro.acceleration
    gyro = accel_gyro.gyro
    magnetic = mag.magnetic

    print("Acceleration: X:{0:7.2f}, Y:{1:7.2f}, Z:{2:7.2f} m/s^2".format(*acceleration))
    print("Gyro:         X:{0:7.2f}, Y:{1:7.2f}, Z:{2:7.2f} rad/s".format(*gyro))
    print("Magnetic:     X:{0:7.2f}, Y:{1:7.2f}, Z:{2:7.2f} uT".format(*magnetic))

    print("")
    time.sleep(1)