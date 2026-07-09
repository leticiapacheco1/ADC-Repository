import numpy as np
import pandas as pd
import time
import os
import RPi.GPIO as GPIO

import updated_imu_test as imu


class SimpleKalman:
    def __init__(self, q=0.01, r=0.1):
        self.q = q
        self.r = r
        self.x = None
        self.p = 1

    def update(self, measurement):
        if self.x is None:
            self.x = measurement

        self.p += self.q
        k = self.p / (self.p + self.r)

        self.x += k * (measurement - self.x)
        self.p *= (1 - k)

        return self.x


class GyroController3D:
    def __init__(self,
                 K=100000.0,
                 M_MAX=0.423):

        self.K = K
        self.M_MAX = M_MAX

    def update(self, B, B_prev, dt):

        if dt <= 0:
            return np.zeros(3)

        B_dot = (B - B_prev) / dt

        m_cmd = -self.K * B_dot

        return np.clip(m_cmd, -self.M_MAX, self.M_MAX)


PWM_PIN = 18   #1 DIR pin + 1 PWM pin -> controlling ALL axes combined. This will change later because we will have 3 PWM and 3 DIR pins for each axis.
DIR_PIN = 12

GPIO.setmode(GPIO.BCM)
GPIO.setup(PWM_PIN, GPIO.OUT)
GPIO.setup(DIR_PIN, GPIO.OUT)

pwm = GPIO.PWM(PWM_PIN, 1000)
pwm.start(0)

I_MAX = 0.34
M_MAX = 0.423


def send_to_hbridge(m_cmd):

    m_mag = np.linalg.norm(m_cmd)
    axis = np.argmax(np.abs(m_cmd))

    sign = np.sign(m_cmd[axis])

    m_mag = np.clip(m_mag, -M_MAX, M_MAX)

    cmd_norm = m_mag / M_MAX

    current = cmd_norm * I_MAX
    current = np.clip(current, -I_MAX, I_MAX)

    if sign >= 0:
        GPIO.output(DIR_PIN, GPIO.HIGH)
    else:
        GPIO.output(DIR_PIN, GPIO.LOW)

    duty = (abs(current) / I_MAX) * 100

    if abs(current) < 0.01:
        duty = 0

    pwm.ChangeDutyCycle(duty)


def read_mag():
    bx, by, bz = imu.read_magnetometer_microT()
    return np.array([bx, by, bz]) * 1e-6


def read_gyro():
    return np.degrees(np.array(imu.read_gyro()))


if __name__ == "__main__":

    print("Starting ADCS B-dot detumbling system")
    print("State tracking enabled (5 deg/s threshold)\n")

    controller = GyroController3D()
    kalman = SimpleKalman()

    dt = 0.05
    time.sleep(1)

    B_prev = read_mag()

    start_time = time.perf_counter()
    last_print = start_time

    try:
        while True:

            loop_start = time.perf_counter()

            B = read_mag()
            gyro = read_gyro()

            wx, wy, wz = gyro

            if (abs(wx) > 5.0) or (abs(wy) > 5.0) or (abs(wz) > 5.0):
                state = "SPINNING"
            else:
                state = "STABLE"

            m_cmd = controller.update(B, B_prev, dt)

            send_to_hbridge(m_cmd)

            t = loop_start - start_time

            if loop_start - last_print >= 1.0:
                last_print = loop_start

                print(f"Time (s): {t:.2f}")
                print(f"State: {state}")

                print(f"ωx: {wx:8.3f} deg/s | ωy: {wy:8.3f} | ωz: {wz:8.3f}")

                print(
                    f"mx: {m_cmd[0]:.6f} | "
                    f"my: {m_cmd[1]:.6f} | "
                    f"mz: {m_cmd[2]:.6f} A·m²"
                )

                print("-----------------------------")

            B_prev = B

            elapsed = time.perf_counter() - loop_start
            if dt - elapsed > 0:
                time.sleep(dt - elapsed)

    except KeyboardInterrupt:
        print("Stopping system")

    finally:
        pwm.stop()
        GPIO.cleanup()