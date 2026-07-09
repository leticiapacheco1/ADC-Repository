import numpy as np
import time
import RPi.GPIO as GPIO
import iimu_test as imu

PWM_PIN = 18
DIR_PIN = 12

GPIO.setmode(GPIO.BCM)
GPIO.setup(PWM_PIN, GPIO.OUT)
GPIO.setup(DIR_PIN, GPIO.OUT)

pwm = GPIO.PWM(PWM_PIN, 1000)
pwm.start(0)

M_MAX = 0.4
I_MAX = 0.34

def send_to_hbridge(m_cmd):
    m_cmd = np.clip(m_cmd, -M_MAX, M_MAX)
    current = (m_cmd / M_MAX) * I_MAX
    GPIO.output(DIR_PIN, GPIO.HIGH if current >= 0 else GPIO.LOW)
    duty = abs(current) / I_MAX * 100
    pwm.ChangeDutyCycle(duty)

m_cmd = 0.4

dt = 0.01
time.sleep(1)

start = time.perf_counter()

while True:
    t = time.perf_counter() - start

    wz = imu.read_gyro()[2]

    send_to_hbridge(m_cmd)

    print(f"t={t:.2f}s | wz={wz:.4f} rad/s | m={m_cmd:.2f} A·m^2")

    time.sleep(dt)