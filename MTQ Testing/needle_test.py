import numpy as np
import time
import iimu_test as imu

def average_mag(n=50):
    readings = []
    for _ in range(n):
        readings.append(imu.read_magnetometer_microT())
        time.sleep(0.02)
    return np.mean(readings, axis=0)

B_EARTH_HORIZONTAL = 20.0  # µT, St. Louis horizontal component

print("=== SETUP ===")
print("1. Place IMU flat on table")
print("2. Place iPhone compass next to IMU")
print("3. Align MTQ rod East-West (perpendicular to compass needle)")
print("4. Make sure compass reads 0° North before starting")
print("5. MTQ is OFF right now")
input("\nPress Enter to record baseline...")

baseline = average_mag(n=50)
print(f"Baseline recorded: Bx={baseline[0]:.2f} By={baseline[1]:.2f} Bz={baseline[2]:.2f} µT")

results = []

while True:
    print("\n--- New measurement ---")
    d_str = input("Enter distance d in cm (or 'done' to finish): ")
    if d_str.lower() == 'done':
        break

    d_cm = float(d_str)

    input(f"Hold MTQ at {d_cm}cm above IMU, power it ON, then press Enter...")
    reading = average_mag(n=50)

    delta = reading - baseline
    B_mtq = np.linalg.norm(delta)

    theta_predicted = np.degrees(np.arctan(B_mtq / B_EARTH_HORIZONTAL))

    print(f"\nIMU measured ΔB = {B_mtq:.2f} µT")
    print(f"Predicted compass deflection = {theta_predicted:.1f}°")

    print("\nNow read the iPhone compass.")
    print("Enter the angle shown (0-360°): ", end="")
    angle_measured = float(input())

    print("Enter needle quadrant (NE, NW, SE, SW): ", end="")
    quadrant = input().strip().upper()

    # Deflection from North (0°)
    # NE/NW means needle swung toward East or West from North
    if quadrant in ["NE", "NW"]:
        theta_measured = angle_measured  # already deflection from North
    else:
        theta_measured = 180 - angle_measured  # SE/SW, deflected past 90

    error = abs(theta_measured - theta_predicted)

    print(f"\n=== RESULT at d={d_cm}cm ===")
    print(f"  ΔB from IMU:          {B_mtq:.2f} µT")
    print(f"  Predicted θ:          {theta_predicted:.1f}°")
    print(f"  Measured θ:           {theta_measured:.1f}°")
    print(f"  Error:                {error:.1f}°")

    # Compute dipole moment from dipole field formula
    # B = (µ0/4π) * (2m/r³) along dipole axis
    r = d_cm / 100  # convert to meters
    mu0_over_4pi = 1e-7
    m_measured = (B_mtq * 1e-6) * r**3 / (2 * mu0_over_4pi)
    print(f"  Dipole moment m:      {m_measured:.4f} A·m²")
    print(f"  Expected m:           0.4153 A·m²")

    results.append({
        'd_cm': d_cm,
        'B_mtq_uT': B_mtq,
        'theta_predicted': theta_predicted,
        'theta_measured': theta_measured,
        'error_deg': error,
        'm_measured': m_measured
    })

print("\n=== FULL SUMMARY ===")
print(f"{'d (cm)':>8} | {'ΔB (µT)':>10} | {'θ pred':>8} | {'θ meas':>8} | {'error':>7} | {'m (A·m²)':>10}")
print("-" * 65)
for r in results:
    print(f"{r['d_cm']:>8.1f} | {r['B_mtq_uT']:>10.2f} | "
          f"{r['theta_predicted']:>8.1f} | {r['theta_measured']:>8.1f} | "
          f"{r['error_deg']:>7.1f} | {r['m_measured']:>10.4f}")

m_values = [r['m_measured'] for r in results]
print(f"\nAverage measured dipole moment: {np.mean(m_values):.4f} A·m²")
print(f"Expected from characterization: 0.4153 A·m²")