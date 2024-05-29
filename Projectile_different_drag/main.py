import numpy as np
import matplotlib.pyplot as plt
import os

def projectile_no_drag(v0, theta, dt, g, xmax):
    t = np.arange(0, xmax / (v0 * np.cos(np.radians(theta))), dt)
    x = v0 * np.cos(np.radians(theta)) * t
    y = v0 * np.sin(np.radians(theta)) * t - 0.5 * g * t**2
    return x, y

def calculate_drag_force(v, angle, C_d, A, rho):
    return 0.5 * rho * v**2 * C_d * A

def projectile_with_drag_euler(v0, theta, dt, g, mass, rho, C_d, A, xmax):
    x = [0]
    y = [0]
    vx = v0 * np.cos(np.radians(theta))
    vy = v0 * np.sin(np.radians(theta))
    t = 0
    while y[-1] >= 0 and x[-1] <= xmax:
        v = np.sqrt(vx**2 + vy**2)
        angle = np.arctan2(vy, vx)
        F_drag = calculate_drag_force(v, angle, C_d, A, rho)
        ax = -F_drag * np.cos(angle) / mass
        ay = -g - F_drag * np.sin(angle) / mass
        vx += ax * dt
        vy += ay * dt
        x.append(x[-1] + vx * dt)
        y.append(y[-1] + vy * dt)
        t += dt
    return np.array(x), np.array(y)

# Define parameters
v0 = 100
theta = 45
g = 9.81
mass = 0.1
rho = 1.225
A = 0.01
xmax = 1000
C_d_values = [0.1, 0.2, 0.3, 0.4, 0.5]  # Different drag coefficients
dt_values = [0.01, 0.05, 0.1, 0.5, 1]

# Directory setup for saving results
directory = 'simulation_results'
if not os.path.exists(directory):
    os.makedirs(directory)

# Analyze for different dt values
for dt in dt_values:
    x_no_drag, y_no_drag = projectile_no_drag(v0, theta, dt, g, xmax)
    plt.figure(figsize=(12, 6))
    plt.plot(x_no_drag, y_no_drag, label='No Drag', linestyle='--')

    for C_d in C_d_values:
        x_drag, y_drag = projectile_with_drag_euler(v0, theta, dt, g, mass, rho, C_d, A, xmax)
        plt.plot(x_drag, y_drag, label=f'C_d = {C_d}')

    plt.title(f'Projectile Motion Simulation with dt = {dt}')
    plt.xlabel('Horizontal Distance (m)')
    plt.ylabel('Vertical Distance (m)')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'{directory}/Projectile_Motion_dt_{dt}.png')
    plt.close()

print("Simulations completed and plots saved.")
