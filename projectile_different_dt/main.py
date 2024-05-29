import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def calculate_drag_force(v, angle, C_d, A, rho):
    """Calculate the drag force given velocity, angle, drag coefficient, area, and air density."""
    return 0.5 * rho * v**2 * C_d * A

def projectile_with_drag_rk4(v0, theta, dt, g, mass, rho, C_d, A, xmax):
    """Simulate projectile motion with air drag using RK4 method."""
    def f(t, state):
        """Function defining the differential equations."""
        x, y, vx, vy = state
        v = np.sqrt(vx**2 + vy**2)
        angle = np.arctan2(vy, vx)
        F_drag = calculate_drag_force(v, angle, C_d, A, rho)
        ax = -F_drag * np.cos(angle) / mass
        ay = -g - F_drag * np.sin(angle) / mass
        return [vx, vy, ax, ay]

    x, y, vx, vy = [0], [0], [v0 * np.cos(np.radians(theta))], [v0 * np.sin(np.radians(theta))]
    while y[-1] >= 0 and x[-1] <= xmax:
        state = [x[-1], y[-1], vx[-1], vy[-1]]
        k1 = dt * np.array(f(0, state))
        k2 = dt * np.array(f(0, state + 0.5 * k1))
        k3 = dt * np.array(f(0, state + 0.5 * k2))
        k4 = dt * np.array(f(0, state + k3))
        new_state = state + (1/6) * (k1 + 2*k2 + 2*k3 + k4)
        x.append(new_state[0])
        y.append(new_state[1])
        vx.append(new_state[2])
        vy.append(new_state[3])
    return np.array(x), np.array(y)

# Constants and initial conditions
v0 = 100
theta = 45
g = 9.81
mass = 10
rho = 1.225
C_d = 0.47
A = 0.05
xmax = 1000
dt_values = sorted([0.0001 ,0.001 ,0.01, 0.05, 0.1, 0.5, 1])

trajectories = {}
plt.figure(figsize=(12, 6))

for dt in dt_values:
    x, y = projectile_with_drag_rk4(v0, theta, dt, g, mass, rho, C_d, A, xmax)
    trajectories[dt] = (x, y)
    plt.plot(x, y, label=f'dt = {dt}')

# Calculate the error between each dt and the next finer dt
errors = {}
for i in range(len(dt_values) - 1):
    current_dt = dt_values[i]
    next_finer_dt = dt_values[i+1]
    current_x, current_y = trajectories[current_dt]
    next_x, next_y = trajectories[next_finer_dt]
    
    interp_func = interp1d(current_x, current_y, kind='linear', fill_value="extrapolate", bounds_error=False)
    interpolated_y = interp_func(next_x)
    error = np.mean(np.sqrt((next_x - next_x)**2 + (next_y - interpolated_y)**2))
    errors[f"{current_dt} to {next_finer_dt}"] = error

print("Incremental Average Euclidean Errors:")
for dt_pair, error in errors.items():
    print(f"Error from {dt_pair}: {error}")

plt.title('Projectile Motion with Air Drag for Different dt Values (RK4)')
plt.xlabel('Horizontal Distance (m)')
plt.ylabel('Vertical Distance (m)')
plt.legend()
plt.grid(True)
plt.show()
