#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <filesystem>
namespace fs = std::filesystem;

struct Vector2 {
    double x, y;
};

const double PI = 3.14159265358979323846;
const double G = 9.81; // gravitational acceleration in m/s^2
const double rho = 1.225; // air density in kg/m^3

// Convert degrees to radians
double toRadians(double theta) {
    return theta * PI / 180.0;
}

// Projectile motion without air drag
std::vector<Vector2> projectile_no_drag(double v0, double theta, double dt, double xmax) {
    std::vector<Vector2> trajectory;
    double t = 0;
    double vx = v0 * std::cos(toRadians(theta));
    double vy = v0 * std::sin(toRadians(theta));
    while (vx * t <= xmax) {
        double x = vx * t;
        double y = vy * t - 0.5 * G * t * t;
        if (y < 0) break;
        trajectory.push_back({ x, y });
        t += dt;
    }
    return trajectory;
}

// Compute drag force vector components
Vector2 calculate_drag_force(double v, double angle, double C_d, double A) {
    double drag_magnitude = 0.5 * rho * v * v * C_d * A;
    return { drag_magnitude * std::cos(angle), drag_magnitude * std::sin(angle) };
}

// Projectile motion with air drag using Euler’s method
std::vector<Vector2> projectile_with_drag_euler(double v0, double theta, double dt, double mass, double C_d, double A, double xmax) {
    std::vector<Vector2> trajectory;
    Vector2 position = { 0, 0 };
    Vector2 velocity = { v0 * std::cos(toRadians(theta)), v0 * std::sin(toRadians(theta)) };

    while (position.y >= 0) {
        trajectory.push_back(position);
        double v = std::hypot(velocity.x, velocity.y);
        double angle = std::atan2(velocity.y, velocity.x);
        Vector2 drag_force = calculate_drag_force(v, angle, C_d, A);
        Vector2 acceleration = {
                -drag_force.x / mass,
                -G - drag_force.y / mass
        };
        velocity.x += acceleration.x * dt;
        velocity.y += acceleration.y * dt;
        position.x += velocity.x * dt;
        position.y += velocity.y * dt;
        if (position.x > xmax) break;
    }

    return trajectory;
}

// Projectile motion with air drag using RK2 method
std::vector<Vector2> projectile_with_drag_RK2(double v0, double theta, double dt, double mass, double C_d, double A, double xmax) {
    std::vector<Vector2> trajectory;
    Vector2 position = { 0, 0 };
    Vector2 velocity = { v0 * std::cos(toRadians(theta)), v0 * std::sin(toRadians(theta)) };

    while (position.y >= 0) {
        trajectory.push_back(position);
        double v = std::hypot(velocity.x, velocity.y);
        double angle = std::atan2(velocity.y, velocity.x);
        Vector2 drag_force = calculate_drag_force(v, angle, C_d, A);
        Vector2 acceleration = {
                -drag_force.x / mass,
                -G - drag_force.y / mass
        };
        // First step
        Vector2 k1_vel = { velocity.x + acceleration.x * dt, velocity.y + acceleration.y * dt };
        Vector2 k1_pos = { position.x + velocity.x * dt, position.y + velocity.y * dt };

        // Second step
        double k1_v = std::hypot(k1_vel.x, k1_vel.y);
        double k1_angle = std::atan2(k1_vel.y, k1_vel.x);
        Vector2 k1_drag_force = calculate_drag_force(k1_v, k1_angle, C_d, A);
        Vector2 k1_acceleration = { -k1_drag_force.x / mass, -G - k1_drag_force.y / mass };
        Vector2 k2_vel = { velocity.x + k1_acceleration.x * dt, velocity.y + k1_acceleration.y * dt };
        Vector2 k2_pos = { position.x + k1_vel.x * dt, position.y + k1_vel.y * dt }; // Corrected position calculation using k1 velocity

        // Update position and velocity using the average of k1 and k2
        velocity.x += (k1_acceleration.x + (-calculate_drag_force(std::hypot(k2_vel.x, k2_vel.y), std::atan2(k2_vel.y, k2_vel.x), C_d, A).x / mass)) * 0.5 * dt;
        velocity.y += (k1_acceleration.y + (-G - calculate_drag_force(std::hypot(k2_vel.x, k2_vel.y), std::atan2(k2_vel.y, k2_vel.x), C_d, A).y / mass)) * 0.5 * dt;
        position.x += (k1_vel.x + k2_vel.x) * 0.5 * dt;
        position.y += (k1_vel.y + k2_vel.y) * 0.5 * dt;

        if (position.x > xmax) break;
    }

    return trajectory;
}

// Projectile motion with air drag using RK4 method
std::vector<Vector2> projectile_with_drag_RK4(double v0, double theta, double dt, double mass, double C_d, double A, double xmax) {
    std::vector<Vector2> trajectory;
    Vector2 position = {0, 0};
    Vector2 velocity = {v0 * std::cos(toRadians(theta)), v0 * std::sin(toRadians(theta))};

    while (position.y >= 0) {
        trajectory.push_back(position);
        double v = std::hypot(velocity.x, velocity.y);
        double angle = std::atan2(velocity.y, velocity.x);

        // Compute the drag force at the beginning of the interval
        Vector2 drag_force = calculate_drag_force(v, angle, C_d, A);
        Vector2 k1_acc = {
                -drag_force.x / mass,
                -G - drag_force.y / mass
        };
        Vector2 k1_vel = {
                velocity.x + k1_acc.x * dt,
                velocity.y + k1_acc.y * dt
        };
        Vector2 k1_pos = {
                position.x + velocity.x * dt,
                position.y + velocity.y * dt
        };

        // Midpoint
        v = std::hypot(k1_vel.x, k1_vel.y);
        angle = std::atan2(k1_vel.y, k1_vel.x);
        drag_force = calculate_drag_force(v, angle, C_d, A);
        Vector2 k2_acc = {
                -drag_force.x / mass,
                -G - drag_force.y / mass
        };
        Vector2 k2_vel = {
                velocity.x + k2_acc.x * dt * 0.5,
                velocity.y + k2_acc.y * dt * 0.5
        };
        Vector2 k2_pos = {
                position.x + k2_vel.x * dt * 0.5,
                position.y + k2_vel.y * dt * 0.5
        };

        // Another midpoint using k2 values
        v = std::hypot(k2_vel.x, k2_vel.y);
        angle = std::atan2(k2_vel.y, k2_vel.x);
        drag_force = calculate_drag_force(v, angle, C_d, A);
        Vector2 k3_acc = {
                -drag_force.x / mass,
                -G - drag_force.y / mass
        };
        Vector2 k3_vel = {
                velocity.x + k3_acc.x * dt * 0.5,
                velocity.y + k3_acc.y * dt * 0.5
        };
        Vector2 k3_pos = {
                position.x + k3_vel.x * dt * 0.5,
                position.y + k3_vel.y * dt * 0.5
        };

        // End point using k3 values
        v = std::hypot(k3_vel.x, k3_vel.y);
        angle = std::atan2(k3_vel.y, k3_vel.x);
        drag_force = calculate_drag_force(v, angle, C_d, A);
        Vector2 k4_acc = {
                -drag_force.x / mass,
                -G - drag_force.y / mass
        };
        Vector2 k4_vel = {
                velocity.x + k4_acc.x * dt,
                velocity.y + k4_acc.y * dt
        };
        Vector2 k4_pos = {
                position.x + k4_vel.x * dt,
                position.y + k4_vel.y * dt
        };

        // Combine all increments to get new position and velocity
        position.x += (k1_pos.x + 2*k2_pos.x + 2*k3_pos.x + k4_pos.x) / 6;
        position.y += (k1_pos.y + 2*k2_pos.y + 2*k3_pos.y + k4_pos.y) / 6;
        velocity.x += (k1_vel.x + 2*k2_vel.x + 2*k3_vel.x + k4_vel.x) / 6;
        velocity.y += (k1_vel.y + 2*k2_vel.y + 2*k3_vel.y + k4_vel.y) / 6;

        if (position.x > xmax) break;
    }

    return trajectory;
}

int main() {
    // Example usage
    double v0 = 10.0; // initial velocity (m/s)
    double theta = 45.0; // launch angle (degrees)
    double dt = 0.01; // time step (s)
    double mass = 1.0; // mass of the projectile (kg)
    double C_d = 0.1; // drag coefficient
    double A = 0.01; // cross-sectional area (m^2)
    double xmax = 100.0; // maximum range (m)


    std::vector<Vector2> trajectory_no_drag = projectile_no_drag(v0, theta, dt, xmax);
    std::vector<Vector2> trajectory_drag_euler = projectile_with_drag_euler(v0, theta, dt, mass, C_d, A, xmax);
    std::vector<Vector2> trajectory_drag_RK2 = projectile_with_drag_RK2(v0, theta, dt, mass, C_d, A, xmax);
    std::vector<Vector2> trajectory_drag_RK4 = projectile_with_drag_RK4(v0 , theta ,dt,mass ,C_d ,A , xmax );

    std::cout << "Trajectory with euler method:" << std::endl;
    for (const auto& point : trajectory_drag_euler) {
        std::cout << "(" << point.x << ", " << point.y << ")" << std::endl;
    }


    std::cout << "Trajectory with RK2 method:" << std::endl;
    for (const auto& point : trajectory_drag_RK2) {
        std::cout << "(" << point.x << ", " << point.y << ")" << std::endl;
    }

    std::cout << "Trajectory with RK4 method:" << std::endl;
    for (const auto& point : trajectory_drag_RK4) {
        std::cout << "(" << point.x << ", " << point.y << ")" << std::endl;
    }
    return 0;
}