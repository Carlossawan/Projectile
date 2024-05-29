
# Projectile Motion Simulation

## Project Overview

This project simulates projectile motion considering the effects of air resistance using two numerical methods: Euler and Runge-Kutta (RK4). The primary objectives are to understand the fundamental principles behind these numerical methods, implement projectile motion models with air resistance, and compare their accuracy, stability, and computational efficiency. Visualizations of the trajectories are provided to illustrate the differences between the methods.

## Objectives

- **Understand the Basics**: Grasp the fundamental principles behind the Euler and Runge-Kutta (RK4) numerical methods for solving ordinary differential equations (ODEs).
- **Implement Projectile Motion Models**: Develop simulation models for projectile motion that include the effects of air resistance, using both the Euler and RK4 methods.
- **Incorporate Air Resistance**: Integrate a model of air resistance into the simulations, using both linear and quadratic resistance models.
- **Compare Numerical Methods**: Evaluate and compare the accuracy, stability, and computational efficiency of the Euler and RK4 methods.
- **Visualize Trajectories**: Create visual representations of the projectile trajectories.
- **Analyze the Impact of Air Resistance**: Investigate how varying degrees of air resistance affect the projectile's motion.
- **Evaluate Performance**: Assess the performance of both methods in terms of computation time and accuracy.
- **Explore Parameter Sensitivity**: Examine how changes in initial conditions and model parameters impact the simulations.



## C++ Implementation

### Prerequisites

- A C++ compiler (e.g., g++)

## Python Implementation

### Prerequisites

- Python 3.x
- Anaconda
- Required libraries: NumPy, Matplotlib

### Installing Dependencies

Install the required libraries using pip:
```bash
conda activate projectile
conda install numpy matplotlib
```

### Running the Simulation

1. Navigate to the `python` directory:
    ```bash
    cd python
    ```
2. Run the main script:
    ```bash
    python main.py
    ```



## Results and Analysis

- **Accuracy and Stability**: The RK4 method provides higher accuracy and stability compared to the Euler method, especially for smaller step sizes.
- **Computational Efficiency**: The Euler method is computationally less intensive but less accurate. The RK4 method, while more computationally demanding, offers greater precision.
- **Impact of Air Resistance**: The simulations show significant differences in trajectories when air resistance is considered. The quadratic model provides a more realistic representation of air resistance effects.

## Visualizations

Visualizations of the projectile trajectories generated by both the Euler and RK4 methods are included in the project. These visualizations highlight the differences in the paths and final locations due to the chosen numerical method and the effects of air resistance.
