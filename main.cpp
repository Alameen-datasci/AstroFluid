#include <iostream>
#include <cmath>
#include <vector>
#include "universe.hpp"
#include "gas_material.hpp"
#include "vec3.hpp"

int main() {
    // ---------------------------------------------------------
    // 1. Physical Constants
    // ---------------------------------------------------------
    constexpr double MASS = 1e30;
    constexpr double R = 1e15;
    constexpr double G = 6.67430e-11;
    constexpr double PI = 3.141592653589793;

    // ---------------------------------------------------------
    // 2. Simulation Parameters
    // ---------------------------------------------------------
    constexpr size_t NUM_PARTICLES = 2000;
    constexpr double DT = 1e9;                // Base timestep (seconds)
    constexpr double THETA = 0.5;               // Barnes-Hut opening angle

    std::cout << "--- SPH Engine Physics Test ---" << std::endl;

    // Initialize Universe
    Universe sim(DT, THETA);
    
    // Note: Assuming these are public members in universe.hpp
    sim.enable_vtk_export = true; 
    sim.vtk_export_inteval = 10;

    // ---------------------------------------------------------
    // 3. Calculate SPH Kernel Radius
    // ---------------------------------------------------------
    // The kernel radius needs to be slightly larger than the average
    // particle spacing to ensure sufficient neighbor overlap.
    double volume = (4.0 / 3.0) * PI * std::pow(R, 3);
    double particle_volume = volume / NUM_PARTICLES;
    double spacing = std::cbrt(particle_volume);
    double kernel_radius = 1.5 * spacing; 

    // ---------------------------------------------------------
    // 4. Register Material (Monatomic Hydrogen Gas)
    // ---------------------------------------------------------
    int mat_id = sim.registerMaterial(
        kernel_radius, 
        1,
        1.4,
        0.002
    );

    double initial_temp = 3.211113784;
    
    std::cout << "Smoothing Length (Kernel): " << kernel_radius << " m" << std::endl;

    // ---------------------------------------------------------
    // 6. Generate the Cloud
    // ---------------------------------------------------------
    std::cout << "\nInitializing Cloud..." << std::endl;
    sim.addGasCloud(
        NUM_PARTICLES,
        MASS,
        R,
        {0.0, 0.0, 0.0}, // Center Position
        {0.0, 0.0, 0.0}, // Initial Velocity
        initial_temp,     // Initial Temperature
        mat_id,          // Material ID
        1e12           // Relaxation Time (s)
    );

    // ---------------------------------------------------------
    // 7. Run Main Simulation
    // ---------------------------------------------------------
    double simulation_time = 2e12;
    std::cout << "\nStarting main simulation run for " << simulation_time << " seconds..." << std::endl;
    
    sim.run(simulation_time);

    std::cout << "Simulation Complete! Check VTK outputs." << std::endl;
    return 0;
}