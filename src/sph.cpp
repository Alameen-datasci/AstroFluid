#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include "sph.hpp"
#include "vec3.hpp"
#include "spatial_hash.hpp"


SPHSolver::SPHSolver(
    const std::vector<GasMaterial>& materials_,
    const std::vector<int>& material_ids_,
    const std::vector<double>& masses_,
    const std::vector<size_t>& activeParticles_,
    std::vector<double>& pos_x_, std::vector<double>& pos_y_, std::vector<double>& pos_z_,
    std::vector<double>& vel_x_, std::vector<double>& vel_y_, std::vector<double>& vel_z_,
    std::vector<double>& acc_x_, std::vector<double>& acc_y_, std::vector<double>& acc_z_,
    std::vector<double>& ext_acc_x_, std::vector<double>& ext_acc_y_, std::vector<double>& ext_acc_z_,
    std::vector<double>& density_, std::vector<double>& pressure_,
    std::vector<double>& u_, std::vector<double>& du_dt_)
    : materials(materials_), material_ids(material_ids_), masses(masses_), activeParticles(activeParticles_),
    pos_x(pos_x_), pos_y(pos_y_), pos_z(pos_z_),
    vel_x(vel_x_), vel_y(vel_y_), vel_z(vel_z_),
    acc_x(acc_x_), acc_y(acc_y_), acc_z(acc_z_),
    ext_acc_x(ext_acc_x_), ext_acc_y(ext_acc_y_), ext_acc_z(ext_acc_z_),
    density(density_), pressure(pressure_), u(u_), du_dt(du_dt_) {};

void SPHSolver::computeDensity(SpatialHash& grid) {
    #pragma omp parallel for
    for (size_t l = 0; l < activeParticles.size(); l++) {
        int i = activeParticles[l];
        double local_density = 0.0;
        const GasMaterial& mat_i = materials[material_ids[i]];

        CellCord center = grid.toCell(pos_x[i], pos_y[i], pos_z[i]);

        for (int cx = -1; cx <= 1; cx++) {
            for (int cy = -1; cy <= 1; cy++) {
                for (int cz = -1; cz <= 1; cz++) {
                    CellCord neighbour_cell = {center.x + cx, center.y + cy, center.z + cz};

                    for (int j = grid.getFirst(neighbour_cell); j != -1; j = grid.getNext(j)) {
                        const GasMaterial& mat_j = materials[material_ids[j]];

                        double dx = pos_x[j] - pos_x[i];
                        double dy = pos_y[j] - pos_y[i];
                        double dz = pos_z[j] - pos_z[i];

                        double r2 = dx*dx + dy*dy + dz*dz;
                        double h_ij = (mat_i.kernel_radius + mat_j.kernel_radius) * 0.5;
                        double h_ij2 = h_ij * h_ij;

                        double mask = (r2 < h_ij2) ? 1.0 : 0.0;

                        double h9 = h_ij2 * h_ij2 * h_ij2 * h_ij2 * h_ij;
                        double inv_denom = 64.0 * 3.141592653589793 * h9;
                        double factor = 315.0 / inv_denom;

                        double value = h_ij2 - r2;
                        double value3 = value * value * value;

                        local_density += (masses[j] * factor * value3) * mask;
                    }
                }
            }
        }
        double h_i = mat_i.kernel_radius;
        double volume = (4.0 / 3.0) * 3.141592653589793 * h_i * h_i * h_i;
        double min_density = masses[i] / volume;
        
        density[i] = std::max(local_density, min_density * 0.01);
    }
}

void SPHSolver::computePressure() {
    #pragma omp parallel for simd
    for (size_t i = 0; i < density.size(); i++) {
        const GasMaterial& mat = materials[material_ids[i]];

        // Ideal Gas Equation for SPH: P = (gamma - 1.0) * rho * u
        pressure[i] = (mat.gamma - 1.0) * density[i] * u[i];
    }
}

void SPHSolver::computeAccelerations(SpatialHash& grid) {
    #pragma omp parallel for
    for (size_t l = 0; l < activeParticles.size(); l++) {
        int i = activeParticles[l];
        double p_acc_x = 0.0, p_acc_y = 0.0, p_acc_z = 0.0;
        double v_acc_x = 0.0, v_acc_y = 0.0, v_acc_z = 0.0;
        double du_i = 0.0;

        int id_i = material_ids[i];
        const GasMaterial& mat_i = materials[id_i];

        CellCord center = grid.toCell(pos_x[i], pos_y[i], pos_z[i]);

        for (int cx = -1; cx <= 1; cx++) {
            for (int cy = -1; cy <= 1; cy++) {
                for (int cz = -1; cz <= 1; cz++) {
                    CellCord neighbour_cell = {center.x + cx, center.y + cy, center.z + cz};

                    for (int j = grid.getFirst(neighbour_cell); j != -1; j = grid.getNext(j)) {
                        int id_j = material_ids[j];
                        const GasMaterial& mat_j = materials[id_j];

                        double dx = pos_x[j] - pos_x[i];
                        double dy = pos_y[j] - pos_y[i];
                        double dz = pos_z[j] - pos_z[i];

                        double r2 = dx*dx + dy*dy + dz*dz;
                        double h_ij = (mat_i.kernel_radius + mat_j.kernel_radius) * 0.5;
                        double r = std::max(std::sqrt(r2), 0.5 * h_ij);
                        double h_ij2 = h_ij * h_ij;
                        double h_ij3 = h_ij * h_ij * h_ij;
                        double h_ij6 = h_ij3 * h_ij3;

                        double mask = (r2 < h_ij2 && i != j) ? 1.0 : 0.0;
                        
                        // pressure force calculation
                        double inv_denom_grad = 3.1415926535 * h_ij6 * r;
                        double grad_factor = -45.0 / inv_denom_grad;
                        double value = h_ij - r;
                        double grad_w_mag = grad_factor * value * value;
                        double pressure_term = (pressure[i] / (density[i] * density[i])) +
                        (pressure[j] / (density[j] * density[j]));

                        double p_force = masses[j] * pressure_term * grad_w_mag * mask;

                        // viscosity force calculations
                        double inv_denom_lap = 3.1415926535 * h_ij6;
                        double lap_factor = 45.0 / inv_denom_lap;
                        double laplacian = lap_factor * value;

                        double visc_ij = (mat_i.viscosity + mat_j.viscosity) * 0.5;
                        double dv_x = vel_x[j] - vel_x[i];
                        double dv_y = vel_y[j] - vel_y[i];
                        double dv_z = vel_z[j] - vel_z[i];

                        double v_force = (visc_ij * masses[j] * laplacian / density[j]) * mask;

                        // internal energy calculation
                        double dv_dot_dx = dv_x * dx + dv_y * dy + dv_z * dz;
                        double dv_mag2 = dv_x * dv_x + dv_y * dv_y + dv_z * dv_z;

                        // mechanical work (compression heats, expansion cools)
                        double energy_change_pressure = 0.5 * p_force * dv_dot_dx;

                        // viscosity friction (frictions causes heating)
                        double energy_change_viscosity = 0.5 * v_force * dv_mag2;

                        du_i += energy_change_pressure + energy_change_viscosity;

                        p_acc_x += p_force * dx;
                        p_acc_y += p_force * dy;
                        p_acc_z += p_force * dz;

                        v_acc_x += v_force * dv_x;
                        v_acc_y += v_force * dv_y;
                        v_acc_z += v_force * dv_z;
                    }
                }
            }
        }

        acc_x[i] = p_acc_x + v_acc_x;
        acc_y[i] = p_acc_y + v_acc_y;
        acc_z[i] = p_acc_z + v_acc_z;

        du_dt[i] = du_i;
    }
}

void SPHSolver::step(double dt, double damping_factor) {
    #pragma omp parallel for
    for (size_t k = 0; k < activeParticles.size(); k++) {
        int i = activeParticles[k];
        acc_x[i] = 0.0;
        acc_y[i] = 0.0;
        acc_z[i] = 0.0;
    }

    double h_max = 0.0;
    for (auto& m : materials) {
        h_max = std::max(h_max, m.kernel_radius);
    }

    SpatialHash grid(h_max, pos_x.size());
    grid.build(pos_x, pos_y, pos_z);

    computeDensity(grid);
    computePressure();
    computeAccelerations(grid);
    
    #pragma omp parallel for simd
    for (size_t k = 0; k < activeParticles.size(); k++) {
        int i = activeParticles[k];
        // Kick 2
        vel_x[i] += 0.5 * (acc_x[i] + ext_acc_x[i]) * dt;
        vel_y[i] += 0.5 * (acc_y[i] + ext_acc_y[i]) * dt;
        vel_z[i] += 0.5 * (acc_z[i] + ext_acc_z[i]) * dt;

        vel_x[i] *= damping_factor;
        vel_y[i] *= damping_factor;
        vel_z[i] *= damping_factor;

        u[i] += 0.5 * du_dt[i] * dt;
        const GasMaterial& mat = materials[material_ids[i]];
        double cv = 8.314 / (mat.molar_mass * (mat.gamma - 1.0));
        double u_min = cv * 2.7;
        u[i] = std::max(u[i], u_min);
    }
}