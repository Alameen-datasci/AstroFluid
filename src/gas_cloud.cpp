#include <cmath>
#include<random>
#include "gas_cloud.hpp"


GasCloud::GasCloud(
    size_t num_particles_, double total_mass_, double radius_, std::vector<double>& masses_,
    int k_,
    std::vector<double>& pos_x_, std::vector<double>& pos_y_, std::vector<double>& pos_z_,
    std::vector<double>& vel_x_, std::vector<double>& vel_y_, std::vector<double>& vel_z_
) : num_particles(num_particles_), total_mass(total_mass_), R(radius_), masses(masses_), k(k_),
    pos_x(pos_x_), pos_y(pos_y_), pos_z(pos_z_), vel_x(vel_x_), vel_y(vel_y_), vel_z(vel_z_)
{
    pos_x.reserve(num_particles);
    pos_y.reserve(num_particles);
    pos_z.reserve(num_particles);

    generateParticle();

    num_particles = pos_x.size();

    vel_x.resize(num_particles, 0.0);
    vel_y.resize(num_particles, 0.0);
    vel_z.resize(num_particles, 0.0);

    masses.resize(num_particles, total_mass / num_particles);
}

void GasCloud::generateParticle() {
    double volume = (4.0 / 3.0) * 3.14159265358979323846 * R * R * R;
    double particle_volume = volume / num_particles;
    r = 0.9 * std::cbrt(particle_volume);
    cell_size = r / std::sqrt(3.0);

    Nx = std::ceil(2.0 * R / cell_size);
    Ny = std::ceil(2.0 * R / cell_size);
    Nz = std::ceil(2.0 * R / cell_size);

    grid.resize(Nx * Ny * Nz, -1);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> dist_z(-1.0, 1.0);
    std::uniform_real_distribution<double> dist_theta(0.0, 2.0 * 3.14159265358979323846);
    std::uniform_real_distribution<double> dist_mag(r, 2.0 * r);

    vec3 center = vec3::zero;
    pos_x.push_back(center.x);
    pos_y.push_back(center.y);
    pos_z.push_back(center.z);

    CellCord center_cell = toCell(center.x, center.y, center.z);

    active.push_back(0);
    
    int center_idx = Nx * Ny * center_cell.z + Nx * center_cell.y + center_cell.x;
    grid[center_idx] = 0;

    while (!active.empty()) {
        std::uniform_int_distribution<int> dist_active(0, active.size() - 1);
        int activeIndex = dist_active(gen);
        int host = active[activeIndex];

        vec3 HostPos = {
            pos_x[host],
            pos_y[host],
            pos_z[host]
        };
        
        bool point_found = false;

        CellCord HostCell = toCell(
            HostPos.x,
            HostPos.y,
            HostPos.z
        );

        for (size_t attempt = 0; attempt < k; attempt++) {
            double z_val = dist_z(gen);
            double theta = dist_theta(gen);
            double r_xy = std::sqrt(1.0 - z_val * z_val);
            double mag = dist_mag(gen);
            vec3 unit_dir = {r_xy * std::cos(theta), r_xy * std::sin(theta), z_val};
            vec3 point = HostPos + unit_dir * mag;

            if ((point.x * point.x + point.y * point.y + point.z * point.z > R * R)) continue;

            CellCord CandidateCell = toCell(point.x, point.y, point.z);

            if (
                CandidateCell.x < 0 || CandidateCell.x >= Nx ||
                CandidateCell.y < 0 || CandidateCell.y >= Ny ||
                CandidateCell.z < 0 || CandidateCell.z >= Nz
            ) continue;

            bool valid = true;

            for (int cx = -2; cx <= 2; cx++) {
                for (int cy = -2; cy <= 2; cy++) {
                    for (int cz = -2; cz <= 2; cz++) {
                        int grid_x = CandidateCell.x + cx;
                        int grid_y = CandidateCell.y + cy;
                        int grid_z = CandidateCell.z + cz;

                        if (
                            grid_x < 0 || grid_x >= Nx ||
                            grid_y < 0 || grid_y >= Ny ||
                            grid_z < 0 || grid_z >= Nz
                        ) continue;

                        int idx = Nx * Ny * grid_z + Nx * grid_y + grid_x;
                        int stored_value = grid[idx];

                        if (stored_value == -1) continue;

                        double dx = point.x - pos_x[stored_value];
                        double dy = point.y - pos_y[stored_value];
                        double dz = point.z - pos_z[stored_value];
                        double dist2 = dx * dx + dy * dy + dz * dz;
                        double r2 = r * r;

                        if (dist2 > r2) {
                            valid = false;
                            break;
                        }
                    }
                    if (!valid) break;
                }
                if (!valid) break;
            }
            if (valid) {
                int point_grid_idx = Nx * Ny * CandidateCell.z + Nx * CandidateCell.y + CandidateCell.x;
                int Particle_ID = pos_x.size();

                pos_x.push_back(point.x);
                pos_y.push_back(point.y);
                pos_z.push_back(point.z);

                grid[point_grid_idx] = Particle_ID;
                active.push_back(Particle_ID);

                point_found = true;
                break;
            }
        }
        if (!point_found) {
            active[activeIndex] = active.back();
            active.pop_back();
        }
    }
}

CellCord GasCloud::toCell(double x, double y, double z) {
    CellCord cell;
    cell.x = std::floor((x + R) / cell_size);
    cell.y = std::floor((y + R) / cell_size);
    cell.z = std::floor((z + R) / cell_size);
    return cell;
}