#pragma once
#include <vector>
#include "vec3.hpp"
#include "spatial_hash.hpp"
#include "gas_material.hpp"


class SPHSolver {
    public:
        SPHSolver(
            const std::vector<GasMaterial>& materials_,
            const std::vector<int>& material_ids_,
            const std::vector<double>& masses_,
            const std::vector<size_t>& activeParticles_,
            std::vector<double>& pos_x_, std::vector<double>& pos_y_, std::vector<double>& pos_z_,
            std::vector<double>& vel_x_, std::vector<double>& vel_y_, std::vector<double>& vel_z_,
            std::vector<double>& acc_x_, std::vector<double>& acc_y_, std::vector<double>& acc_z_,
            std::vector<double>& ext_acc_x_, std::vector<double>& ext_acc_y_, std::vector<double>& ext_acc_z_,
            std::vector<double>& density_, std::vector<double>& pressure_,
            std::vector<double>& u_, std::vector<double>& du_dt_
        );

        void computeDensity(SpatialHash& grid);
        void computePressure();
        void computeAccelerations(SpatialHash& grid);
        void step(double dt, double damping_factor);

    private:
        const std::vector<GasMaterial>& materials;
        const std::vector<int>& material_ids;
        const std::vector<double>& masses;
        const std::vector<size_t>& activeParticles;

        std::vector<double>& pos_x; std::vector<double>& pos_y; std::vector<double>& pos_z;
        std::vector<double>& vel_x; std::vector<double>& vel_y; std::vector<double>& vel_z;
        std::vector<double>& acc_x; std::vector<double>& acc_y; std::vector<double>& acc_z;
        std::vector<double>& ext_acc_x; std::vector<double>& ext_acc_y; std::vector<double>& ext_acc_z;

        std::vector<double>& density;
        std::vector<double>& pressure;

        std::vector<double>& u;
        std::vector<double>& du_dt;
};