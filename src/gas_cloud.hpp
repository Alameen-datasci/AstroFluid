#pragma once
#include <vector>
#include "sph.hpp"
#include "vec3.hpp"


class GasCloud {
    public:
        GasCloud(
            size_t num_particles, double total_mass, double radius, std::vector<double>& masses_,
            int k_,
            std::vector<double>& pos_x_, std::vector<double>& pos_y_, std::vector<double>& pos_z_,
            std::vector<double>& vel_x_, std::vector<double>& vel_y_, std::vector<double>& vel_z_
        );

    private:
        void generateParticle();

        CellCord toCell(double x, double y, double z);

        size_t num_particles;
        double total_mass;
        double R;
        double r;
        int k;
        double cell_size;
        int Nx, Ny, Nz;

        std::vector<double>& masses;

        std::vector<double>& pos_x;
        std::vector<double>& pos_y;
        std::vector<double>& pos_z;

        std::vector<double>& vel_x;
        std::vector<double>& vel_y;
        std::vector<double>& vel_z;

        std::vector<int> active;
        std::vector<int> grid;
};