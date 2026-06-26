#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include "vec3.hpp"
#include "gas_cloud.hpp"
#include "sph.hpp"
#include "gravity.hpp"
#include "gas_cloud.hpp"


class Universe {
    public:
        std::vector<GasMaterial> materials;

        std::vector<double> pos_x;
        std::vector<double> pos_y;
        std::vector<double> pos_z;

        std::vector<double> vel_x;
        std::vector<double> vel_y;
        std::vector<double> vel_z;

        std::vector<double> acc_x;
        std::vector<double> acc_y;
        std::vector<double> acc_z;

        std::vector<double> ext_acc_x;
        std::vector<double> ext_acc_y;
        std::vector<double> ext_acc_z;

        std::vector<double> density;
        std::vector<double> pressure;

        std::vector<double> u;
        std::vector<double> du_dt;

        std::vector<double> masses;

        std::vector<int> material_ids;

        std::vector<size_t> timeBins;
        std::vector<size_t> activeParticles;

        std::unique_ptr<SPHSolver> solver;
        std::unique_ptr<Octree> gravity;

        double dt;
        double global_theta;
        double damping_factor = 1.0;
        bool enable_vtk_export = false;
        int vtk_export_inteval = 10;

        Universe(double dt, double theta = 0.5);

        int registerMaterial(double kernel_radius, double viscosity, double gamma, double molar_mass);

        void addGasCloud(size_t num_prticles, double cloud_mass, double radius, vec3 center, vec3 initial_vel, double initial_temp, int material_id, double relaxation_time);
        
        void step(double dt);
        
        void run(double T);

        void exportVTK(const std::string& filename);
        
    private:
        void generateRelaxedGasCloud(size_t num_particles, double cloud_mass, double radius, vec3 center, vec3 initial_vel, double initial_temp, int material_id, double relaxation_time);

        void updateTimeBins(int max_bins, size_t current_ti);
};