#pragma once
#include <vector>
#include <memory>
#include <array>
#include <stdexcept>
#include "vec3.hpp"
#include <algorithm>
#include <cmath>
#include "gas_material.hpp"


struct Node {
    vec3 center;
    double half_side;
    int body_idx;
    int first_child;
    vec3 com = vec3::zero;
    double mass = 0.0;

    // Constructor declaration
    Node(const vec3& c, const double hs);
};

class Octree {
    private:
        const double* masses;
        const std::vector<double>& pos_x;
        const std::vector<double>& pos_y;
        const std::vector<double>& pos_z;
        std::vector<double>& acc_x;
        std::vector<double>& acc_y;
        std::vector<double>& acc_z;
        std::vector<Node> pool;
        int root_idx;
        size_t num_particles;
        double theta;
        const std::vector<GasMaterial>& materials;
        const std::vector<int>& material_ids;

        void insert_impl(int node_idx, int i);
        void compute_node_forces(int node_idx, int i);

    public:
        // Constructor declaration
        Octree(
            const double* masses_,
            const std::vector<double>& pos_x_,
            const std::vector<double>& pos_y_,
            const std::vector<double>& pos_z_,
            std::vector<double>& acc_x_,
            std::vector<double>& acc_y_,
            std::vector<double>& acc_z_,
            const size_t num_particles,
            const double theta_,
            const std::vector<GasMaterial>& materials_,
            const std::vector<int>& material_ids_
        );

        const Node* getRoot() const;

        void build();
        void insert(int i);
        int allocate_8_children(int node_idx);
        int choose_child(int node_idx, double x, double y, double z);
        int get_octant_idx(int node_idx, double x, double y, double z);
        double point_to_node_distance(int node_idx, double x, double y, double z);
        void computeForces(const std::vector<size_t>& activeParticles);
};