#pragma once
#include <vector>
#include <stdexcept>
#include "vec3.hpp"


struct CellCord {
    int x, y, z;

    bool operator==(const CellCord& other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    bool operator!=(const CellCord& other) const {
        return !(*this == other);
    }
};

class SpatialHash {
    private:
        double kernel_radius;
        size_t num_bodies;
        size_t table_size;
        std::vector<int> head;
        std::vector<int> next;
        
    public:
        SpatialHash(double kernel_radius_, size_t num_bodies_)
        : kernel_radius(kernel_radius_), num_bodies(num_bodies_) {
            if (kernel_radius <= 0.0) throw std::invalid_argument("Invalid Radius: kernel_radius should be Positive.");
        };

        void build(const std::vector<double>& pos_x, const std::vector<double>& pos_y, const std::vector<double>& pos_z);

        CellCord toCell(double x, double y, double z);

        size_t hash(const CellCord& cell);

        int getFirst(const CellCord& cell);

        int getNext(int particle_id);
};