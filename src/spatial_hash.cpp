#include <cmath>
#include <algorithm>
#include "spatial_hash.hpp"


void SpatialHash::build(const std::vector<double>& pos_x, const std::vector<double>& pos_y, const std::vector<double>& pos_z) {
    table_size = num_bodies * 10;

    if (head.size() != table_size) head.resize(table_size, -1);
    if (next.size() != num_bodies) next.resize(num_bodies, -1);

    std::fill(head.begin(), head.end(), -1);
    std::fill(next.begin(), next.end(), -1);

    for (size_t i = 0; i < num_bodies; i++) {
        int cell_hash = hash(toCell(pos_x[i], pos_y[i], pos_z[i]));

        next[i] = head[cell_hash];

        head[cell_hash] = i;
    }
}

CellCord SpatialHash::toCell(double x, double y, double z) {
    CellCord cell;
    cell.x = std::floor(x / kernel_radius);
    cell.y = std::floor(y / kernel_radius);
    cell.z = std::floor(z / kernel_radius);
    return cell;
}

size_t SpatialHash::hash(const CellCord& cell) {
    size_t hx = static_cast<size_t>(cell.x) * 73856093;
    size_t hy = static_cast<size_t>(cell.y) * 19349663;
    size_t hz = static_cast<size_t>(cell.z) * 83492791;
    return (hx ^ hy ^ hz) % table_size;
}

int SpatialHash::getFirst(const CellCord& cell) {
    return head[hash(cell)];
}

int SpatialHash::getNext(int particle_id) {
    return next[particle_id];
}