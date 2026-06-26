#include "gravity.hpp"

constexpr double G = 6.67430e-11;
constexpr double EPS = 1e-10;
constexpr double EPS2 = EPS * EPS;

// Node constructor
Node::Node(const vec3& c, const double hs) : center(c), half_side(hs), body_idx(-1), first_child(-1) {
    if (hs <= 0) throw std::invalid_argument("Half Side should be Positive.");
}

// Octree constructor
Octree::Octree(
    const double* masses_,
    const std::vector<double>& pos_x_,
    const std::vector<double>& pos_y_,
    const std::vector<double>& pos_z_,
    std::vector<double>& acc_x_,
    std::vector<double>& acc_y_,
    std::vector<double>& acc_z_,
    const size_t num_particles_,
    const double theta_,
    const std::vector<GasMaterial>& materials_,
    const std::vector<int>& material_ids_
) : masses(masses_), pos_x(pos_x_), pos_y(pos_y_), pos_z(pos_z_),
    acc_x(acc_x_), acc_y(acc_y_), acc_z(acc_z_),
    root_idx(-1),
    num_particles(num_particles_), theta(theta_),
    materials(materials_), material_ids(material_ids_) {
    if (num_particles_ == 0) throw std::invalid_argument("Cannot build octree: positions array is empty.");
}

// getRoot
const Node* Octree::getRoot() const {
    if (root_idx == -1) return nullptr;
    return &pool[root_idx];
}

// build
void Octree::build() {
    pool.clear();
    pool.reserve(num_particles * 4);

    double xmin = pos_x[0], ymin = pos_y[0], zmin = pos_z[0];
    double xmax = pos_x[0], ymax = pos_y[0], zmax = pos_z[0];

    for (size_t i = 1; i < num_particles; i++) {

        if (pos_x[i] < xmin) xmin = pos_x[i];
        if (pos_y[i] < ymin) ymin = pos_y[i];
        if (pos_z[i] < zmin) zmin = pos_z[i];

        if (pos_x[i] > xmax) xmax = pos_x[i];
        if (pos_y[i] > ymax) ymax = pos_y[i];
        if (pos_z[i] > zmax) zmax = pos_z[i];
    }
    
   double half_side = std::max({xmax - xmin, ymax - ymin, zmax - zmin}) * 0.5 + 1e-5;

   double cx = (xmax + xmin) * 0.5, cy = (ymax + ymin) * 0.5, cz = (zmax + zmin) * 0.5;

   vec3 center{cx, cy, cz};
   Node root(center, half_side);
   root_idx = pool.size();
   pool.push_back(root);

   for (size_t i = 0; i < num_particles; i++) insert(i);
}

void Octree::insert(int i) {
    if (root_idx == -1) throw std::invalid_argument("Call build() before insert().");

    insert_impl(root_idx, i);
}

int Octree::allocate_8_children(int node_idx) {
    vec3 node_center = pool[node_idx].center;
    double child_half = pool[node_idx].half_side * 0.5;

    int first_child_idx = pool.size();

    for (double dx : {-child_half, child_half}) {
        for (double dy : {-child_half, child_half}) {
            for (double dz : {-child_half, child_half}) {

                vec3 child_center = {
                    node_center.x + dx,
                    node_center.y + dy,
                    node_center.z + dz
                };

                Node child(child_center, child_half);

                pool.push_back(child);
            }
        }
    }

    return first_child_idx;
}

int Octree::choose_child(int node_idx, double x, double y, double z) {
    int octant = get_octant_idx(node_idx, x, y, z);
    return pool[node_idx].first_child + octant;
}

int Octree::get_octant_idx(int node_idx, double x, double y, double z) {
    int ix = (x >= pool[node_idx].center.x) ? 1 : 0;
    int iy = (y >= pool[node_idx].center.y) ? 1 : 0;
    int iz = (z >= pool[node_idx].center.z) ? 1 : 0;
    return ix * 4 + iy * 2 + iz;
}

double Octree::point_to_node_distance(int node_idx, double x, double y, double z) {
    double dx, dy, dz;

    Node& node = pool[node_idx];

    if (x < (node.center.x - node.half_side)) dx = (node.center.x - node.half_side) - x;
    else if (x > (node.center.x + node.half_side)) dx = x - (node.center.x + node.half_side);
    else dx = 0.0;

    if (y < (node.center.y - node.half_side)) dy = (node.center.y - node.half_side) - y;
    else if (y > (node.center.y + node.half_side)) dy = y - (node.center.y + node.half_side);
    else dy = 0.0;

    if (z < (node.center.z - node.half_side)) dz = (node.center.z - node.half_side) - z;
    else if (z > (node.center.z + node.half_side)) dz = z - (node.center.z + node.half_side);
    else dz = 0.0;

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void Octree::insert_impl(int node_idx, int i) {
    Node& node = pool[node_idx];
    double old_mass = node.mass;
    double new_mass = old_mass + masses[i];
    vec3 r = {pos_x[i], pos_y[i], pos_z[i]};

    if (old_mass == 0.0) {
        node.com = r;
    } else {
        node.com = (node.com * old_mass + r * masses[i]) / new_mass;
    }
    node.mass = new_mass;

    // case one: Empty Leaf
    if (node.body_idx == -1 && node.first_child == -1) {
        node.body_idx = i;
        return;
    }

    // case two: Leaf with One Body
    else if (node.body_idx != -1 && node.first_child == -1) {
        if (node.half_side < 1e-5) return;

        int old_body_idx = node.body_idx;
        node.body_idx = -1;
        pool[node_idx].first_child = allocate_8_children(node_idx);
        insert_impl(choose_child(node_idx, pos_x[old_body_idx], pos_y[old_body_idx], pos_z[old_body_idx]), old_body_idx);
        insert_impl(choose_child(node_idx, pos_x[i], pos_y[i], pos_z[i]), i);
        return;
    }

    // case three: Internal Node
    else {
        insert_impl(choose_child(node_idx, pos_x[i], pos_y[i], pos_z[i]), i);
    }
}

void Octree::compute_node_forces(int node_idx, int i) {
    const GasMaterial& mat_i = materials[material_ids[i]];
    double epsilon = mat_i.kernel_radius * 0.1;
    double eps2 = epsilon * epsilon;

    Node& node = pool[node_idx];

    if (node.mass < EPS) return;

    if (node.first_child == -1 && node.body_idx == i) return;

    if (node.first_child == -1 && node.body_idx != i) {
        int j = node.body_idx;
        double dx = pos_x[j] - pos_x[i], dy = pos_y[j] - pos_y[i], dz = pos_z[j] - pos_z[i];
        double inv_dist = 1.0 / std::sqrt(dx*dx + dy*dy + dz*dz + eps2);
        double inv_dist3 = inv_dist * inv_dist * inv_dist;
        double f_common = G * masses[j] * inv_dist3;

        acc_x[i] += f_common * dx,
        acc_y[i] += f_common * dy,
        acc_z[i] += f_common * dz;
        return;
    }

    double s = node.half_side * 2.0;
    double rx = node.com.x - pos_x[i], ry = node.com.y - pos_y[i], rz = node.com.z - pos_z[i];
    double d = std::sqrt(rx*rx + ry*ry + rz*rz + eps2);
    if (s / d < theta) {
        double inv_d = 1.0 / d;
        double inv_d3 = inv_d * inv_d * inv_d;
        double factor = G * node.mass * inv_d3;

        acc_x[i] += factor * rx,
        acc_y[i] += factor * ry,
        acc_z[i] += factor * rz;
        return;
    }

    else {
        for (int c = 0; c < 8; c++) {
            int child_idx = pool[node_idx].first_child + c;

            if (pool[child_idx].mass > 0.0) {
                compute_node_forces(child_idx, i);
            }
        }
    }
}

void Octree::computeForces(const std::vector<size_t>& activeParticles) {
    #pragma omp parallel for
    for (size_t k = 0; k < activeParticles.size(); k++) {
        int i = activeParticles[k];
        acc_x[i] = 0.0;
        acc_y[i] = 0.0;
        acc_z[i] = 0.0;

        compute_node_forces(root_idx, i);
    }
}