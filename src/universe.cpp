#include "universe.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>


Universe::Universe(double dt_, double theta) : dt(dt_), global_theta(theta) {};

int Universe::registerMaterial(double kernel_radius, double viscosity, double gamma, double molar_mass) {
    GasMaterial new_mat = {kernel_radius, viscosity, gamma, molar_mass};
    materials.push_back(new_mat);

    return materials.size() - 1;
}

void Universe::generateRelaxedGasCloud(size_t num_particles, double cloud_mass, double radius, vec3 center, vec3 initial_vel, double initial_temp, int material_id, double relaxation_time) {

    std::cout << "Engine: Generating Raw lattice..." << std::endl;

    // Initialize a temporary universe to relax the gas cloud
    Universe temp_universe(this->dt, this->global_theta);
    // temp_universe.damping_factor = 0.95;

    // Initialize gas cloud at (0, 0, 0) center with no velocity
    // This will arrange masses, positions and velocities
    int k = 30;
    GasCloud raw_cloud(
        num_particles, cloud_mass, radius, temp_universe.masses, k,
        temp_universe.pos_x, temp_universe.pos_y, temp_universe.pos_z,
        temp_universe.vel_x, temp_universe.vel_y, temp_universe.vel_z
    );

    // Get the actual number of particles
    size_t actual_particles = temp_universe.pos_x.size();

    // Register materialsin temporary universe
    const GasMaterial& mat = this->materials[material_id];
    int temp_mat_id = temp_universe.registerMaterial(mat.kernel_radius, mat.viscosity, mat.gamma, mat.molar_mass);

    double cv = 8.314 / (mat.molar_mass * (mat.gamma - 1.0));
    double u0 = cv * initial_temp;

    temp_universe.u.assign(actual_particles, u0);
    temp_universe.du_dt.assign(actual_particles, 0.0);
    temp_universe.material_ids.assign(actual_particles, temp_mat_id);

    temp_universe.timeBins.assign(actual_particles, 0);

    // resize physical arrays for the temporary universe
    temp_universe.acc_x.resize(actual_particles, 0.0);
    temp_universe.acc_y.resize(actual_particles, 0.0);
    temp_universe.acc_z.resize(actual_particles, 0.0);
    temp_universe.ext_acc_x.resize(actual_particles, 0.0);
    temp_universe.ext_acc_y.resize(actual_particles, 0.0);
    temp_universe.ext_acc_z.resize(actual_particles, 0.0);
    temp_universe.density.resize(actual_particles, 0.0);
    temp_universe.pressure.resize(actual_particles, 0.0);

    // run: the relaxation phase
    std::cout << "Relaxing Cloud for " << relaxation_time << " seconds..." << std::endl;
    temp_universe.run(relaxation_time);

    // port the engine particle into the main universe
    std::cout << "Engine: Stamping relaxed cloud into main simulation..." << std::endl;

    for (size_t i = 0; i < actual_particles; i++) {
        this->pos_x.push_back(temp_universe.pos_x[i] + center.x);
        this->pos_y.push_back(temp_universe.pos_y[i] + center.y);
        this->pos_z.push_back(temp_universe.pos_z[i] + center.z);

        this->vel_x.push_back(initial_vel.x);
        this->vel_y.push_back(initial_vel.y);
        this->vel_z.push_back(initial_vel.z);

        this->masses.push_back(temp_universe.masses[i]);
        this->material_ids.push_back(material_id);

        this->u.push_back(u0);
        this->du_dt.push_back(0.0);

        this->timeBins.push_back(0);
    }

    // resize the main universe physical arrays
    size_t new_total = this->pos_x.size();
    this->acc_x.resize(new_total, 0.0);
    this->acc_y.resize(new_total, 0.0);
    this->acc_z.resize(new_total, 0.0);

    this->ext_acc_x.resize(new_total, 0.0);
    this->ext_acc_y.resize(new_total, 0.0);
    this->ext_acc_z.resize(new_total, 0.0);

    this->density.resize(new_total, 0.0);
    this->pressure.resize(new_total, 0.0);
}

void Universe::addGasCloud(size_t num_particles, double cloud_mass, double radius, vec3 center, vec3 initial_vel, double initial_temp, int material_id, double relaxation_time) {
    generateRelaxedGasCloud(num_particles, cloud_mass, radius, center, initial_vel, initial_temp, material_id, relaxation_time);
    }

void Universe::step(double dt) {
    size_t total_particles = pos_x.size();
    if (total_particles == 0) return;

    // Kick 1 and Drift
    #pragma omp parallel for simd
    for (size_t k = 0; k < activeParticles.size(); k++) {
        int i = activeParticles[k];
        // Kick 1
        vel_x[i] += 0.5 * (ext_acc_x[i] + acc_x[i]) * dt;
        vel_y[i] += 0.5 * (ext_acc_y[i] + acc_y[i]) * dt;
        vel_z[i] += 0.5 * (ext_acc_z[i] + acc_z[i]) * dt;

        vel_x[i] *= damping_factor;
        vel_y[i] *= damping_factor;
        vel_z[i] *= damping_factor;

        u[i] += 0.5 * du_dt[i] * dt;
        const GasMaterial& mat = materials[material_ids[i]];
        double cv = 8.314 / (mat.molar_mass * (mat.gamma - 1.0));
        double u_min = cv * 2.7;
        u[i] = std::max(u[i], u_min);

        // Drift
        pos_x[i] += vel_x[i] * dt;
        pos_y[i] += vel_y[i] * dt;
        pos_z[i] += vel_z[i] * dt;

        // clear acceleration to prepare for the new force calculation
        acc_x[i] = 0.0;
        acc_y[i] = 0.0;
        acc_z[i] = 0.0;
        ext_acc_x[i] = 0.0;
        ext_acc_y[i] = 0.0;
        ext_acc_z[i] = 0.0;

        du_dt[i] = 0.0;
    }

    if (!solver) {
        solver = std::make_unique<SPHSolver>(
            materials, material_ids, masses,
            activeParticles,
            pos_x, pos_y, pos_z,
            vel_x, vel_y, vel_z,
            acc_x, acc_y, acc_z,
            ext_acc_x, ext_acc_y, ext_acc_z,
            density, pressure, u, du_dt
        );
    }

    gravity = std::make_unique<Octree>(
        masses.data(), pos_x, pos_y, pos_z,
        ext_acc_x, ext_acc_y, ext_acc_z, total_particles, global_theta, materials, material_ids
    );
    gravity->build();
    gravity->computeForces(activeParticles);
    solver->step(dt, damping_factor);
}

void Universe::updateTimeBins(int max_bins, size_t current_ti) {
    constexpr double C_CFL = 0.2; // Courant safety factor
    constexpr double ETA = 0.15; // Acceleration safety factor

    #pragma omp parallel for
    for (size_t k = 0; k < activeParticles.size(); k++) {
        int i = activeParticles[k];
        const GasMaterial& mat = materials[material_ids[i]];

        // 1. Hydrodynamic Timestep (CFL condition)
        double c_s = std::sqrt(mat.gamma * (mat.gamma - 1.0) * u[i]); // NOTE: speed of sound in gas
        double dt_hydro = C_CFL * (mat.kernel_radius / (c_s + 1e-8));

        // 2. Kinamatic Timestep (Acceleration condition)
        double ax = acc_x[i] + ext_acc_x[i];
        double ay = acc_y[i] + ext_acc_y[i];
        double az = acc_z[i] + ext_acc_z[i];
        double a_mag = std::sqrt(ax*ax + ay*ay + az*az);
        double dt_acc = ETA * std::sqrt(mat.kernel_radius / (a_mag + 1e-8));

        // 3. Ideal Timestep
        double dt_ideal = std::min(dt_hydro , dt_acc);

        int target_bin = 0;
        double dt_block = this->dt;
        while (dt_block > dt_ideal && target_bin < max_bins) {
            dt_block *= 0.5;
            target_bin++;
        }

        int current_bin = timeBins[i];

        if (target_bin > current_bin) {
            timeBins[i] = target_bin;
        } else if (target_bin < current_bin) {
            size_t target_ti_step = 1 << (max_bins - target_bin);

            if (current_ti % target_ti_step == 0) timeBins[i] = target_bin;
        }
    }
}

void Universe::run(double T) {
    constexpr int MAX_BINS = 10;
    constexpr size_t TIME_BASE = 1 << MAX_BINS;
    size_t current_ti = 0;
    double current_time = 0.0;
    int frame = 0;
    int sync_step_count = 0;

    while (current_time < T) {
        int min_bin_in_system = 0;
        for (int bin : timeBins) {
            if (bin > min_bin_in_system) min_bin_in_system = bin;
        }
        
        size_t ti_step = 1 << (MAX_BINS - min_bin_in_system);
        
        current_ti += ti_step;
        
        double physical_dt_step = dt * (static_cast<double>(ti_step) / TIME_BASE);
        current_time += physical_dt_step;
        
        activeParticles.clear();
        
        for (size_t i = 0; i < pos_x.size(); i++) {
            int particle_bin = timeBins[i];
            size_t particle_ti_step = 1 << (MAX_BINS - particle_bin);
            if (current_ti % particle_ti_step == 0) activeParticles.push_back(i);
        }
        
        step(physical_dt_step);
        updateTimeBins(MAX_BINS, current_ti);

        if (current_ti % TIME_BASE == 0) {
            sync_step_count++;
            if (enable_vtk_export && (sync_step_count % vtk_export_inteval == 0)) {
                std::stringstream ss;
                ss << "sim_" << std::setfill('0') << std::setw(4) << frame++ << ".vtk";
                exportVTK(ss.str());
                std::cout << "Engine: Exported " << ss.str() << " at t = " << current_time << std::endl;
            }
        }
    }
}

void Universe::exportVTK(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to Open " << filename << " for writing.\n";
        return;
    }
    
    size_t n = pos_x.size();
    
    // VTK Header
    file << "# vtk DataFile Version 3.0\n";
    file << "SPH Simulation Data\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    // Positions
    file << "POINTS " << n << " float\n";
    for (size_t i = 0; i < n; i++) {
        file << pos_x[i] << " " << pos_y[i] << " " << pos_z[i] << "\n";
    }

    file << "POINT_DATA " << n << "\n";

    // Density
    file << "SCALARS density float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < n; i++) file << density[i] << "\n";

    // Pressure
    file << "SCALARS pressure float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < n; i++) file << pressure[i] << "\n";

    // Internal Energy
    file << "SCALARS internal_energy float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < n; i++) file << u[i] << "\n";

    // Velocities
    file << "VECTORS velocity float\n";
    for (size_t i = 0; i < n; i++) {
        file << vel_x[i] << " " << vel_y[i] << " " << vel_z[i] << "\n";
    }

    file.close();
}