# AstroFluid

**AstroFluid** is a high-performance 3D astrophysical gas cloud simulator written in modern C++. It combines **Smoothed Particle Hydrodynamics (SPH)** with an **Octree-based Barnes–Hut gravity solver** to simulate the evolution of self-gravitating gaseous systems.

The project is built with a strong focus on **scientific computing**, **numerical methods**, and **high-performance computing (HPC)**. It is currently optimized for multi-core CPU execution using **OpenMP**, **SIMD vectorization**, and a **Structure of Arrays (SoA)** memory layout. The long-term goal is to scale the simulation to GPU accelerators and distributed HPC clusters.

---

## Features

* **Smoothed Particle Hydrodynamics (SPH)**

  * Gas dynamics using the Ideal Gas Equation of State.
  * Density, pressure, and acceleration calculations.

* **Barnes–Hut Gravity**

  * Efficient N-body gravity using a dynamic octree.
  * Reduces computational complexity from **O(N²)** to **O(N log N)**.

* **Hierarchical Block Time-Stepping**

  * Adaptive time integration based on the CFL condition and particle acceleration.
  * Improves stability while avoiding unnecessary computation.

* **Spatial Hashing**

  * Uniform grid spatial hashing for efficient neighbor searches.
  * Average **O(1)** particle lookup.

* **Kick–Drift–Kick (Leapfrog) Integration**

  * Symplectic integration scheme for improved long-term energy conservation.

* **Poisson Disk Sampling**

  * Generates well-distributed initial particle configurations while avoiding particle overlap.

* **OpenMP Parallelization**

  * Parallel density, pressure, gravity, and particle update loops across CPU cores.

* **SIMD Vectorization**

  * Data-oriented memory layout (Structure of Arrays) designed for compiler auto-vectorization.

* **VTK Output**

  * Native export to `.vtk` format for visualization with ParaView or VisIt.

---

## Algorithms

The simulator currently implements:

* Smoothed Particle Hydrodynamics (SPH)
* Barnes–Hut Octree Gravity
* Hierarchical Block Time-Stepping
* Spatial Hashing
* Poisson Disk Sampling
* Kick–Drift–Kick (Leapfrog) Integration
* OpenMP Parallelization
* SIMD Vectorization

---

## Project Structure

```text
AstroFluid/
├── CMakeLists.txt
├── LICENSE
├── README.md
├── main.cpp
├── src/
│   ├── gas_cloud.cpp
│   ├── gas_cloud.hpp
│   ├── gas_material.hpp
│   ├── gravity.cpp
│   ├── gravity.hpp
│   ├── spatial_hash.cpp
│   ├── spatial_hash.hpp
│   ├── sph.cpp
│   ├── sph.hpp
│   ├── universe.cpp
│   ├── universe.hpp
│   ├── vec3.cpp
│   └── vec3.hpp
└── images/
```

---

## Example Output

Simulation frames are exported as VTK files and can be visualized using **ParaView** or **VisIt**.

> *(Add screenshots or GIFs here once available.)*

---

## Requirements

* C++20 compatible compiler

  * GCC
  * Clang
  * MSVC

* CMake 3.15 or newer

* OpenMP

---

## Building

Clone the repository:

```bash
git clone https://github.com/yourusername/AstroFluid.git
cd AstroFluid
```

Configure the project:

```bash
cmake -B build
```

Build:

```bash
cmake --build build -j
```

Alternatively:

```bash
mkdir build
cd build
cmake ..
make -j
```

---

## Running

After building:

```bash
./build/simulation
```

Simulation parameters such as particle count, cloud radius, temperature, mass, and output frequency can be configured in `main.cpp`.

Generated VTK files can be opened directly in **ParaView** for visualization.

---

## Performance

Current CPU optimizations include:

* OpenMP multi-threading
* SIMD-friendly Structure of Arrays (SoA)
* Barnes–Hut gravity
* Spatial hashing for neighbor searches
* Adaptive hierarchical time-stepping

---

## Roadmap

### In Progress

* Improved SPH validation
* Additional benchmark cases
* Performance profiling
* Code refactoring

### Planned

* CUDA GPU acceleration
* MPI domain decomposition
* Distributed-memory execution
* Additional equations of state
* Individual adaptive particle time-stepping
* HDF5 checkpointing and restart support
* Continuous benchmarking

---

## Motivation

AstroFluid is a personal project created to explore the intersection of:

* Scientific Computing
* High-Performance Computing (HPC)
* Astrophysical Simulation
* Computational Physics
* Numerical Methods
* Modern C++

The project serves as a platform for learning and implementing professional techniques used in large-scale particle simulations.

---

## License

This project is licensed under the MIT License.

See the **LICENSE** file for details.
