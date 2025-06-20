# Orbital Simulation with REBOUND

This project simulates the orbital dynamics of particles (e.g., asteroids or small bodies) around a central star using the REBOUND N-body simulation package. It employs a custom velocity Verlet integrator, combined with REBOUND's WHFast integrator, to compute gravitational interactions. This generates a 2D (or optionally 3D) animation of particle orbits with trails, which is saved as an MP4 video using FFmpeg. The simulation includes stability checks for escapes or near-collisions and calculates particle densities. Units are astronomical: AU (astronomical units) for distance, years for time, and solar masses (Msun) for mass.


# Features

Simulates 50 particles orbiting a 1 Msun central star with random orbital elements (semi-major axis, eccentricity, inclination, etc.).

Uses velocity Verlet integration with a small timestep (10⁻⁵ yr) for accuracy.

Checks for unbound particles (energy ≥ 0), escapes (r > 100 AU), or near-collisions (r < 0.01 AU).

Visualizes orbits in 2D (default) or 3D with particle trails (up to 20 frames long).

Outputs particle data (ID, density, mass, semi-major axis, radius) to the console and an MP4 animation.

Configurable parameters via a CONFIG dictionary for easy customization.


# Requirements

Python: 3.8–3.11


# Dependencies:

rebound (3.24.0 recommended): N-body simulation library.

numpy: Numerical computations.

matplotlib: Plotting and animation.

ffmpeg-python: Interface for FFmpeg (optional, for MP4 output).

FFmpeg: Required for saving the animation as an MP4 video. Install via:

Windows: Download from ffmpeg.org and add to PATH.

macOS: brew install ffmpeg

Linux: sudo apt-get install ffmpeg

Hardware: 8-core CPU, 16GB RAM recommended for smooth performance.


# Installation

Clone or download this repository:

git clone https://github.com/LuceilIsRight/Star-System
cd orbital-simulation


# Create a virtual environment and activate it:

python -m venv venv
source venv/bin/activate  # Linux/macOS
venv\Scripts\activate     # Windows


# Install Python dependencies:

pip install rebound==3.24.0 numpy matplotlib ffmpeg-python

# Verify FFmpeg installation:

ffmpeg -version

If not found, install FFmpeg and ensure it’s in your PATH.

# Usage

Save the code as orbital_simulation.py.

# Run the simulation:

python orbital_simulation.py


# The script will:

Initialize the simulation with 50 particles.

Print a table of particle properties (ID, density, mass, semi-major axis, radius).

Generate an MP4 animation named orbit_simulation_YYYYMMDD_HHMMSS.mp4.

View the output video using any media player.


# Configuration

Edit the CONFIG dictionary in orbital_simulation.py to customize the simulation:

central_mass: Mass of the central star (default: 1.0 Msun).

N_particles: Number of particles (default: 50).

a_range: Semi-major axis range in AU (default: 0.5–5.0).

e_range: Eccentricity range (default: 0.0–0.05).

inc_range: Inclination range in radians (default: 0.0–0.1).

time_span: Simulation duration in years (default: 1000.0).

N_outputs: Number of animation frames (default: 300).

mass_range: Particle mass range in Msun (default: 10⁻⁷–10⁻⁴; uncomment for 10⁻¹⁰–10⁻⁷ for realistic ~1–10 g/cm³).

particle_radius_range: Particle radius range in AU (default: 10⁻⁶–5×10⁻⁶ AU, ~150–750 km).

plot_3d: Enable 3D plotting (default: False).

trail_length: Length of particle trails in frames (default: 20).

fps: Frames per second for the animation (default: 15).

dpi: Resolution of the output video (default: 300).

dt: Integration timestep in years (default: 10⁻⁵).

# Example modification for a shorter simulation with fewer particles:

CONFIG = {
    'central_mass': 1.0,
    'N_particles': 20,
    'time_span': 500.0,
    'N_outputs': 150,
    # ... other parameters
}

# Output

Console Output: A table of particle properties, e.g.:

Particle ID | Density (g/cm^3) | Mass (Msun) | Semi-major Axis (AU) | Radius (AU)
---------------------------------------------------------------------------
          1 |         1.23e+04 |   1.00e-07 |                 2.35 |   1.00e-06
...



Animation: An MP4 file (orbit_simulation_YYYYMMDD_HHMMSS.mp4) showing:

Cyan particles orbiting a yellow central star.

Particle trails (cyan, semi-transparent) showing recent paths.

Time and particle count in the title.

2D view (default) or 3D view if plot_3d is True.

Warnings: Console warnings for unbound particles, escapes (r > 100 AU), or near-collisions (r < 0.01 AU).


# Performance

Runtime: ~5–10 minutes on an 8-core CPU with 16GB RAM for default settings.

Memory: ~1–2 GB RAM.

Bottlenecks: Animation rendering and FFmpeg encoding. Reduce N_outputs, dpi, or N_particles for faster execution.

# Troubleshooting

FFmpeg Not Found:

Ensure FFmpeg is installed and in your PATH.

Alternatively, use a different Matplotlib writer (e.g., pillow for GIF):

ani.save('orbit_simulation.gif', writer='pillow', dpi=CONFIG['dpi'])



# REBOUND Errors:

Verify REBOUND version:

python -c "import rebound; print(rebound.__version__)"

Reinstall: pip install rebound==3.24.0

High Density Values:

The default mass_range (10⁻⁷–10⁻⁴ Msun) yields high densities (~10⁴ g/cm³). For realistic asteroid-like densities (~1–10 g/cm³), uncomment:

'mass_range': (1e-10, 1e-7),

# Slow Performance:

Reduce N_particles, N_outputs, or dpi in CONFIG.

Example:

CONFIG['N_particles'] = 20
CONFIG['N_outputs'] = 100
CONFIG['dpi'] = 150

# Notes

The simulation uses a small timestep (10⁻⁵ yr) for accuracy, which may slow down large simulations. Adjust dt cautiously to balance speed and stability.

The WHFast integrator is optimized for near-Keplerian orbits, suitable for the low-eccentricity (e_range: 0.0–0.05) particles in this setup.

For 3D visualization, enable plot_3d in CONFIG, but note increased rendering time.

Random seed (np.random.seed(42)) ensures reproducible results.

The derivations of equations will be posted on ArXiv.com.

# License

This project is unlicensed. Feel free to use, modify, and distribute the code as needed.

# Contact

For issues or questions, please open an issue on the repository or contact the maintainer.

# Authors

@LuceilIsRight(https://github.com/LuceilIsRight)

@K-Suvan(https://github.com/K-Suvan)
