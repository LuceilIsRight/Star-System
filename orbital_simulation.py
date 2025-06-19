import rebound
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from datetime import datetime
import warnings
import shutil

# Configuration dictionary
CONFIG = {
    'central_mass': 1.0,  # Msun
    'N_particles': 50,
    'a_range': (0.5, 5.0),  # AU
    'e_range': (0.0, 0.05),  # Tighter eccentricity
    'inc_range': (0.0, 0.1),  # radians
    'time_span': 1000.0,  # years
    'N_outputs': 300,
    'mass_range': (1e-7, 1e-4),  # Current range (high density)
    # 'mass_range': (1e-10, 1e-7),  # Uncomment for realistic density (~1–10 g/cm³)
    'particle_radius_range': (1e-6, 5e-6),  # AU (~150–750 km)
    'plot_3d': False,
    'trail_length': 20,
    'fps': 15,
    'dpi': 300,
    'dt': 1e-5,  # Smaller timestep
}


def initialize_simulation():
    """Initialize REBOUND simulation with a central star and particles."""
    sim = rebound.Simulation()
    sim.units = ('AU', 'yr', 'Msun')
    sim.add(m=CONFIG['central_mass'])

    np.random.seed(42)

    radius_array = np.linspace(CONFIG['particle_radius_range'][0],
                               CONFIG['particle_radius_range'][1],
                               CONFIG['N_particles'])
    mass_array = np.logspace(np.log10(CONFIG['mass_range'][0]),
                             np.log10(CONFIG['mass_range'][1]),
                             CONFIG['N_particles'])

    G = 4 * np.pi ** 2
    particle_data = []

    for i in range(CONFIG['N_particles']):
        a = np.random.uniform(*CONFIG['a_range'])
        e = np.random.uniform(*CONFIG['e_range'])
        inc = np.random.uniform(*CONFIG['inc_range'])
        Omega = np.random.uniform(0, 2 * np.pi)
        omega = np.random.uniform(0, 2 * np.pi)
        f = np.random.uniform(0, 2 * np.pi)
        mass_msun = mass_array[i]
        radius_au = radius_array[i]

        radius_cm = radius_au * 1.496e13
        volume_cm3 = (4 / 3) * np.pi * radius_cm ** 3
        mass_g = mass_msun * 1.989e33
        density = mass_g / volume_cm3

        sim.add(m=mass_msun, a=a, e=e, inc=inc, Omega=Omega, omega=omega, f=f)

        p = sim.particles[-1]
        r = np.sqrt(p.x ** 2 + p.y ** 2 + p.z ** 2)
        v2 = p.vx ** 2 + p.vy ** 2 + p.vz ** 2
        energy = 0.5 * v2 - G * CONFIG['central_mass'] / r
        if energy >= 0:
            warnings.warn(f"Particle {i + 1} unbound: energy={energy:.2e}, a={a:.2f}, e={e:.3f}")
            sim.remove(sim.N - 1)
            continue

        particle_data.append((i + 1, density, mass_msun, a, radius_au))

    sim.move_to_com()
    sim.t = 0.0
    return sim, particle_data


def get_accelerations(sim):
    """Compute accelerations for all particles."""
    sim.integrator = "whfast"
    sim.dt = 1e-6
    sim.integrate(sim.t + 1e-6, exact_finish_time=0)
    acc = np.zeros((sim.N, 3))
    for i in range(sim.N):
        acc[i, 0] = sim.particles[i].ax
        acc[i, 1] = sim.particles[i].ay
        acc[i, 2] = sim.particles[i].az
    return acc


def verlet_step(sim, dt):
    """Perform one velocity Verlet step."""
    N = sim.N
    pos = np.array([[p.x, p.y, p.z] for p in sim.particles])
    vel = np.array([[p.vx, p.vy, p.vz] for p in sim.particles])
    acc = get_accelerations(sim)

    new_pos = pos + vel * dt + 0.5 * acc * dt ** 2
    for i in range(N):
        sim.particles[i].x, sim.particles[i].y, sim.particles[i].z = new_pos[i]

    new_acc = get_accelerations(sim)
    new_vel = vel + 0.5 * (acc + new_acc) * dt

    for i in range(N):
        sim.particles[i].vx, sim.particles[i].vy, sim.particles[i].vz = new_vel[i]
        sim.particles[i].ax, sim.particles[i].ay, sim.particles[i].az = new_acc[i]

    sim.t += dt


def integrate_verlet(sim, target_time):
    """Integrate to target_time using velocity Verlet."""
    dt = CONFIG['dt']
    current_time = sim.t
    if current_time > target_time:
        sim.t = 0.0
        current_time = 0.0
    steps = int(np.ceil((target_time - current_time) / dt))
    for _ in range(steps):
        if current_time + dt > target_time:
            dt_adjusted = target_time - current_time
            verlet_step(sim, dt_adjusted)
            break
        verlet_step(sim, dt)
        current_time += dt


def check_simulation_stability(sim):
    """Check for escapes or collisions, return indices to remove."""
    to_remove = []
    for i in range(1, sim.N):
        p = sim.particles[i]
        r = np.sqrt(p.x ** 2 + p.y ** 2 + p.z ** 2)
        if r > 100:
            warnings.warn(f"Particle {i} escaped: r={r:.2f} AU, t={sim.t:.2f} yr, a={p.a:.2f}, e={p.e:.3f}")
            to_remove.append(i)
        elif r < 0.01:
            warnings.warn(f"Particle {i} near collision: r={r:.2f} AU, t={sim.t:.2f} yr")
    return to_remove


def setup_plot():
    """Set up Matplotlib figure and axes."""
    plt.style.use('dark_background')
    fig = plt.figure(figsize=(6, 6))
    if CONFIG['plot_3d']:
        from mpl_toolkits.mplot3d import Axes3D
        ax = fig.add_subplot(111, projection='3d')
        ax.set_zlabel("z [AU]", color='white')
        ax.set_zlim(-6, 6)
    else:
        ax = fig.add_subplot(111)
    ax.set_facecolor('black')
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_xlabel("x [AU]", color='white')
    ax.set_ylabel("y [AU]", color='white')
    ax.tick_params(colors='white')
    ax.set_title(f"Orbits Around {CONFIG['central_mass']} Msun Star", color='white')
    return fig, ax


def update(frame, sim, ax, trails):
    """Update animation frame."""
    time = frame
    integrate_verlet(sim, time)
    to_remove = check_simulation_stability(sim)

    for idx in sorted(to_remove, reverse=True):
        if idx < sim.N:
            sim.remove(idx)

    ps = sim.particles
    x = [p.x for p in ps[1:]]
    y = [p.y for p in ps[1:]]
    z = [p.z for p in ps[1:]] if CONFIG['plot_3d'] else None

    trails.append((x, y, z))
    if len(trails) > CONFIG['trail_length']:
        trails.pop(0)

    ax.clear()
    ax.set_facecolor('black')

    for trail_x, trail_y, trail_z in trails:
        if CONFIG['plot_3d']:
            ax.plot(trail_x, trail_y, trail_z, 'cyan', alpha=0.3, linewidth=0.5)
        else:
            ax.plot(trail_x, trail_y, 'cyan', alpha=0.3, linewidth=1)

    if CONFIG['plot_3d']:
        ax.scatter(x, y, z, c='cyan', s=1)
        ax.scatter([0], [0], [0], c='yellow', s=60)
        ax.set_zlim(-6, 6)
        ax.set_zlabel("z [AU]", color='white')
    else:
        ax.scatter(x, y, c='cyan', s=2)
        ax.plot([0], [0], 'yo', markersize=10)

    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_xlabel("x [AU]", color='white')
    ax.set_ylabel("y [AU]", color='white')
    ax.tick_params(colors='white')
    ax.set_title(f"Time = {time:.2f} yr, {sim.N - 1} Particles", color='white')
    return [ax]


def main():
    """Run the orbital simulation and animation."""
    sim, particle_data = initialize_simulation()

    print("\nParticle Density, Mass, Orbital Radius, and Size:")
    print("Particle ID | Density (g/cm^3) | Mass (Msun) | Semi-major Axis (AU) | Radius (AU)")
    print("-" * 75)
    for pid, density, mass_msun, a, radius_au in particle_data:
        print(f"{pid:11d} | {density:15.2f} | {mass_msun:.2e} | {a:20.2f} | {radius_au:.2e}")

    fig, ax = setup_plot()
    times = np.linspace(0, CONFIG['time_span'], CONFIG['N_outputs'])
    trails = []

    try:
        if shutil.which('ffmpeg') is None:
            raise RuntimeError("FFmpeg not found. Install FFmpeg or use a different writer (e.g., pillow for GIF).")
        interval = 1000 / CONFIG['fps']
        ani = FuncAnimation(fig, update, frames=times, fargs=(sim, ax, trails),
                            interval=interval, blit=False)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f'orbit_simulation_{timestamp}.mp4'
        ani.save(output_file, writer='ffmpeg', dpi=CONFIG['dpi'])
        print(f"\nAnimation saved as {output_file}")
    except Exception as e:
        print(f"Simulation failed: {e}")
    finally:
        plt.close(fig)


if __name__ == "__main__":
    main()
