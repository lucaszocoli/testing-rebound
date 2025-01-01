import numpy as np
import matplotlib.pyplot as plt
import rebound
import time

sim = rebound.Simulation()
sim.add(m=1.0)

planets = [
    {"name": "Mercúrio", "m": 0.166e-6, "a": 0.39, "e": 0.205, "inc": 7.0, "omega": 48.3},
    {"name": "Vênus", "m": 2.45e-6, "a": 0.72, "e": 0.007, "inc": 3.39, "omega": 76.7},
    {"name": "Terra", "m": 3.0e-6, "a": 1.0, "e": 0.017, "inc": 0.0, "omega": 0.0},
    {"name": "Marte", "m": 0.32e-6, "a": 1.52, "e": 0.093, "inc": 1.85, "omega": 49.5},
    {"name": "Júpiter", "m": 954.79e-6, "a": 5.2, "e": 0.049, "inc": 1.3, "omega": 273.0},
    {"name": "Saturno", "m": 285.88e-6, "a": 9.58, "e": 0.056, "inc": 2.49, "omega": 336.0},
    {"name": "Urano", "m": 43.66e-6, "a": 19.22, "e": 0.046, "inc": 0.77, "omega": 98.0},
    {"name": "Netuno", "m": 51.51e-6, "a": 30.05, "e": 0.010, "inc": 1.77, "omega": 131.8},
]

for planet in planets:
    sim.add(
        m=planet["m"],
        a=planet["a"],
        e=planet["e"],
        inc=np.deg2rad(planet["inc"]),
        omega=np.deg2rad(planet["omega"])
    )

sim.move_to_com()


p = sim.particles
times = np.linspace(0, 1, 1000)
smas = np.full((len(p)-1, len(times)), np.nan)
ds = np.full((len(p)-1, len(times)), np.nan)
time0 = time.time()

for i, t in enumerate(times):
    sun_loc = [p[0].x, p[0].y, p[0].z]
    plt.plot(p[0].x, p[0].y, '.', color='C0', alpha= i/len(times))

    for j in range(1, len(p)):
        plt.plot(p[j].x, p[j].y, '.', color='C'+str(j), alpha= i/len(times))
        ds[j-1, i] = np.sqrt(
            (p[j].x - sun_loc[0])**2 +
            (p[j].y - sun_loc[1])**2 +
            (p[j].z - sun_loc[2])**2
        )
        smas[j-1, i] = p[j].a

    sim.integrate(t)

lim = 35
plt.xlim(-lim, lim)
plt.ylim(-lim, lim)
plt.show()

print(time.time() - time0)

plt.figure()
plt.plot(times, smas.T)

plt.figure()
plt.plot(times, ds.T)

plt.show()