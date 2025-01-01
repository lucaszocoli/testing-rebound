import numpy as np
import matplotlib.pyplot as plt
import rebound
from astropy import units as u
import time

sim = rebound.Simulation()

sim.units = ('yr', 'AU', 'Msun')

## Add stars
sim.add(m=0.867)
sim.add(m = (1*u.Mearth).to_value(u.Msun),
        a = .7,
        e = .05,
        omega = 6.01)

sim.add(m = (5*u.Mearth).to_value(u.Msun),
        a = 2.1,
        e = .1,
        omega = 2.85,
        inc = (6.5*u.deg).to_value(u.rad),
        Omega = 4.05)

sim.add(m = (93*u.Mearth).to_value(u.Msun),
        a = 6.3,
        e = .3,
        omega = 5.86,
        inc = (21*u.deg).to_value(u.rad),
        Omega = 5.72)

sim.add(m = (.6*u.Mearth).to_value(u.Msun),
        a = .51,
        e = .12,
        omega = 5.37,
        inc = (-12*u.deg).to_value(u.rad),
        Omega = 2.14)

sim.add(primary = sim.particles[1],
        m = (2.1e17*u.kg).to_value(u.Msun),
        a = (27000*u.km).to_value(u.AU),
        e = .09,
        omega = 2.35,
        inc = (5.2*u.deg).to_value(u.rad),
        Omega = 4.87)

sim.add(primary = sim.particles[1],
        m = (4.7e15*u.kg).to_value(u.Msun),
        a = (150000*u.km).to_value(u.AU),
        e = .23,
        omega = 3.78,
        inc = (-17*u.deg).to_value(u.rad),
        Omega = 3.98)

sim.add(primary = sim.particles[1],
        m = (6.1e20*u.kg).to_value(u.Msun),
        a = (232000*u.km).to_value(u.AU),
        e = .04,
        omega = 5.10,
        inc = (1.4*u.deg).to_value(u.rad),
        Omega = 0.77)

fig = rebound.OrbitPlot(sim)
op1 = rebound.OrbitPlot(sim, particles = [1,2,3,4])
op2 = rebound.OrbitPlot(sim, particles = [5,6,7], primary = 1, show_primary = False)

#op = rebound.OrbitPlot(sim, color=True, periastron=True)
plt.show()

p = sim.particles

times = np.linspace(0, 1000)
smas = np.full((len(p)-1,len(times)), np.nan)
ds = np.full((len(p)-1,len(times)), np.nan)

time0 = time.time()
for i, t in enumerate(times):
    print(t)
    star_loc = [p[0].x, p[0].y, p[0].z]
    for j in range(1, len(p)):
        ds[j-1, i] = np.sqrt(
            (p[j].x-star_loc[0])**2 +
            (p[j].y-star_loc[1])**2 +
            (p[j].z-star_loc[2])**2 
        )
        if j<=4:
            smas[j-1,i] = p[j].a
        else:
            orb = p[j].orbit(primary=p[1])
            smas[j-1,i] = orb.a

    sim.integrate(t)

print(time.time() - time0)

plt.figure()
plt.plot(times, smas.T)

plt.figure()
plt.plot(times, ds.T)

plt.show()