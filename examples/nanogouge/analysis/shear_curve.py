#!/usr/bin/env python3
"""Friction coefficient from shear_<sigma>.dat (columns: step, gamma, shear stress MPa, normal stress MPa, mu, Lz, rho_b).
    python shear_curve.py shear_5.6.dat shear_56.dat            -> shear_curve.png and the mean mu of the last third of the strain
"""
import sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for fn in sys.argv[1:]:
    d = np.loadtxt(fn, comments='#')
    if d.ndim < 2 or len(d) < 3: print(fn, 'too few rows'); continue
    g, ss, sn, mu, lz, rho = d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5], d[:, 6]
    m = g > g.max() * 2 / 3; print(f'{fn}: gamma up to {g.max():.2f}; last third: sigma_n = {sn[m].mean():.1f} MPa, tau = {ss[m].mean():.2f} MPa, mu = {ss[m].mean() / sn[m].mean():.4f}; peak mu {mu.max():.4f} at gamma {g[mu.argmax()]:.2f}')
    ax[0].plot(g, ss, label=fn); ax[1].plot(g, mu); ax[2].plot(g, rho)
ax[0].set_xlabel('shear strain'); ax[0].set_ylabel('shear stress (MPa)'); ax[0].legend(fontsize=7)
ax[1].set_xlabel('shear strain'); ax[1].set_ylabel('friction coefficient mu')
ax[2].set_xlabel('shear strain'); ax[2].set_ylabel('dry bulk density (g/cm3)')
plt.tight_layout(); plt.savefig('shear_curve.png', dpi=110)
