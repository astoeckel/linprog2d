#!/usr/bin/env python3

#   linprog2d --- Two-dimensional linear programming solver
#   Copyright (C) 2018 Andreas Stöckel
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

import cvxopt
cvxopt.solvers.options['show_progress'] = False

import linprog2d
import numpy as np
import time

import matplotlib.pyplot as plt

def generate_2d_linprog_problem(n, seed):
    rnd = np.random.RandomState(seed)

    c = rnd.uniform(-1, 1, 2) # Direction
    x = rnd.uniform(-1, 1, 2) # Random point

    # Generate a random constraint, if the constraint is not compatible with the
    # random point, flip the constraint
    G = np.empty((2, n))
    h = np.empty(n)
    for i in range(n):
        Gi = rnd.uniform(-100, 100, 2) # Direction
        hi = rnd.uniform(-100, 100) # Offset
        if Gi @ x < hi:
            Gi = -Gi
            hi = -hi
        G[:, i] = Gi
        h[i] = hi

    return c, G, h

def benchmark(ns, repeat=10):
    ts = np.empty((2, len(ns), repeat))
    diff = np.empty((len(ns), repeat))
    for i, n in enumerate(ns):
        for j in range(repeat):
            # Generate a problem
            seed = (i + 6328) * 527 * (j + 383)
            c, G, h = generate_2d_linprog_problem(n, seed)

            # linprog2d
            t10 = time.perf_counter()
            res_lp2 = linprog2d.solve(c[0], c[1], G[0], G[1], h)
            t11 = time.perf_counter()
            ts[0, i, j] = t11 - t10

            # cvxopt
            ccvx = cvxopt.matrix(c)
            Gcvx = cvxopt.matrix(-G.T)
            hcvx = cvxopt.matrix(-h)
            t20 = time.perf_counter()
            res_cvx = cvxopt.solvers.lp(ccvx, Gcvx, hcvx)['x']
            t21 = time.perf_counter()

            ts[1, i, j] = t21 - t20
            diff[i, j] = (res_lp2.x1 - res_cvx[0]) ** 2 + (res_lp2.y1 - res_cvx[1]) ** 2
    return ts, diff

# Gatheer data

ns = np.array(np.logspace(2, 6, 10), dtype=np.int32)
ts, diff = benchmark(ns)

# Plot timing data

fig, ax = plt.subplots()
ys0 = np.mean(ts[0], axis=1)
m0, b0 = np.polyfit(ns, ys0, deg=1)

ys1 = np.mean(ts[1], axis=1)
m1, b1 = np.polyfit(ns, ys1, deg=1)

b0, b1 = np.clip(b0, 0, None), np.clip(b1, 0, None)

ax.loglog(ns, ns * m0 + b0, '--', color='k')
ax.loglog(ns, ys0, 'd-', label='linprog2d')

ax.loglog(ns, ns * m1 + b1, '--', color='k')
ax.loglog(ns, ys1, 'o-', label='cvxopt')

ax.set_ylabel('Processing time in seconds (◀ better)')
ax.set_xlabel('Problem size $n$ (number of constraints)')

ax.legend()

fig.savefig('linprog2d_performance.png', dpi=200)


# Print error data

fig, ax = plt.subplots()
ys0 = np.sqrt(np.max(diff, axis=1))
ys1 = np.sqrt(np.median(diff, axis=1))

ax.loglog(ns, ys0, '+-', color='k', label='Max')
ax.loglog(ns, ys1, 'o-', color='k', label='Median')

ax.set_ylabel('Max difference (◀ better)')
ax.set_xlabel('Problem size $n$ (number of constraints)')

ax.legend()

fig.savefig('linprog2d_difference.png', dpi=200)
