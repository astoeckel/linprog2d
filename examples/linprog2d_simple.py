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

# This program maximizes 5 * x + 10 * y such that
#          x          >=   0,
#                   y >=   0,
#          x          <   15,
#      8 * x +  8 * y <  160,
#      4 * x + 12 * y <  180

import linprog2d

Gx = [1.0, 0.0,  -1.0,   -8.0,   -4.0];
Gy = [0.0, 1.0,   0.0,   -8.0,  -12.0];
h  = [0.0, 0.0, -15.0, -160.0, -180.0];
res = linprog2d.solve(-5.0, -10.0, Gx, Gy, h)
if res.status == linprog2d.STATUS_POINT:
    print("x =", res.x1, "y =", res.y1)
