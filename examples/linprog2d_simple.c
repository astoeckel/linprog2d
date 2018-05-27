/*
 *  linprog2d --- Two-dimensional linear programming solver
 *  Copyright (C) 2018 Andreas Stöckel
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <linprog2d.h>
#include <stdio.h>

/* This program maximizes 5 * x + 10 * y such that
          x          >=   0,
                   y >=   0,
          x          <   15,
      8 * x +  8 * y <  160,
      4 * x + 12 * y <  180                         */
int main() {
	/* Input arrays */
	const double Gx[5U] = {1.0, 0.0, -1.0, -8.0, -4.0};
	const double Gy[5U] = {0.0, 1.0, 0.0, -8.0, -12.0};
	const double h[5U] = {0.0, 0.0, -15.0, -160.0, -180.0};

	/* linprog2d_solve_simple allocated memory for the solver, solves the
	   problem, and frees the memory it allocated. linprog2d provides functions
	   that allow to re-use the same memory for multiple problems, as well as
	   functions that allow to perform manual memory management. */
	const double cx = -5.0;
	const double cy = -10.0;
	linprog2d_result_t res = linprog2d_solve_simple(cx, cy, Gx, Gy, h, 5U);

	/* Print the solution */
	if (res.status == LP2D_POINT) {
		printf("x=%0.2f y=%0.2f\n", res.x1, res.y1);
		return 0;
	}
	printf("Problem is infeasible, unbounded, or not a single point.");
	return 1;
}

