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

/**
 * @file linprog2d.h
 *
 * C and C++ header for the 2D linear programming solver.
 *
 * @author Andreas Stöckel
 */

#ifndef LINPROG_2D_H_
#define LINPROG_2D_H_

#ifdef __cplusplus
extern "C" {
#endif
/**
 * Enum describing the type of the solution to a linear programming problem.
 */
enum linprog2d_status {
	/**
	 * There is not enough memory in the given linprog2d instance to solve the
	 * problem or the gradient is zero.
	 */
	LP2D_ERROR,

	/**
	 * There is no solution to this problem, i.e. the solution space is empty.
	 * This may happen if there are constraints that contradict each other.
	 */
	LP2D_INFEASIBLE,

	/**
	 * The problem is unbounded. The result should be interpreted as: any point
	 * starting at (x1, y1) in direction of the vector (x2, y2) is a valid
	 * solution.
	 */
	LP2D_UNBOUNDED,

	/**
	 * There is an entire edge along which the solution is optimal. The edge is
	 * described by the two points (x1, y1) - (x2, y2).
	 */
	LP2D_EDGE,

	/**
	 * The solution is a single point stored in (x1, y1).
	 */
	LP2D_POINT
};

/**
 * Result of the linear programming algorithm.
 */
struct linprog2d_result {
	/**
	 * Enum describing how the fields of this structure should be interpreted.
	 */
	enum linprog2d_status status;

	/**
	 * The result is encoded as two points. If the optimum is a single point,
	 * this point is stored as (x1, y1). Otherwise, if the optimum is a line or
	 * a ray, the result is encoded in the pair (x1, y1), (x2, y2). Also see the
	 * documentation for linprog2d_result_type for more details.
	 */
	double x1, y1, x2, y2;
};

/**
 * Typedef of the above structure.
 */
typedef struct linprog2d_result linprog2d_result_t;

/**
 * Opaque type used to represent a linprog2d instance.
 */
typedef void linprog2d_t;

/**
 * Computes the number of bytes required to store a Linprog2DSolver instance
 * with the given capacity.
 */
unsigned int linprog2d_mem_size(unsigned int capacity);

/**
 * Constructs a linprog2d instance with the given capacity inplace at the
 * given memory location. The required size of the memory region can be computed
 * by calling linprog2d_mem_size(). You may want to use linprog2d_create()
 * instead, which automatically allocates the required amount of memory on the
 * stack.
 *
 * @param capacity is the number of constraints that the newly created linprog2d
 * instance should be able to handle.
 * @param mem is a pointer at the memory region the linprog2d instance should be
 * written to. The caller must make sure that enough space is available in this
 * memory region. The required space can be computed by calling
 * linprog2d_create.
 */
linprog2d_t *linprog2d_init(unsigned int capacity, char *mem);

/**
 * Creates a new linprog2d instance that is able to represent at least n
 * constraints. The returned pointer must be freed using linprog2d_free. If a
 * failure occurs (out of memory) or the library has been compiled without
 * linking to the C standard library (the LINPROG2D_NO_ALLOC flag is defined),
 * returns null.
 */
linprog2d_t *linprog2d_create(unsigned int capacity);

/**
 * Frees a previously created linprog2d instance. If the LINPROG2D_NO_ALLOC flag
 * is defined, this function does nothing.
 */
void linprog2d_free(linprog2d_t *prog);

/**
 * Returns the maximum number of constraints in a problem that can be solved
 * with this linprog2d_t instance.
 */
unsigned int linprog2d_capacity(const linprog2d_t *prog);

/**
 * Solves a two-dimensional linear programming problem.
 */
linprog2d_result_t linprog2d_solve(linprog2d_t *prog, double cx, double cy,
                                   const double *Gx, const double *Gy,
                                   const double *h, unsigned int n);

/**
 * Convenience function which allocates a new linprog2d_t instance, calls
 * its solve function, destroys the instance and returns the result. If you
 * want to reuse the same linprog2d_t instance, use linprog2d_create,
 * linprog2d_solve, and linprog2d_free. Solves a two-dimensional linear
 * programming problem of the form.
 *
 * minimize cx * x + cy * y
 *   w.r.t. Gx[i] * x + Gy[i] * y >= h[i] for all 0 <= i < n
 *
 * The result is encoded in the returned linprog2d_result structure.
 */
linprog2d_result_t linprog2d_solve_simple(double cx, double cy,
                                          const double *Gx, const double *Gy,
                                          const double *h, unsigned int n);

#ifdef __cplusplus
}
#endif

#endif /* LINPROG_2D_H_ */
