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

#include "linprog2d.h"

#include <math.h>

#ifndef LINPROG2D_NO_ALLOC
#include <stdlib.h>
#endif

/******************************************************************************
 * PRIVATE HELPER FUNCTIONS                                                   *
 ******************************************************************************/

/******************************************************************************
 * Constants/math functions missing in ANSI C                                 *
 ******************************************************************************/

typedef int bool_t;

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef NULL
#define NULL ((void *)0)
#endif

static double fmax_(double x, double y) { return (x > y) ? x : y; }
static double fmin_(double x, double y) { return (x < y) ? x : y; }
static double hypot_(double x, double y) { return sqrt(x * x + y * y); }

/******************************************************************************
 * Floating point comparison                                                  *
 ******************************************************************************/

/* There is no obvious reason why I chose these values. They seem fairly
   reasonable, but depending on the problem domain they may just be wrong. */
#define MAX_EPS_ABS 1e-30 /* maximum absolute difference */
#define MAX_EPS_REL 1e-15 /* maximum relative difference */

static bool_t feq_(double x, double y) {
	const double dlt = fabs(x - y);
	return (dlt < MAX_EPS_ABS) || (dlt < MAX_EPS_REL * fmax_(fabs(x), fabs(y)));
}

/******************************************************************************
 * 2D vector and matrix code                                                  *
 ******************************************************************************/

/**
 * Vector in 2D space.
 */
struct vec2 {
	double x, y;
};

/**
 * Returns a vec2 instance that is filled with the given values.
 */
static struct vec2 vec2_create(double x, double y) {
	struct vec2 res;
	res.x = x;
	res.y = y;
	return res;
}

/**
 * Structure used internally to represent a 2x2 matrix of doubles. This is used
 * to rotate the entire problem space to a canonical form.
 */
struct mat22 {
	/**
	 * Matrix entries, a[i, j], where i is the row, j the column.
	 */
	double a11, a12, a21, a22;
};

/**
 * Fills a mat22 instance with the given values.
 */
static struct mat22 mat22_create(double a11, double a12, double a21,
                                 double a22) {
	struct mat22 res;
	res.a11 = a11;
	res.a12 = a12;
	res.a21 = a21;
	res.a22 = a22;
	return res;
}

/**
 * Returns a rotation matrix that aligns the given (x, y) vector with the
 * y axis; i.e. (x, y) is rotated to some (0, -y').
 */
static struct mat22 mat22_rot(double x, double y) {
	const double h = hypot_(x, y);
	return mat22_create(y / h, -x / h, x / h, y / h);
}

/******************************************************************************
 * Memory helper functions                                                    *
 ******************************************************************************/

/**
 * Aligns the given pointer plus the offset in bytes to the next cache-line
 * boundary. On most architectures 64 is the length of a cacheline. As a
 * side-effect this means that all arrays are aligned to a 16-byte boundary,
 * which allows parallelisation in some circumstances.
 *
 * TODO: use __builtin_assume_aligned on GCC in the appropriate places.
 */
static void *mem_align64(void *p, unsigned long int offs) {
	return (void *)(((unsigned long int)p + offs + 63UL) & (~63UL));
}

/******************************************************************************
 * Result datastructure helper functions                                      *
 ******************************************************************************/

/**
 * Fills the result structure with the given values.
 */
static linprog2d_result_t linprog2d_result_create(enum linprog2d_status status,
                                                  double x1, double y1,
                                                  double x2, double y2) {
	linprog2d_result_t res;
	res.status = status, res.x1 = x1, res.y1 = y1, res.x2 = x2, res.y2 = y2;
	return res;
}

/**
 * Creates a linprog2d_result instance representing an error.
 */
static linprog2d_result_t linprog2d_result_err() {
	return linprog2d_result_create(LP2D_ERROR, 0.0, 0.0, 0.0, 0.0);
}

/**
 * Creates a linprog2d_result instance representing an infeasible problem.
 */
static linprog2d_result_t linprog2d_result_infeasible() {
	return linprog2d_result_create(LP2D_INFEASIBLE, 0.0, 0.0, 0.0, 0.0);
}

static linprog2d_result_t linprog2d_result_unbounded() {
	return linprog2d_result_create(LP2D_UNBOUNDED, 0.0, 0.0, 0.0, 0.0);
}

/**
 * Used internally by linprog2d_result_point() and linprog2d_result_edge() to
 * transform the result back to the original coordinate system.
 */
static void linprog2d_result_transform_back(const struct mat22 *R,
                                            const struct vec2 *o, double *x,
                                            double *y) {
	double xt, yt;
	xt = *x + o->x, yt = *y + o->y; /* Undo the offset */

	/* Rotate back by multiplying with the inverse of R, which happens to be the
	   transpose (since R is a rotation matrix) */
	*x = R->a11 * xt + R->a21 * yt;
	*y = R->a12 * xt + R->a22 * yt;
}

/**
 * Transforms the given point (x, y) back to the original coordinate system and
 * returns a corresponding linprog2d_result instance.
 */
static linprog2d_result_t linprog2d_result_point(const struct mat22 *R,
                                                 const struct vec2 *o, double x,
                                                 double y) {
	linprog2d_result_transform_back(R, o, &x, &y);
	return linprog2d_result_create(LP2D_POINT, x, y, 0.0, 0.0);
}

/**
 * Transforms the given point (x, y) back to the original coordinate system and
 * returns a corresponding linprog2d_result instance.
 */
static linprog2d_result_t linprog2d_result_edge(const struct mat22 *R,
                                                const struct vec2 *o, double x1,
                                                double y1, double x2,
                                                double y2) {
	linprog2d_result_transform_back(R, o, &x1, &y1);
	linprog2d_result_transform_back(R, o, &x2, &y2);
	return linprog2d_result_create(LP2D_EDGE, x1, y1, x2, y2);
}

/******************************************************************************
 * Linear time median algorithm                                               *
 ******************************************************************************/

/* See http://pages.ripco.net/~jgamble/nw.html as well as
   https://stackoverflow.com/a/2789530 */

/* Function prototype */
static double median(double *d, unsigned int len);

#define SWAP(x, y)         \
	{                      \
		double tmp = d[x]; \
		d[x] = d[y];       \
		d[y] = tmp;        \
	}

/**
 * Sorts a list of up to five elements in "constant time" (well, virtually ANY
 * algorithm sorts a size-constrained list in "constant time").
 */
static void sort(double *d, unsigned int len) {
#define SWAP_IF_GT(x, y) \
	if (d[y] < d[x])     \
	SWAP(x, y)
	switch (len) {
		case 0U:
		case 1U:
			break;
		case 2U:
			SWAP_IF_GT(0U, 1U);
			break;
		case 3U:
			SWAP_IF_GT(1U, 2U);
			SWAP_IF_GT(0U, 2U);
			SWAP_IF_GT(0U, 1U);
			break;
		case 4U:
			SWAP_IF_GT(0U, 1U);
			SWAP_IF_GT(2U, 3U);
			SWAP_IF_GT(0U, 2U);
			SWAP_IF_GT(1U, 3U);
			SWAP_IF_GT(1U, 2U);
			break;
		case 5U:
			SWAP_IF_GT(0U, 1U);
			SWAP_IF_GT(3U, 4U);
			SWAP_IF_GT(2U, 4U);
			SWAP_IF_GT(2U, 3U);
			SWAP_IF_GT(0U, 3U);
			SWAP_IF_GT(0U, 2U);
			SWAP_IF_GT(1U, 4U);
			SWAP_IF_GT(1U, 3U);
			SWAP_IF_GT(1U, 2U);
			break;
	}
#undef SWAP_IF_GT
}

/**
 * Partitions the given list in such a way, that all values smaller than the
 * given piviot are at the beginning of the list, all values larger than the
 * piviot are at the end of the list, and the piviot itself is between the two
 * lists. Returns the number of values smaller than the piviot.
 */
static unsigned int partition(double *d, unsigned int len, double piviot) {
	unsigned int i, l = 0, r = len - 1;
	for (i = 0; i <= r;) {
		if (d[i] < piviot) {
			SWAP(l, i);
			l++;
			i++;
		} else if (d[i] > piviot) {
			SWAP(r, i);
			r--;
		} else {
			i++;
		}
	}
	return l;
}

/**
 * Computes the kth-smallest element in the list d with length len. Operates
 * inline on d.
 */
static double kth_smallest(double *d, unsigned int len, unsigned int k) {
	/* See http://www-di.inf.puc-rio.br/~laber/median-lineartime.pdf */
	unsigned int i, j, l;
	double piviot; /* median-of-medians */

	/* If the list has less than five entries, just sort the list and pick the
	   entry referenced by k */
	if (len <= 5) {
		sort(d, len);
		return d[k];
	}

	/* Iterate over the list in groups of five; we can ignore the elements at
	   the end if the list size is not perfectly divisible by five. */
	j = 0;
	for (i = 0; i + 5 <= len; i += 5, j++) {
		/* Compute the median of the 5-element sub-list in constant time. */
		median(d + i, 5);

		/* Move the median of this sub-list (which is at i + 2) to the beginning
		   of the list. */
		SWAP(i + 2, j);
	}

	/* Compute the median of the medians at the beginning of the list. */
	piviot = median(d, j);

	/* Compute the median by pivioting around the median-of-medians; partition
	   the array in such a way that all entries smaller than the piviot are in
	   the left half of the list, and all entries larger than the piviot are in
	   the right half of the list. */
	l = partition(d, len, piviot);
	if (l == k) {
		return piviot; /* the piviot happens to be the median */
	} else if (l > k) {
		return kth_smallest(d, l, k);
	} else {
		return kth_smallest(d + l + 1, len - l - 1, k - l - 1);
	}
}

/**
 * Returns the element which, if the list d were sorted, was at position len / 2
 */
static double median(double *d, unsigned int len) {
	return kth_smallest(d, len, len / 2);
}

#undef SWAP

/******************************************************************************
 * Actual implementation of the 2D linprog algorithm                          *
 ******************************************************************************/

/**
 * Internally used structure holding all the data associated with a linprog2d
 * instance.
 */
struct linprog2d_data {
	/**
	 * Pointer at the x-part of the LHS of the constraints.
	 */
	double *Gx;

	/**
	 * Pointer at the y-part of the LHS of the constraints.
	 */
	double *Gy;

	/**
	 * Pointer at the RHS of the constraints.
	 */
	double *h;

	/**
	 * Slopes of the individual constraints.
	 */
	double *dx;

	/**
	 * y-axis offset of the individual constraints.
	 */
	double *y0;

	/**
	 * x-coordinates of the constraint intersections. This list has only
	 * capacity / 2 entries. There can only be capacity / 2 intersections.
	 */
	double *x_intersect;

	/**
	 * Array of indices corresponding to the ceiling constraints.
	 */
	unsigned int *ceil;

	/**
	 * Array of indices corresponding to the floor constraints.
	 */
	unsigned int *floor;

	/**
	 * Temporarily used memory for storing the new ceil/floor constraints in
	 * linprog2d_calculate_intersects().
	 */
	unsigned int *tmp;

	/**
	 * Current left and right boundaries. Solutions must be to the left/right of
	 * these boundaries.
	 */
	double x0, x1;

	/**
	 * Number of valid constraints in the individual lists.
	 */
	unsigned int ceil_len, floor_len, intersect_len;

	/**
	 * Number of elements that can be stored in the arrays Gx, Gy, and h.
	 */
	unsigned int capacity;

	/**
	 * Rotation matrix in the current problem.
	 */
	struct mat22 R;

	/**
	 * Offset vector for the current problem.
	 */
	struct vec2 o;

	/**
	 * Number of constraints in the current problem.
	 */
	unsigned int n;
};

typedef struct linprog2d_data linprog2d_data_t;

/**
 * Function that clears problem-specific data from the linprog2d_data structure.
 * This must be called at the beginning of the solution process.
 */
static void linprog2d_reset(linprog2d_data_t *prog, unsigned int n) {
	prog->ceil_len = 0;
	prog->floor_len = 0;
	prog->x0 = -HUGE_VAL;
	prog->x1 = HUGE_VAL;
	prog->ceil_len = 0;
	prog->floor_len = 0;
	prog->intersect_len = 0;
	prog->R = mat22_create(0.0, 0.0, 0.0, 0.0);
	prog->o = vec2_create(0.0, 0.0);
	prog->n = n;
}

static linprog2d_t *linprog2d_init_internal(linprog2d_data_t *prog,
                                            unsigned int capacity, char *mem) {
#define SD sizeof(double)
#define SU sizeof(unsigned int)
	if (!prog) {
		return NULL;
	}

	/* Calculate the offsets for the individual arrays from the continuous
	   piece of memory passed to this function */
	prog->Gx = (double *)mem_align64(mem, 0U);
	prog->Gy = (double *)mem_align64(prog->Gx, SD * capacity);
	prog->h = (double *)mem_align64(prog->Gy, SD * capacity);
	prog->dx = (double *)mem_align64(prog->h, SD * capacity);
	prog->y0 = (double *)mem_align64(prog->dx, SD * capacity);
	prog->x_intersect = (double *)mem_align64(prog->y0, SD * capacity);
	prog->ceil =
	    (unsigned int *)mem_align64(prog->x_intersect, SD * capacity / 2);
	prog->floor = (unsigned int *)mem_align64(prog->ceil, SU * capacity);
	prog->tmp = (unsigned int *)mem_align64(prog->floor, SU * capacity);
	prog->capacity = capacity;

	/* Reset all other fields to their initial values */
	linprog2d_reset(prog, 0U);

	return prog;
#undef SD
#undef SU
}

/**
 * Computes the value that should be used to normalize the given constraint.
 * Currently this normalisation is only based on the direction (Gx, Gy) of the
 * constraint, not on the offset h.
 */
static double linprog2d_normalization_coeff(double Gx, double Gy) {
	return fmax_(fabs(Gx), fabs(Gy));
}

/**
 * Takes the linear program provided by the user and rotates it such that the
 * gradient is aligned with the y-axis. Additionally, normalizes all constraints
 * such that the maximum coefficient is one in each line. Furthermore shifts the
 * entire problem space such that all constraints are centered around the orign.
 * This centering is performed by finding an offset vector (o.x, o.y) s.t.
 * the following expression is minimized
 *
 * sum_i (h[i] - Gx[i] * o.x - Gy[i] * o.y)^2
 *
 * The closed-form solution to this least-squares optimization problem is
 *
 * o = (G.T * G)^-1 * G.T * h,
 *
 * which can be computed in linear time with constant memory, since G.T * G is a
 * 2x2 matrix.
 */
static int linprog2d_condition_problem(linprog2d_data_t *prog, double cx,
                                       double cy, const double *src_Gx,
                                       const double *src_Gy,
                                       const double *src_h) {
	struct mat22 R = mat22_rot(cx, cy);
	struct vec2 o = vec2_create(0.0, 0.0);               /* Offset vector */
	struct mat22 GTG = mat22_create(0.0, 0.0, 0.0, 0.0); /* Matrix G.T G */
	struct vec2 GTc = vec2_create(0.0, 0.0);             /* Vector G.T c */
	double Gx, Gy, h, norm, GTG_det;                     /* Temp variables */
	unsigned int i_tar = 0, i = 0;
	double *tar_Gx = prog->Gx, *tar_Gy = prog->Gy, *tar_h = prog->h;

	/* Copy the memory from the source to the target location; rotate all the
	   source vectors. At the same time normalize the problem such that the
	   coefficient with the largest absolute value is scaled to +-1. */
	for (i = 0; i < prog->n; i++) {
		/* Rotate the constraint direction on the left-hand side */
		Gx = R.a11 * src_Gx[i] + R.a12 * src_Gy[i];
		Gy = R.a21 * src_Gx[i] + R.a22 * src_Gy[i];
		h = src_h[i];

		/* Skip invalid constraints */
		if (feq_(Gx, 0.0) && feq_(Gy, 0.0)) {
			if (h <= 0.0) {
				/* Constraint of the form 0 >= h is always true for h <= 0.0 */
				continue;
			} else {
				/* This constraint is always false. Abort. */
				return FALSE;
			}
		}

		/* Normalize the constraints by dividing both the right- and left-hand
		   side by the largest direction coefficient. */
		norm = linprog2d_normalization_coeff(Gx, Gy);
		Gx /= norm, Gy /= norm, h /= norm;

		/* Update the matrix G.T * G */
		GTG.a11 += Gx * Gx;
		GTG.a12 += Gx * Gy; /* Same as a21 */
		GTG.a22 += Gy * Gy;

		/* Update the matrix G.T * h */
		GTc.x += Gx * h;
		GTc.y += Gy * h;

		/* Write the result to memory and increment the write pointer */
		tar_Gx[i_tar] = Gx, tar_Gy[i_tar] = Gy, tar_h[i_tar] = h;
		i_tar++;
	}

	/* Invert the GTG matrix (if possible) and compute o. The GTG is not
	   invertible if there is an infinite number of possible offsets that
	   minimize the error function. This is for example the case if there is
	   only one constraint. We just don't do the offsetting in this case, which
	   is only meant to help with numerical stability. */
	GTG_det = GTG.a11 * GTG.a22 - GTG.a12 * GTG.a12;
	if (GTG_det != 0.0) {
		o.x = (GTG.a22 * GTc.x - GTG.a12 * GTc.y) / GTG_det;
		o.y = (-GTG.a12 * GTc.x + GTG.a22 * GTc.y) / GTG_det;
	}

	/* Update the linear program data */
	prog->n = i_tar; /* Constraints may have been eliminated */
	prog->R = R;
	prog->o = o;

	/* Update h according to the computed offset vector */
	for (i = 0; i < prog->n; i++) {
		tar_h[i] -= o.x * tar_Gx[i] + o.y * tar_Gy[i];
	}

	return TRUE; /* Success */
}

#define CAT_VERT_LEFT 0
#define CAT_VERT_RIGHT 1
#define CAT_CEIL 2
#define CAT_FLOOR 3

/**
 * Calculates the category of a constraint based on its direction Gx, Gy.
 */
static int linprog2d_constraint_category(double Gx, double Gy) {
	if (feq_(Gy, 0.0)) {
		if (Gx > 0.0) {
			return CAT_VERT_LEFT;
		} else {
			return CAT_VERT_RIGHT;
		}
	} else if (Gy > 0.0) {
		return CAT_FLOOR;
	} else /* if (Gy < 0.0) */ {
		return CAT_CEIL;
	}
}

/**
 * Sorts the constraints into the ceil and floor lists and updates the left
 * and right boundary if a constraint is perfectly vertical.
 */
static int linprog2d_categorize_constraints(linprog2d_data_t *prog) {
	unsigned int i;
	const double *Gx = prog->Gx, *Gy = prog->Gy, *h = prog->h;
	for (i = 0; i < prog->n; i++) {
		switch (linprog2d_constraint_category(Gx[i], Gy[i])) {
			case CAT_VERT_LEFT:
				prog->x0 = fmax_(prog->x0, h[i] / Gx[i]);
				break;
			case CAT_VERT_RIGHT:
				prog->x1 = fmin_(prog->x1, h[i] / Gx[i]);
				break;
			case CAT_CEIL:
				prog->ceil[prog->ceil_len++] = i;
				break;
			case CAT_FLOOR:
				prog->floor[prog->floor_len++] = i;
				break;
		}
	}
	return prog->x0 <= prog->x1;
}

/**
 * For each non-vertical constraint in the given list computes the slope.
 */
static void linprog2d_calculate_yoffset_form(const unsigned int *idcs,
                                             unsigned int idcs_len,
                                             const double *Gx, const double *Gy,
                                             const double *h, double *dx,
                                             double *y0) {
	unsigned int i;
	for (i = 0; i < idcs_len; i++) {
		dx[idcs[i]] = -Gx[idcs[i]] / Gy[idcs[i]];
		y0[idcs[i]] = h[idcs[i]] / Gy[idcs[i]];
	}
}

/**
 * Calculates the intersection point between two constraints.
 */
static int linprog2d_calculate_intersect(double Gx1, double Gy1, double h1,
                                         double Gx2, double Gy2, double h2,
                                         double *x, double *y) {
	const double num_x = h1 * Gy2 - h2 * Gy1;
	const double num_y = h2 * Gx1 - h1 * Gx2;
	const double den = Gx1 * Gy2 - Gx2 * Gy1;

	if (feq_(den, 0.0)) {
		return FALSE; /* Lines are parallel */
	}

	*x = num_x / den, *y = num_y / den;
	return TRUE;
}

/**
 * Of the two constraint indices ci0, ci1 returns the index of the constraint
 * that is not redundant.
 */
static unsigned int linprog2d_eliminate_constraint(
    const double *h, const double *dx, unsigned int ci0, unsigned int ci1,
    bool_t is_ceil, bool_t is_parallel, bool_t optimum_is_left) {
	/* Get the constraint types, vertical constraints have already been filtered
	   out. */
	if (is_parallel) {
		/* The following check only works because the normal vector of the
		   constraint has been normalised. */
		if (h[ci0] >= h[ci1]) {
			return ci0;
		} else {
			return ci1;
		}
	} else {
		/* The greater than/less than relation changes depending on whether the
		   optimum is on the right/left side or this is a ceil/floor
		   constraint. */
		int dir = (optimum_is_left ? 1 : -1) * (is_ceil ? 1 : -1);
		if (dir * dx[ci0] >= dir * dx[ci1]) {
			return ci0;
		} else {
			return ci1;
		}
	}
}

/**
 * Calculates intersections for pairs of constraints. If two constraints happen
 * to be parallel or the intersection point lies outside the current left/right
 * boundaries, one of the constraints in the pair is redundant. This code
 * removes one of those contraints.
 */
static void linprog2d_calculate_intersects(linprog2d_data_t *prog,
                                           unsigned int *idcs,
                                           unsigned int *idcs_len,
                                           bool_t is_ceil, bool_t has_median,
                                           double mx, bool_t optimum_is_left) {
	/* We have two write pointers. One at the beginning of the temporary list
	   and one at the end of the list. The first pointer is used to write the
	   indices of pairs that may be feasible. We must ensure that we do not
	   change these pairs between iterations to ensure that thay are eliminated
	   after x0/x1 have been updated in the main loop. The second pointer holds
	   the indices of single constraints, i.e. stemming from those constraint
	   pairs for which we eliminated a constraint. */
	unsigned int i_tar_pair = 0U, i_tar_single = prog->n - 1U;
	unsigned int i, ci0, ci1;
	double x, y;
	const double *Gx = prog->Gx, *Gy = prog->Gy, *h = prog->h, *dx = prog->dx;
	unsigned int *tmp = prog->tmp;

	/* Iterate over pairs of constraints, for each pair compute the intersect */
	for (i = 0U; i < (*idcs_len) / 2U; i++) {
		ci0 = idcs[2 * i + 0], ci1 = idcs[2 * i + 1];
		if (!linprog2d_calculate_intersect(Gx[ci0], Gy[ci0], h[ci0], Gx[ci1],
		                                   Gy[ci1], h[ci1], &x, &y)) {
			tmp[i_tar_single--] = linprog2d_eliminate_constraint(
			    h, dx, ci0, ci1, is_ceil, TRUE, FALSE);
		} else if (x < prog->x0 ||
		           (has_median && feq_(x, mx) && !optimum_is_left)) {
			tmp[i_tar_single--] = linprog2d_eliminate_constraint(
			    h, dx, ci0, ci1, is_ceil, FALSE, FALSE);
		} else if (x > prog->x1 ||
		           (has_median && feq_(x, mx) && optimum_is_left)) {
			tmp[i_tar_single--] = linprog2d_eliminate_constraint(
			    h, dx, ci0, ci1, is_ceil, FALSE, TRUE);
		} else {
			/* As far as we know, the point may lie in the feasible range.
			   Remember the intersection point and store the indicies of the
			   constraints this intersection point belongs to. */
			prog->x_intersect[prog->intersect_len++] = x;
			tmp[i_tar_pair++] = ci0, tmp[i_tar_pair++] = ci1;
		}
	}

	/* The previous number of constraints was uneven. Correspondingly, we didn't
	   test the last constraint. Make sure to add this constraint to the updated
	   constraint list. */
	if ((*idcs_len) & 1U) {
		tmp[i_tar_single--] = idcs[(*idcs_len) - 1U];
	}

	/* Reassmble idcs from target */
	*idcs_len = 0U;
	for (i = 0U; i < i_tar_pair; i++) {
		idcs[(*idcs_len)++] = tmp[i];
	}
	for (i = prog->n - 1U; i > i_tar_single; i--) { /* Note: i_tar_single > 0 */
		idcs[(*idcs_len)++] = tmp[i];
	}
}

/**
 * Structure containing an extreme value and the minimum/maximum slope.
 */
struct linprog2d_extremum {
	/**
	 * Extreme value, either the minimum or the maximum, depending on the value
	 * of the "compute_min" flag passed to linprog2d_track_extrema().
	 */
	double y;

	/**
	 * Minimum/maximum slope.
	 */
	double min_dx, max_dx;

	/**
	 * True if the data was extracted from at least one constraint, false if
	 * there was no constraint to base the data on.
	 */
	bool_t valid;
};

/**
 * For the given constraints, computes the minimum/maximum at the given
 * x-coordinates and tracks the minimum/maximum slope at that point.
 */
static struct linprog2d_extremum linprog2d_track_extrema(
    double x, const double *dx, const double *y0, const unsigned int *idcs,
    unsigned int idcs_len, bool_t compute_min) {
	unsigned int i, j;
	double y;
	struct linprog2d_extremum e;
	e.y = compute_min ? HUGE_VAL : -HUGE_VAL;
	e.min_dx = HUGE_VAL, e.max_dx = -HUGE_VAL;
	e.valid = idcs_len > 0;

	for (i = 0; i < idcs_len; i++) {
		j = idcs[i];           /* actual constraint we're testing */
		y = y0[j] + dx[j] * x; /* evaluate the constraint at x */
		if (feq_(y, e.y)) {
			/* Another constraint going through the same point with a different
			   slope */
			e.max_dx = fmax_(dx[j], e.max_dx);
			e.min_dx = fmin_(dx[j], e.min_dx);
		} else if ((compute_min && y < e.y) || (!compute_min && y > e.y)) {
			e.y = y; /* this is the new extremum, reset the min/max slope */
			e.min_dx = e.max_dx = dx[j];
		}
	}

	return e;
}

#define LOC_INFEASIBLE 0
#define LOC_LEFT 1
#define LOC_RIGHT 2
#define LOC_HERE 3
#define LOC_HERE_EDGE 4

/**
 * Determines where the optimum is w.r.t. the given median mx. This function
 * assumes that there is at least one floor constraint.
 */
static int linprog2d_locate_optimum(linprog2d_data_t *prog, double mx,
                                    double *y) {
	/* Compute the value of the ceil/floor constraints at mx and track their
	   slope. Since multiple constraints may go through exactly the same point,
	   we need to track both the minimum and the maximum slope for all
	   constraints that go through the same extreme point. */
	struct linprog2d_extremum e_ceil, e_floor;
	e_ceil = linprog2d_track_extrema(mx, prog->dx, prog->y0, prog->ceil,
	                                 prog->ceil_len, TRUE);
	e_floor = linprog2d_track_extrema(mx, prog->dx, prog->y0, prog->floor,
	                                  prog->floor_len, FALSE);

	if (e_ceil.valid && e_ceil.y < e_floor.y) {
		/* mx is outside the feasible region, (implicitly) evaluate
		   d/dx f(x) - g(x) */
		if (e_floor.min_dx > e_ceil.max_dx) {
			return LOC_LEFT;
		} else if (e_floor.max_dx < e_ceil.min_dx) {
			return LOC_RIGHT;
		}
		return LOC_INFEASIBLE;
	}

	if (feq_(e_floor.min_dx, 0.0) && !feq_(e_floor.max_dx, 0.0)) {
		/* Solution is an edge, but this is the right-most point. */
		return LOC_LEFT;
	} else if (feq_(e_floor.max_dx, 0.0) && !feq_(e_floor.min_dx, 0.0)) {
		/* Solution is an edge, but this is the left-most point. */
		return LOC_RIGHT;
	} else if (feq_(e_floor.max_dx, 0.0) && feq_(e_floor.min_dx, 0.0)) {
		/* This one is tough. The floor is horizontal, which means that the
		   solution is an edge, but there is no intersection with another
		   floor constraint that would allow us to progress naturally. We
		   must compute the intersection between the horizontal floor and
		   all other floors/ceils and return the min/max. Signal this by
		   returning LOC_HERE_EDGE. */
		return LOC_HERE_EDGE;
	} else if (e_floor.min_dx < 0.0 && e_floor.max_dx > 0.0) {
		/* Vee-shape. This is the solution */
		*y = e_floor.y;
		return LOC_HERE;
	} else if (e_floor.min_dx > 0.0) {
		return LOC_LEFT;
	} else {
		return LOC_RIGHT;
	}
}

/**
 * Used internally in linprog2d_calculate_edge to check intersections between
 * the top-most horizontal floor constraint and all other ceil/floor
 * constraints.
 */
static void linprog2d_calculate_edge_intersections(linprog2d_data_t *prog,
                                                   const unsigned int *idcs,
                                                   unsigned int idcs_len,
                                                   unsigned int if0, double mx,
                                                   bool_t is_ceil) {
	const double *Gx = prog->Gx, *Gy = prog->Gy, *h = prog->h, *dx = prog->dx;
	double rx1, ry1;
	unsigned int i;

	/* Iterate over all floor and ceiling constraints and calculate the
	   intersection with our ceiling constraint. */
	for (i = 0; i < idcs_len; i++) {
		unsigned int j = idcs[i];
		if (j == if0) { /* Skip self-itersections */
			continue;
		}
		if (linprog2d_calculate_intersect(Gx[if0], Gy[if0], h[if0], Gx[j],
		                                  Gy[j], h[j], &rx1, &ry1)) {
			if (((is_ceil && dx[j] > 0.0) || (!is_ceil && dx[j] < 0.0)) &&
			    rx1 > prog->x0) {
				prog->x0 = rx1;
			}
			if (((is_ceil && dx[j] < 0.0) || (!is_ceil && dx[j] > 0.0)) &&
			    rx1 < prog->x1) {
				prog->x1 = rx1;
			}
		}
	}
}

/**
 * We know that mx is optimal, but it is part of an entire edge. This function
 * computes the beginning and end of the edge and returns it.
 */
static linprog2d_result_t linprog2d_calculate_edge(linprog2d_data_t *prog,
                                                   double mx) {
	unsigned int i, j, if0 = 0;
	const double *dx = prog->dx, *y0 = prog->y0;
	double ry0 = -HUGE_VAL;

	/* Find the top-most horizontal floor constraint. This must also be the
	   top-most horizontal floor constraint at mx. This function will only be
	   called if such a constraint exists. */
	for (i = 0; i < prog->floor_len; i++) {
		j = prog->floor[i];
		if (feq_(dx[j], 0.0) && y0[j] > ry0) {
			ry0 = y0[j];
			if0 = j;
		}
	}

	/* Calculate all intersections between if0 and the ceil/floor constraints,
	   update prog->x0, prog->x1 accordingly */
	linprog2d_calculate_edge_intersections(prog, prog->ceil, prog->ceil_len,
	                                       if0, mx, TRUE);
	linprog2d_calculate_edge_intersections(prog, prog->floor, prog->floor_len,
	                                       if0, mx, FALSE);

	/* Check whether the result is just a point on the edge */
	if (feq_(prog->x0, prog->x1)) {
		return linprog2d_result_point(&(prog->R), &(prog->o), prog->x0, ry0);
	} else {
		/* Return the actual edge */
		return linprog2d_result_edge(&(prog->R), &(prog->o), prog->x0, ry0,
		                             prog->x1, ry0);
	}
}

/**
 * Calculates the optimal point for a single remaining floor and ceil
 * constraint. This is the last step in the linprog2d_solve() function.
 */
static linprog2d_result_t linprog2d_calculate_result(linprog2d_data_t *prog) {
	/* Aliases */
	const unsigned int ic0 = prog->ceil[0], if0 = prog->floor[0];
	const double *Gx = prog->Gx, *Gy = prog->Gy, *h = prog->h;
	const double *dx = prog->dx, *y0 = prog->y0;
	double x0 = prog->x0, x1 = prog->x1, ry0, ry1;

	/* There is no floor constraint. The problem is unbounded. */
	if (prog->floor_len == 0U) {
		return linprog2d_result_unbounded();
	}

	/* If there is a single ceiling constraint left, compute the intersection
	   point with the floor constraint and adapt the left or right bound
	   accordingly. */
	if (prog->ceil_len > 0U) {
		double ix, iy;
		if (linprog2d_calculate_intersect(Gx[ic0], Gy[ic0], h[ic0], Gx[if0],
		                                  Gy[if0], h[if0], &ix, &iy)) {
			if (dx[if0] > dx[ic0]) {
				x1 = fmin_(x1, ix); /* optimum is on the left side, move x1 */
			} else {
				x0 = fmax_(x0, ix); /* optimum is on the right side, move x0 */
			}
		} else {
			/* Ceil and floor are parallel. Abort if problem is not feasible. */
			if (!feq_(y0[if0], y0[ic0]) && y0[if0] > y0[ic0]) {
				/* y-offset of floor constraint above the ceiling constraint. */
				return linprog2d_result_infeasible();
			}
		}
	}

	/* Return the lowest point on the remaining floor constraint. */
	ry0 = y0[if0] + x0 * dx[if0], ry1 = y0[if0] + x1 * dx[if0];
	if (feq_(dx[if0], 0.0)) { /* Floor is horizontal. Result may be a line. */
		if (x0 > -HUGE_VAL && x1 < HUGE_VAL) {
			/* Result is a line. Return this line. */
			return linprog2d_result_edge(&(prog->R), &(prog->o), x0, ry0, x1,
			                             ry1);
		} else {
			return linprog2d_result_unbounded();
		}
	} else if (dx[if0] > 0.0) { /* Minimum is on the left */
		if (x0 <= -HUGE_VAL) {
			return linprog2d_result_unbounded();
		}
		return linprog2d_result_point(&(prog->R), &(prog->o), x0, ry0);
	} else /* if (dx[if0] < 0.0) */ { /* Minimum is on the right */
		if (x1 >= HUGE_VAL) {
			return linprog2d_result_unbounded();
		}
		return linprog2d_result_point(&(prog->R), &(prog->o), x1, ry1);
	}
}

/******************************************************************************
 * EXTERNAL API                                                               *
 ******************************************************************************/

linprog2d_t *linprog2d_init(unsigned int capacity, char *mem) {
	return linprog2d_init_internal((linprog2d_data_t *)mem, capacity,
	                               mem + sizeof(linprog2d_data_t));
}

linprog2d_result_t linprog2d_solve(linprog2d_t *prog_, double cx, double cy,
                                   const double *Gx, const double *Gy,
                                   const double *h, unsigned int n) {
	double x = 0.0, y = 0.0; /* result x, y */
	linprog2d_data_t *prog = (linprog2d_data_t *)prog_;
	bool_t optimum_is_left = FALSE, has_median = FALSE;

	/* Make sure the given linprog2d instance has sufficient memory to solve
	   the problem. If not, return with an error. */
	if (!prog || prog->capacity < n) {
		return linprog2d_result_err();
	}

	/* Copy the problem to the program storage and condition it. */
	linprog2d_reset(prog, n);
	linprog2d_condition_problem(prog, cx, cy, Gx, Gy, h);

	/* Categorize the constraints into ceil, floor, and vertical constraints. */
	if (!linprog2d_categorize_constraints(prog)) {
		return linprog2d_result_infeasible();
	}

	/* Calculate the slope for the ceil and floor constraints */
	linprog2d_calculate_yoffset_form(prog->ceil, prog->ceil_len, prog->Gx,
	                                 prog->Gy, prog->h, prog->dx, prog->y0);
	linprog2d_calculate_yoffset_form(prog->floor, prog->floor_len, prog->Gx,
	                                 prog->Gy, prog->h, prog->dx, prog->y0);

	/* Repeat until there is at most one floor and ceil constraint left or the
	   left and right bounds are invalid. */
	while ((prog->floor_len != 0U) &&
	       (prog->floor_len > 1U || prog->ceil_len > 1U) &&
	       ((prog->x1 > prog->x0) || feq_(prog->x1, prog->x0))) {
		/* Calculate constraint intersection points. Of those constraints that
		   are parallel or have an intersection point outside of [x0, x1], throw
		   one away. Furthermore, if we calculated a median in the last round
		   and know its location w.r.t. the optimum, check whether the
		   intersection point is on that median. Note that the two functions
		   below edit the ceil and floor list inplace. */
		prog->intersect_len = 0U; /* number of intersections */
		linprog2d_calculate_intersects(prog, prog->ceil, &(prog->ceil_len),
		                               TRUE, has_median, x, optimum_is_left);
		linprog2d_calculate_intersects(prog, prog->floor, &(prog->floor_len),
		                               FALSE, has_median, x, optimum_is_left);

		/* If we have no intersections, then the above code must have eliminated
		   some constraints. This will give us new pairs to try. */
		if (prog->intersect_len == 0U) {
			continue;
		}

		/* Compute the median of the x-coordinates of the intersection points
		   and update the left/right boundary. */
		x = median(prog->x_intersect, prog->intersect_len);
		switch (linprog2d_locate_optimum(prog, x, &y)) {
			case LOC_INFEASIBLE:
				return linprog2d_result_infeasible();
			case LOC_LEFT:
				prog->x1 = fmin_(prog->x1, x);
				optimum_is_left = TRUE;
				has_median = TRUE;
				break;
			case LOC_RIGHT:
				prog->x0 = fmax_(prog->x0, x);
				optimum_is_left = FALSE;
				has_median = TRUE;
				break;
			case LOC_HERE:
				return linprog2d_result_point(&prog->R, &prog->o, x, y);
			case LOC_HERE_EDGE:
				return linprog2d_calculate_edge(prog, x);
		}
	}

	/* Compute the results from the remaining floor and ceil constraint */
	return linprog2d_calculate_result(prog);
}

#ifndef LINPROG2D_REDUCED_INTERFACE
linprog2d_size_t linprog2d_mem_size(unsigned int capacity) {
	linprog2d_size_t res = 0UL;

	/* Main datastructure plus alignment */
	res += sizeof(linprog2d_data_t) + 64UL;

	/* Space for the Gx, Gy, h, dx, y0, x_intersect lists plus alignment. The
	   x_intersect list only has half the length. */
	res +=
	    (sizeof(double) * 5UL + sizeof(double) / 2UL) * capacity + 64UL * 6UL;

	/* Space for the ceil, floor, tmp lists plus alignment. */
	res += sizeof(unsigned int) * 3UL * capacity + 64UL * 3UL;

	return res;
}

linprog2d_t *linprog2d_create(unsigned int capacity) {
#ifndef LINPROG2D_NO_ALLOC
	return linprog2d_init(capacity,
	                      (char *)malloc(linprog2d_mem_size(capacity)));
#else
	return NULL;
#endif
}

void linprog2d_free(linprog2d_t *prog) {
#ifndef LINPROG2D_NO_ALLOC
	free(prog); /* Free the previously allocated memory */
#endif
}

unsigned int linprog2d_capacity(const linprog2d_t *prog) {
	return ((linprog2d_data_t *)prog)->capacity;
}

linprog2d_result_t linprog2d_solve_simple(double cx, double cy,
                                          const double *Gx, const double *Gy,
                                          const double *h, unsigned int n) {
#ifndef LINPROG2D_NO_ALLOC
	linprog2d_t *prog = linprog2d_create(n);
	if (prog) {
		linprog2d_result_t res = linprog2d_solve(prog, cx, cy, Gx, Gy, h, n);
		linprog2d_free(prog);
		return res;
	}
#endif /* LINPROG2D_NO_ALLOC */
	return linprog2d_result_err();
}
#endif /* LINPROG2D_REDUCED_INTERFACE */
