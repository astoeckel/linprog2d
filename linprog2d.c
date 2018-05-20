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
#include <stdlib.h>

/******************************************************************************
 * PRIVATE HELPER FUNCTIONS                                                   *
 ******************************************************************************/

#define TRUE 1
#define FALSE 0

/**
 * Implementation of the C99 fmax.
 */
static double fmax_(double x, double y) { return (x >= y) ? x : y; }

/**
 * Implementation of the C99 fmin.
 */
static double fmin_(double x, double y) { return (x <= y) ? x : y; }

/**
 * Implementation of the C99 hypot.
 */
static double hypot_(double x, double y) { return sqrt(x * x + y * y); }

/**
 * Helper function; aligns the given pointer at the next 64-byte boundary; on
 * most systems 64 bytes is the size of a cache-line.
 */
static void *mem_align64(void *p, unsigned int offs)
{
	return (void *)(((unsigned long int)p + offs + 63) & (~64UL));
}

/**
 * Vector in 2D space.
 */
struct vec2 {
	double x, y;
};

/**
 * Returns a vec2 instance that is filled with the given values.
 */
static struct vec2 vec2_create(double x, double y)
{
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
static struct mat22 mat22_create(double a11, double a12, double a21, double a22)
{
	struct mat22 res;
	res.a11 = a11;
	res.a12 = a12;
	res.a21 = a21;
	res.a22 = a22;
	return res;
}

/**
 * Returns a rotation matrix that aligns the given (x, y) vector with the
 * y axis.
 */
static struct mat22 mat22_rot(double x, double y)
{
	const double h = hypot_(x, y);
	return mat22_create(-y / h, -x / h, x / h, -y / h);
}

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

/**
 * Typedef allowing to easily
 */
typedef struct linprog2d_data linprog2d_data_t;

/**
 * Computes the value the entire expression a, b, c should be divided by in
 * order to get a unique constraint.
 */
static double linprog2d_normalization_coeff(double a, double b, double c)
{
	a = fabs(a), b = fabs(b), c = fabs(c);
	return fmax_(fmax_(a, b), fmax_(a, c));
}

/**
 * Takes the linear program provided by the user and rotates it such that the
 * gradient is aligned with the y-axis. Additionally, normalises all constraints
 * such that the maximum coefficient is one in each line. Furthermore shifts the
 * entire problem space such that all constraints are centered around the orign.
 * This centering is performed by finding an offset vector (o.x, o.y) s.t.
 * the following expression is minimized:
 *
 * sum_i (h_i - Gx_i * o.x - Gy_i * o.y)^2
 *
 * The closed-form solution to this least-squares optimization problem is
 *
 * o = (G.T * G)^-1 * G.T * h
 *
 * Which can be computed in linear time, since G.T * G is a 2x2 matrix.
 */
static int linprog2d_condition_problem(double cx, double cy,
                                       const double *src_Gx,
                                       const double *src_Gy,
                                       const double *src_h,
                                       linprog2d_data_t *prog)
{
	struct mat22 R = mat22_rot(cx, cy);
	struct vec2 o = vec2_create(0.0, 0.0);               /* Offset vector */
	struct mat22 GTG = mat22_create(0.0, 0.0, 0.0, 0.0); /* Matrix G.T G */
	struct vec2 GTc = vec2_create(0.0, 0.0);             /* Vector G.T c */
	double Gx, Gy, h, norm, GTG_det;                     /* Temp variables */
	unsigned int i_tar = 0, i_src = 0;
	double *tar_Gx = prog->Gx, *tar_Gy = prog->Gy, *tar_h = prog->h;

	/* Copy the memory from the source to the target location; rotate all the
	   source vectors. At the same time normalise the problem such that the
	   coefficient with the largest absolute value is scaled to +-1. */
	for (i_src = 0; i_src < prog->n; i_src++) {
		/* Rotate the constraint direction on the left-hand side */
		Gx = R.a11 * src_Gx[i_src] + R.a12 * src_Gy[i_src];
		Gy = R.a21 * src_Gx[i_src] + R.a22 * src_Gy[i_src];
		h = src_h[i_src];

		/* Skip invalid constraints */
		if (Gx == 0.0 && Gy == 0.0) {
			if (h <= 0.0) {
				/* Constraint of the form 0 >= h is always true for h <= 0.0 */
				continue;
			}
			else {
				/* This constraint is always false. Abort. */
				return FALSE;
			}
		}

		/* Normalize the constraints by dividing both the right- and left-hand
		   side by the largest coefficient. Due to the above check we know that
		   we will never divide by zero. */
		norm = linprog2d_normalization_coeff(Gx, Gy, h);
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

	/* Invert the GTG matrix (if possible) and compute o */
	GTG_det = GTG.a11 * GTG.a22 - GTG.a12 * GTG.a12;
	if (GTG_det != 0.0) {
		o.x = (GTG.a22 * GTc.x - GTG.a12 * GTc.y) / GTG_det;
		o.y = (-GTG.a12 * GTc.x + GTG.a22 * GTc.y) / GTG_det;
	}

	/* Update the linear program data */
	prog->n = i_tar; /* Constraints may have been eliminated */
	prog->R = R;
	prog->o = o;

	return TRUE; /* Success */
}

/******************************************************************************
 * EXTERNAL API                                                               *
 ******************************************************************************/

unsigned int linprog2d_mem_size(unsigned int capacity)
{
	return sizeof(linprog2d_data_t) + sizeof(double) * 3 * capacity + 64 * 3;
}

linprog2d_t *linprog2d_init(unsigned int capacity, char *mem)
{
	linprog2d_data_t *prog = (linprog2d_data_t *)mem;
	if (prog) {
		prog->Gx = mem_align64(mem, 0U);
		prog->Gy = mem_align64(prog->Gx, sizeof(double) * capacity);
		prog->h = mem_align64(prog->Gy, sizeof(double) * capacity);
		prog->capacity = capacity;
		prog->R = mat22_create(0.0, 0.0, 0.0, 0.0);
		prog->o = vec2_create(0.0, 0.0);
		prog->n = 0U;
	}
	return mem;
}

linprog2d_t *linprog2d_create(unsigned int capacity)
{
	return linprog2d_init(capacity, malloc(linprog2d_mem_size(capacity)));
}

void linprog2d_free(linprog2d_t *prog)
{
	free(prog); /* Free the previously allocated memory */
}

unsigned int linprog2d_capacity(const linprog2d_t *prog)
{
	return ((linprog2d_data_t*)prog)->capacity;
}

linprog2d_result_t linprog2d_solve(linprog2d_t *prog, double cx, double cy,
                                   const double *Gx, const double *Gy,
                                   const double *h, unsigned int n)
{
	/* Make sure the given linprog2d instance has sufficient memory to solve
	   the problem. If not, return with an error. */
	if (!prog || linprog2d_capacity(prog) < n) {
		/*		return linprog2d_result_t{LP2D_ERROR, 0.0, 0.0, 0.0, 0.0}; */
	}
	/*	return (linprog2d_d*)prog->solve(cx, cy, Gx, Gy, h,
	                                                            n);*/
}

linprog2d_result_t linprog2d_solve_simple(double cx, double cy,
                                          const double *Gx, const double *Gy,
                                          const double *h, unsigned int n)
{
	linprog2d_t *prog = linprog2d_create(n);
	if (prog) {
		linprog2d_result_t res = linprog2d_solve(prog, cx, cy, Gx, Gy, h, n);
		linprog2d_free(prog);
		return res;
	}
	else {
		/*		return linprog2d_result_t{LP2D_ERROR, 0.0, 0.0, 0.0, 0.0}; */
	}
}

/******************************************************************************
 * UNIT TESTS                                                                 *
 ******************************************************************************/

#ifdef LINPROG_2D_TEST

/*
 * The following code is only compiled if the LINPROG_2D_TEST prepocessor flag
 * exists. In this case, this translation unit exports a main function that
 * executes the unit tests defined below. I decided to move unit tests here so
 * I can test private functions without having to expose them through the public
 * interface.
 */

#include <setjmp.h>
#include <stdio.h>

/******************************************************************************
 * Minimal, yet nicely coloured unit testing framework                        *
 ******************************************************************************/

static volatile int failed = FALSE;
static jmp_buf assert_buffer;

#define ANSI_RED "\33[31;1m"
#define ANSI_GRAY "\33[37;2m"
#define ANSI_GREEN "\33[32;1m"
#define ANSI_RESET "\33[0m"

#define STR_DETAIL(X) #X
#define STR(X) STR_DETAIL(X)

#define EXPECT(should, is, rel)                                       \
	do {                                                              \
		if (!((should)rel(is))) {                                     \
			fprintf(stderr, ANSI_RED "[ERR]" ANSI_RESET               \
			                         " Assertion failed in " __FILE__ \
			                         ", line " STR(__LINE__) "\n");   \
			failed = TRUE;                                            \
		}                                                             \
	} while (0)

#define ASSERT(should, is, rel)        \
	do {                               \
		EXPECT(should, is, rel);       \
		if (failed) {                  \
			longjmp(assert_buffer, 1); \
		}                              \
	} while (0)

#define RUN(test)                                                              \
	do {                                                                       \
		int old_failed = failed;                                               \
		failed = FALSE;                                                        \
		fprintf(stderr,                                                        \
		        ANSI_GRAY "----> " ANSI_RESET "Running test \"" #test "\"\n"); \
		if (!setjmp(assert_buffer)) {                                          \
			test();                                                            \
		}                                                                      \
		if (failed) {                                                          \
			fprintf(stderr, ANSI_RED "[ERR] " ANSI_RESET                       \
			                         "Failed test \"" #test "\"\n");           \
		}                                                                      \
		else {                                                                 \
			fprintf(stderr, ANSI_GREEN "[OK!] " ANSI_RESET "Test \"" #test     \
			                           "\" successful\n");                     \
		}                                                                      \
		failed = failed || old_failed;                                         \
	} while (0)

#define EXPECT_EQ(should, is) EXPECT(should, is, ==)
#define ASSERT_EQ(should, is) ASSERT(should, is, ==)
#define EXPECT_GT(should, is) EXPECT(should, is, >)
#define ASSERT_GT(should, is) ASSERT(should, is, >)
#define EXPECT_GE(should, is) EXPECT(should, is, >=)
#define ASSERT_GE(should, is) ASSERT(should, is, >=)
#define EXPECT_LT(should, is) EXPECT(should, is, <)
#define ASSERT_LT(should, is) ASSERT(should, is, <)
#define EXPECT_LE(should, is) EXPECT(should, is, <=)
#define ASSERT_LE(should, is) ASSERT(should, is, <=)
#define EXPECT_NE(should, is) EXPECT(should, is, !=)
#define ASSERT_NE(should, is) ASSERT(should, is, !=)
#define EXPECT_NEAR(should, is, eps) EXPECT(eps, fabs(should - is), >=)
#define ASSERT_NEAR(should, is, eps) ASSERT(eps, fabs(should - is), >=)

/******************************************************************************
 * Actual unit tests                                                          *
 ******************************************************************************/

void test_linprog2d_normalization_coeff()
{
	EXPECT_EQ(0.0, linprog2d_normalization_coeff(0.0, 0.0, 0.0));
	EXPECT_EQ(1.0, linprog2d_normalization_coeff(1.0, 0.0, 0.0));
	EXPECT_EQ(1.0, linprog2d_normalization_coeff(0.0, 1.0, 0.0));
	EXPECT_EQ(1.0, linprog2d_normalization_coeff(0.0, 0.0, 1.0));
	EXPECT_EQ(1.0, linprog2d_normalization_coeff(-1.0, 0.0, 0.0));
	EXPECT_EQ(1.0, linprog2d_normalization_coeff(0.0, -1.0, 0.0));
	EXPECT_EQ(1.0, linprog2d_normalization_coeff(0.0, 0.0, -1.0));
	EXPECT_EQ(2.0, linprog2d_normalization_coeff(0.0, -2.0, 1.0));
	EXPECT_EQ(3.0, linprog2d_normalization_coeff(3.0, -2.0, 1.0));
	EXPECT_EQ(3.0, linprog2d_normalization_coeff(-2.0, 3.0, -1.0));
	EXPECT_EQ(3.0, linprog2d_normalization_coeff(2.0, 1.0, -3.0));
}

void test_linprog2d_create_and_capacity()
{
	{
		linprog2d_t *prog = linprog2d_create(128U);
		ASSERT_NE(NULL, prog);
		EXPECT_EQ(128U, linprog2d_capacity(prog));
		linprog2d_free(prog);
	}

	{
		linprog2d_t *prog = linprog2d_create(0U);
		ASSERT_NE(NULL, prog);
		EXPECT_EQ(0U, linprog2d_capacity(prog));
		linprog2d_free(prog);
	}
}

void test_linprog2d_problem_too_large()
{
	linprog2d_t *prog = linprog2d_create(128U);
	ASSERT_EQ(LP2D_ERROR,
	          linprog2d_solve(prog, 0.0, 0.0, NULL, NULL, NULL, 129U).status);
	linprog2d_free(prog);
}

/******************************************************************************
 * Main function                                                              *
 ******************************************************************************/

int main()
{
	RUN(test_linprog2d_normalization_coeff);
	RUN(test_linprog2d_create_and_capacity);
	RUN(test_linprog2d_problem_too_large);
	return failed ? 1 : 0;
}

#endif
