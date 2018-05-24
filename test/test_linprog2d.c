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
 * @file test/test_linprog2d.c
 *
 * Tests the two-dimensional linear programming algorithm including its internal
 * code and data structures.
 *
 * @author Andreas Stöckel
 */

/* We're testing all internals as well, so directly include the source code */
#include "../linprog2d.c"

/******************************************************************************
 * Minimal, yet nicely coloured unit testing framework                        *
 ******************************************************************************/

#include <setjmp.h> /* Jikes. Required as an exception replacement in ASSERT. */
#include <stdio.h>

static volatile int n_failed = 0;
static volatile int n_success = 0;
static volatile bool_t failed = FALSE;
static jmp_buf assert_buffer;

#define ANSI_RED "\33[38;5;166;1m"
#define ANSI_GRAY "\33[37;2m"
#define ANSI_GREEN "\33[38;5;40;1m"
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
		failed = FALSE;                                                        \
		fprintf(stderr,                                                        \
		        ANSI_GRAY "---->" ANSI_RESET " Running test \"" #test "\"\n"); \
		if (!setjmp(assert_buffer)) {                                          \
			test();                                                            \
		}                                                                      \
		if (failed) {                                                          \
			n_failed++;                                                        \
			fprintf(stderr, ANSI_RED "[ERR]" ANSI_RESET " Test \"" #test       \
			                         "\" failed!\n");                          \
		} else {                                                               \
			n_success++;                                                       \
			fprintf(stderr, ANSI_GREEN "[OK!]" ANSI_RESET " Test \"" #test     \
			                           "\" successful\n");                     \
		}                                                                      \
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
#define EXPECT_TRUE(is) EXPECT(TRUE, (!!(is)) ? TRUE : FALSE, ==)
#define ASSERT_TRUE(is) ASSERT(TRUE, (!!(is)) ? TRUE : FALSE, ==)
#define EXPECT_FALSE(is) EXPECT(FALSE, (!!(is)) ? TRUE : FALSE, ==)
#define ASSERT_FALSE(is) ASSERT(FALSE, (!!(is)) ? TRUE : FALSE, ==)
#define EXPECT_NEAR(should, is) EXPECT_TRUE(feq_(should, is))
#define ASSERT_NEAR(should, is) ASSERT_TRUE(feq_(should, is))

/******************************************************************************
 * Actual unit tests                                                          *
 ******************************************************************************/

void test_feq() {
	EXPECT_TRUE(feq_(0.0, -0.0));
	EXPECT_TRUE(feq_(0.0, 0.0));
	EXPECT_TRUE(feq_(0.0, 1e-31));
	EXPECT_FALSE(feq_(0.0, 1.0));
	EXPECT_FALSE(feq_(1e-15, 1.01e-15));
	EXPECT_TRUE(feq_(1.0, 1.0 + 0.1e-15));
	EXPECT_FALSE(feq_(0.0, -1.0));
	EXPECT_FALSE(feq_(-1e-15, -1.01e-15));
	EXPECT_TRUE(feq_(-1.0, -1.0 + 0.1e-15));
	EXPECT_TRUE(feq_(1e-15, 1e-15 + 0.1e-30));
}

void test_memalign64() {
	unsigned long int i;
	EXPECT_EQ((void *)0x0, mem_align64((void *)0x0, 0U));
	for (i = 1; i <= 64; i++) {
		EXPECT_EQ((void *)0x40, mem_align64((void *)i, 0U));
	}
	for (i = 65; i < 128; i++) {
		EXPECT_EQ((void *)0x80, mem_align64((void *)i, 0U));
	}
}

void test_sort() {
	double d[5];

	sort(d, 0); /* Should do nothing. Check valgrind log. */

	d[0] = 1;
	sort(d, 1);
	EXPECT_EQ(1.0, d[0]);

	d[1] = -1;
	sort(d, 2);
	EXPECT_EQ(-1.0, d[0]);
	EXPECT_EQ(1.0, d[1]);

	d[0] = -2;
	sort(d, 2);
	EXPECT_EQ(-2.0, d[0]);
	EXPECT_EQ(1.0, d[1]);

	d[0] = 4;
	d[1] = 3;
	d[2] = 7;
	sort(d, 3);
	EXPECT_EQ(3.0, d[0]);
	EXPECT_EQ(4.0, d[1]);
	EXPECT_EQ(7.0, d[2]);

	d[0] = 7;
	d[1] = 3;
	d[2] = 4;
	sort(d, 3);
	EXPECT_EQ(3.0, d[0]);
	EXPECT_EQ(4.0, d[1]);
	EXPECT_EQ(7.0, d[2]);

	d[0] = 4;
	d[1] = 7;
	d[2] = 3;
	sort(d, 3);
	EXPECT_EQ(3.0, d[0]);
	EXPECT_EQ(4.0, d[1]);
	EXPECT_EQ(7.0, d[2]);

	d[0] = 4;
	d[1] = 7;
	d[2] = 3;
	d[3] = 5;
	sort(d, 4);
	EXPECT_EQ(3.0, d[0]);
	EXPECT_EQ(4.0, d[1]);
	EXPECT_EQ(5.0, d[2]);
	EXPECT_EQ(7.0, d[3]);

	d[0] = 3;
	d[1] = 5;
	d[2] = 4;
	d[3] = 7;
	sort(d, 4);
	EXPECT_EQ(3.0, d[0]);
	EXPECT_EQ(4.0, d[1]);
	EXPECT_EQ(5.0, d[2]);
	EXPECT_EQ(7.0, d[3]);

	d[0] = 3;
	d[1] = 5;
	d[2] = 4;
	d[3] = 7;
	d[4] = 1;
	sort(d, 5);
	EXPECT_EQ(1.0, d[0]);
	EXPECT_EQ(3.0, d[1]);
	EXPECT_EQ(4.0, d[2]);
	EXPECT_EQ(5.0, d[3]);
	EXPECT_EQ(7.0, d[4]);
}

void test_partition() {
	{
		double d[20] = {5,  13, 13, 8, 9, 12, 19, 2,  1,  13,
		                14, 10, 6,  3, 2, 3,  7,  16, 17, 16};

		EXPECT_EQ(8U, partition(d, 20, 8.0));
		EXPECT_EQ(d[8U], 8.0);

		EXPECT_EQ(12U, partition(d, 20, 13.0));
		EXPECT_EQ(d[12U], 13.0);

		EXPECT_EQ(1U, partition(d, 20, 2.0));
		EXPECT_EQ(d[1U], 2.0);

		EXPECT_EQ(0U, partition(d, 20, 1.0));
		EXPECT_EQ(d[0U], 1.0);

		EXPECT_EQ(19U, partition(d, 20, 19.0));
		EXPECT_EQ(d[19U], 19.0);
	}

	{
		double d[8] = {6, 4, 16, 7, 1, 6, 6, 14};

		EXPECT_EQ(2U, partition(d, 8, 6.0));
		EXPECT_EQ(d[2U], 6.0);

		EXPECT_EQ(0U, partition(d, 8, 1.0));
		EXPECT_EQ(d[0U], 1.0);

		EXPECT_EQ(7U, partition(d, 8, 16.0));
		EXPECT_EQ(d[7U], 16.0);
	}

	{
		double d[51] = {4, 15, 1,  3,  16, 0,  9, 0,  8, 11, 14, 15, 12,
		                8, 13, 10, 17, 7,  17, 7, 19, 2, 19, 19, 11, 10,
		                8, 7,  5,  19, 10, 18, 6, 12, 2, 9,  10, 18, 2,
		                5, 8,  6,  19, 7,  5,  9, 17, 1, 5,  2,  12};

		EXPECT_EQ(14U, partition(d, 51, 6.0));
		EXPECT_EQ(d[14U], 6.0);

		EXPECT_EQ(46U, partition(d, 51, 19.0));
		EXPECT_EQ(d[46U], 19.0);

		EXPECT_EQ(0U, partition(d, 51, 0.0));
		EXPECT_EQ(d[0U], 0.0);

		EXPECT_EQ(20U, partition(d, 51, 8.0));
		EXPECT_EQ(d[20U], 8.0);
	}

	{
		double d[1] = {6};

		EXPECT_EQ(0U, partition(d, 1, 6.0));
		EXPECT_EQ(d[0U], 6.0);
	}
}

void test_kth_smallest() {
	{
		double d[1] = {2.2};

		EXPECT_EQ(2.2, kth_smallest(d, 1U, 0U));
	}

	{
		double d[2] = {3.3, 2.1};

		EXPECT_EQ(2.1, kth_smallest(d, 2U, 0U));
		EXPECT_EQ(3.3, kth_smallest(d, 2U, 1U));
	}

	{
		double d[3] = {3.3, 2.1, 4.4};

		EXPECT_EQ(2.1, kth_smallest(d, 3U, 0U));
		EXPECT_EQ(3.3, kth_smallest(d, 3U, 1U));
		EXPECT_EQ(4.4, kth_smallest(d, 3U, 2U));
	}

	{
		double d[4] = {3.3, 2.1, 4.4, 5.2};

		EXPECT_EQ(2.1, kth_smallest(d, 4U, 0U));
		EXPECT_EQ(3.3, kth_smallest(d, 4U, 1U));
		EXPECT_EQ(4.4, kth_smallest(d, 4U, 2U));
		EXPECT_EQ(5.2, kth_smallest(d, 4U, 3U));
	}

	{
		double d[5] = {3.3, 2.1, 4.4, 5.2, 1.2};

		EXPECT_EQ(1.2, kth_smallest(d, 5U, 0U));
		EXPECT_EQ(2.1, kth_smallest(d, 5U, 1U));
		EXPECT_EQ(3.3, kth_smallest(d, 5U, 2U));
		EXPECT_EQ(4.4, kth_smallest(d, 5U, 3U));
		EXPECT_EQ(5.2, kth_smallest(d, 5U, 4U));
	}

	{
		double d[6] = {3.3, 2.1, 4.4, 5.2, 1.2, 2.3};

		EXPECT_EQ(1.2, kth_smallest(d, 6U, 0U));
		EXPECT_EQ(2.1, kth_smallest(d, 6U, 1U));
		EXPECT_EQ(2.3, kth_smallest(d, 6U, 2U));
		EXPECT_EQ(3.3, kth_smallest(d, 6U, 3U));
		EXPECT_EQ(4.4, kth_smallest(d, 6U, 4U));
		EXPECT_EQ(5.2, kth_smallest(d, 6U, 5U));
	}

	{
		double d[51] = {4, 15, 1,  3,  16, 0,  9, 0,  8, 11, 14, 15, 12,
		                8, 13, 10, 17, 7,  17, 7, 19, 2, 19, 19, 11, 10,
		                8, 7,  5,  19, 10, 18, 6, 12, 2, 9,  10, 18, 2,
		                5, 8,  6,  19, 7,  5,  9, 17, 1, 5,  2,  12};

		EXPECT_EQ(0.0, kth_smallest(d, 51U, 0U));
		EXPECT_EQ(0.0, kth_smallest(d, 51U, 1U));
		EXPECT_EQ(1.0, kth_smallest(d, 51U, 2U));
		EXPECT_EQ(1.0, kth_smallest(d, 51U, 3U));
		EXPECT_EQ(2.0, kth_smallest(d, 51U, 4U));
		EXPECT_EQ(2.0, kth_smallest(d, 51U, 5U));
		EXPECT_EQ(2.0, kth_smallest(d, 51U, 6U));
		EXPECT_EQ(2.0, kth_smallest(d, 51U, 7U));
		EXPECT_EQ(3.0, kth_smallest(d, 51U, 8U));
		EXPECT_EQ(4.0, kth_smallest(d, 51U, 9U));
		EXPECT_EQ(5.0, kth_smallest(d, 51U, 10U));
		EXPECT_EQ(5.0, kth_smallest(d, 51U, 11U));
		EXPECT_EQ(5.0, kth_smallest(d, 51U, 12U));
		EXPECT_EQ(5.0, kth_smallest(d, 51U, 13U));
		EXPECT_EQ(6.0, kth_smallest(d, 51U, 14U));
		EXPECT_EQ(6.0, kth_smallest(d, 51U, 15U));
		EXPECT_EQ(7.0, kth_smallest(d, 51U, 16U));
		EXPECT_EQ(7.0, kth_smallest(d, 51U, 17U));
		EXPECT_EQ(7.0, kth_smallest(d, 51U, 18U));
		EXPECT_EQ(7.0, kth_smallest(d, 51U, 19U));
		EXPECT_EQ(8.0, kth_smallest(d, 51U, 20U));
		EXPECT_EQ(8.0, kth_smallest(d, 51U, 21U));
		EXPECT_EQ(8.0, kth_smallest(d, 51U, 22U));
		EXPECT_EQ(8.0, kth_smallest(d, 51U, 23U));
		EXPECT_EQ(9.0, kth_smallest(d, 51U, 24U));
		EXPECT_EQ(9.0, kth_smallest(d, 51U, 25U));
		EXPECT_EQ(9.0, kth_smallest(d, 51U, 26U));
		EXPECT_EQ(10.0, kth_smallest(d, 51U, 27U));
		EXPECT_EQ(10.0, kth_smallest(d, 51U, 28U));
		EXPECT_EQ(10.0, kth_smallest(d, 51U, 29U));
		EXPECT_EQ(10.0, kth_smallest(d, 51U, 30U));
		EXPECT_EQ(11.0, kth_smallest(d, 51U, 31U));
		EXPECT_EQ(11.0, kth_smallest(d, 51U, 32U));
		EXPECT_EQ(12.0, kth_smallest(d, 51U, 33U));
		EXPECT_EQ(12.0, kth_smallest(d, 51U, 34U));
		EXPECT_EQ(12.0, kth_smallest(d, 51U, 35U));
		EXPECT_EQ(13.0, kth_smallest(d, 51U, 36U));
		EXPECT_EQ(14.0, kth_smallest(d, 51U, 37U));
		EXPECT_EQ(15.0, kth_smallest(d, 51U, 38U));
		EXPECT_EQ(15.0, kth_smallest(d, 51U, 39U));
		EXPECT_EQ(16.0, kth_smallest(d, 51U, 40U));
		EXPECT_EQ(17.0, kth_smallest(d, 51U, 41U));
		EXPECT_EQ(17.0, kth_smallest(d, 51U, 42U));
		EXPECT_EQ(17.0, kth_smallest(d, 51U, 43U));
		EXPECT_EQ(18.0, kth_smallest(d, 51U, 44U));
		EXPECT_EQ(18.0, kth_smallest(d, 51U, 45U));
		EXPECT_EQ(19.0, kth_smallest(d, 51U, 46U));
		EXPECT_EQ(19.0, kth_smallest(d, 51U, 47U));
		EXPECT_EQ(19.0, kth_smallest(d, 51U, 48U));
		EXPECT_EQ(19.0, kth_smallest(d, 51U, 49U));
		EXPECT_EQ(19.0, kth_smallest(d, 51U, 50U));
	}
}

void test_median() {
	{
		double d[1] = {1.2};
		EXPECT_EQ(1.2, median(d, 1U));
	}

	{
		double d[2] = {1.2, 2.4};
		EXPECT_EQ(2.4, median(d, 2U));
	}

	{
		double d[3] = {3.5, 1.2, 2.4};
		EXPECT_EQ(2.4, median(d, 3U));
	}

	{
		double d[4] = {6.8, 3.5, 1.2, 2.4};
		EXPECT_EQ(3.5, median(d, 4U));
	}

	{
		double d[5] = {6.8, 2.9, 3.5, 1.2, 2.4};
		EXPECT_EQ(2.9, median(d, 5U));
	}

	{
		double d[6] = {6.8, 5.6, 2.9, 3.5, 1.2, 2.4};
		EXPECT_EQ(3.5, median(d, 6U));
	}

	{
		double d[7] = {7.0, 6.8, 5.6, 2.9, 3.5, 1.2, 2.4};
		EXPECT_EQ(3.5, median(d, 7U));
	}
	{
		double d[51] = {4, 15, 1,  3,  16, 0,  9, 0,  8, 11, 14, 15, 12,
		                8, 13, 10, 17, 7,  17, 7, 19, 2, 19, 19, 11, 10,
		                8, 7,  5,  19, 10, 18, 6, 12, 2, 9,  10, 18, 2,
		                5, 8,  6,  19, 7,  5,  9, 17, 1, 5,  2,  12};
		EXPECT_EQ(9.0, median(d, 51U));
	}

	{
		double d[193] = {
		    56, 77, 40, 23, 40, 20, 76, 17, 69, 29, 84, 1,  4,  27, 43, 55, 60,
		    3,  73, 0,  15, 61, 1,  21, 78, 47, 22, 19, 94, 67, 78, 83, 47, 45,
		    2,  98, 17, 63, 44, 44, 81, 62, 53, 86, 65, 15, 21, 39, 53, 72, 51,
		    63, 28, 54, 29, 2,  69, 83, 68, 86, 8,  32, 6,  43, 45, 62, 6,  60,
		    2,  64, 77, 28, 67, 31, 59, 1,  63, 46, 31, 67, 51, 31, 45, 47, 55,
		    19, 98, 14, 38, 73, 44, 94, 84, 64, 67, 65, 70, 93, 96, 7,  6,  96,
		    53, 87, 90, 43, 56, 19, 88, 41, 75, 15, 80, 71, 26, 35, 35, 28, 65,
		    22, 30, 52, 51, 73, 24, 69, 19, 87, 7,  94, 25, 98, 32, 1,  24, 10,
		    36, 52, 80, 77, 20, 0,  37, 59, 6,  55, 31, 4,  60, 17, 13, 27, 27,
		    93, 59, 26, 45, 29, 92, 2,  78, 32, 61, 0,  79, 83, 49, 49, 67, 14,
		    76, 58, 50, 11, 2,  46, 76, 21, 66, 67, 21, 26, 50, 38, 86, 98, 3,
		    71, 92, 57, 90, 73, 82};
		EXPECT_EQ(49.0, median(d, 193U));
	}
}

void test_linprog2d_normalization_coeff() {
	EXPECT_EQ(0.0, linprog2d_normalization_coeff(0.0, 0.0));
	EXPECT_EQ(1.0, linprog2d_normalization_coeff(1.0, 0.0));
	EXPECT_EQ(1.0, linprog2d_normalization_coeff(0.0, 1.0));
	EXPECT_EQ(1.0, linprog2d_normalization_coeff(-1.0, 0.0));
	EXPECT_EQ(1.0, linprog2d_normalization_coeff(0.0, -1.0));
	EXPECT_EQ(2.0, linprog2d_normalization_coeff(1.0, -2.0));
	EXPECT_EQ(2.0, linprog2d_normalization_coeff(-2.0, -1.0));
}

void test_linprog2d_create_and_capacity() {
	{
		unsigned int i;
		linprog2d_data_t *prog = linprog2d_create(128U);
		ASSERT_NE(NULL, prog);
		EXPECT_EQ(128U, linprog2d_capacity(prog));

		/* Fill the individual lists with data */
		for (i = 0; i < 128U; i++) {
			prog->Gx[i] = 10 * i + 0;
			prog->Gy[i] = 10 * i + 1;
			prog->h[i] = 10 * i + 2;
			prog->dx[i] = 10 * i + 3;
			prog->y0[i] = 10 * i + 4;
			if (i < 64U) {
				prog->x_intersect[i] = 10 * i + 5;
			}
			prog->ceil[i] = 10 * i + 6;
			prog->floor[i] = 10 * i + 7;
			prog->tmp[i] = 10 * i + 8;
		}

		/* Try to read the data back */
		for (i = 0; i < 128U; i++) {
			EXPECT_EQ(10 * i + 0, (unsigned int)(prog->Gx[i]));
			EXPECT_EQ(10 * i + 1, (unsigned int)(prog->Gy[i]));
			EXPECT_EQ(10 * i + 2, (unsigned int)(prog->h[i]));
			EXPECT_EQ(10 * i + 3, (unsigned int)(prog->dx[i]));
			EXPECT_EQ(10 * i + 4, (unsigned int)(prog->y0[i]));
			if (i < 64U) {
				EXPECT_EQ(10 * i + 5, (unsigned int)(prog->x_intersect[i]));
			}
			EXPECT_EQ(10 * i + 6, (unsigned int)(prog->ceil[i]));
			EXPECT_EQ(10 * i + 7, (unsigned int)(prog->floor[i]));
			EXPECT_EQ(10 * i + 8, (unsigned int)(prog->tmp[i]));
		}

		linprog2d_free(prog);
	}

	{
		linprog2d_t *prog = linprog2d_create(0U);
		ASSERT_NE(NULL, prog);
		EXPECT_EQ(0U, linprog2d_capacity(prog));
		linprog2d_free(prog);
	}
}

void test_linprog2d_problem_too_large() {
	linprog2d_data_t prog;
	linprog2d_init_internal(&prog, 128U, NULL);

	EXPECT_EQ(128U, linprog2d_capacity(&prog));
	EXPECT_EQ(LP2D_ERROR,
	          linprog2d_solve(&prog, 0.0, 0.0, NULL, NULL, NULL, 129U).status);
}

void test_linprog2d_condition_problem_rotation() {
	/* Manually setup the linprog2d_data_t structure */
	linprog2d_data_t prog;
	linprog2d_reset(&prog, 0U);

	EXPECT_EQ(TRUE,
	          linprog2d_condition_problem(&prog, 0.0, 1.0, NULL, NULL, NULL));
	EXPECT_EQ(1.0, prog.R.a11);
	EXPECT_EQ(0.0, prog.R.a12);
	EXPECT_EQ(0.0, prog.R.a21);
	EXPECT_EQ(1.0, prog.R.a22);

	EXPECT_EQ(TRUE,
	          linprog2d_condition_problem(&prog, 0.0, 2.0, NULL, NULL, NULL));
	EXPECT_EQ(1.0, prog.R.a11);
	EXPECT_EQ(0.0, prog.R.a12);
	EXPECT_EQ(0.0, prog.R.a21);
	EXPECT_EQ(1.0, prog.R.a22);

	EXPECT_EQ(TRUE,
	          linprog2d_condition_problem(&prog, 1.0, 0.0, NULL, NULL, NULL));
	EXPECT_EQ(0.0, prog.R.a11);
	EXPECT_EQ(-1.0, prog.R.a12);
	EXPECT_EQ(1.0, prog.R.a21);
	EXPECT_EQ(0.0, prog.R.a22);
}

void test_linprog2d_condition_problem_eliminate_invalid() {
	linprog2d_data_t prog;
	double Gx, Gy, h;
	double Gx_tar, Gy_tar, h_tar;

	/* Manually setup the linprog2d_data_t structure */
	prog.capacity = 1U;
	prog.Gx = &Gx_tar, prog.Gy = &Gy_tar, prog.h = &h_tar;

	Gx = 0.0, Gy = 0.0, h = 0.0;
	prog.n = 1U;
	EXPECT_EQ(TRUE, linprog2d_condition_problem(&prog, 1.0, 0.0, &Gx, &Gy, &h));
	EXPECT_EQ(0U, prog.n);

	Gx = 0.0, Gy = 0.0, h = -1.0;
	prog.n = 1U;
	EXPECT_EQ(TRUE, linprog2d_condition_problem(&prog, 1.0, 0.0, &Gx, &Gy, &h));
	EXPECT_EQ(0U, prog.n);

	Gx = 0.0, Gy = 0.0, h = 1.0;
	prog.n = 1U;
	EXPECT_EQ(FALSE,
	          linprog2d_condition_problem(&prog, 1.0, 0.0, &Gx, &Gy, &h));
}

void test_linprog2d_condition_problem_offset1() {
	/* Setup a set of constraints that form a box from (3, 4) to (5, 8). This
	   box should be shifted to the origin by linprog2d_condition_problem,
	   resulting in a box from (-1, -2) to (1, 2). */
	linprog2d_data_t prog;
	double Gx[4] = {1.0, -1.0, 0.0, 0.0};
	double Gy[4] = {0.0, 0.0, 1.0, -1.0};
	double h[4] = {3.0, -5.0, 4.0, -8.0};
	double Gx_tar[4], Gy_tar[4], h_tar[4];

	/* Manually setup the linprog2d_data_t structure */
	linprog2d_reset(&prog, 4U);
	prog.Gx = Gx_tar, prog.Gy = Gy_tar, prog.h = h_tar;

	EXPECT_EQ(TRUE, linprog2d_condition_problem(&prog, 0.0, 1.0, Gx, Gy, h));

	EXPECT_EQ(1.0, Gx_tar[0]);
	EXPECT_EQ(-1.0, Gx_tar[1]);
	EXPECT_EQ(0.0, Gx_tar[2]);
	EXPECT_EQ(0.0, Gx_tar[3]);

	EXPECT_EQ(0.0, Gy_tar[0]);
	EXPECT_EQ(0.0, Gy_tar[1]);
	EXPECT_EQ(1.0, Gy_tar[2]);
	EXPECT_EQ(-1.0, Gy_tar[3]);

	EXPECT_EQ(-1.0, h_tar[0]);
	EXPECT_EQ(-1.0, h_tar[1]);
	EXPECT_EQ(-2.0, h_tar[2]);
	EXPECT_EQ(-2.0, h_tar[3]);

	EXPECT_EQ(4.0, prog.o.x);
	EXPECT_EQ(6.0, prog.o.y);
}

void test_linprog2d_condition_problem_offset2() {
	/* Setup a set of constraints that form a box rotated by 45°. The centre
	   of this box is at (4.5, 4.5). */
	linprog2d_data_t prog;
	double Gx[4] = {1.0, -1.0, 1.0, -1.0};
	double Gy[4] = {1.0, 1.0, -1.0, -1.0};
	double h[4] = {6.0, -6.0, -6.0, -12.0};
	double Gx_tar[4], Gy_tar[4], h_tar[4];

	/* Manually setup the linprog2d_data_t structure */
	linprog2d_reset(&prog, 4U);
	prog.Gx = Gx_tar, prog.Gy = Gy_tar, prog.h = h_tar;

	EXPECT_EQ(TRUE, linprog2d_condition_problem(&prog, 0.0, 1.0, Gx, Gy, h));

	EXPECT_EQ(1.0, Gx_tar[0]);
	EXPECT_EQ(-1.0, Gx_tar[1]);
	EXPECT_EQ(1.0, Gx_tar[2]);
	EXPECT_EQ(-1.0, Gx_tar[3]);

	EXPECT_EQ(1.0, Gy_tar[0]);
	EXPECT_EQ(1.0, Gy_tar[1]);
	EXPECT_EQ(-1.0, Gy_tar[2]);
	EXPECT_EQ(-1.0, Gy_tar[3]);

	EXPECT_EQ(-3.0, h_tar[0]);
	EXPECT_EQ(-6.0, h_tar[1]);
	EXPECT_EQ(-6.0, h_tar[2]);
	EXPECT_EQ(-3.0, h_tar[3]);

	EXPECT_EQ(4.5, prog.o.x);
	EXPECT_EQ(4.5, prog.o.y);
}

void test_linprog2d_condition_problem_offset_and_rescale_single() {
	linprog2d_data_t prog;
	double Gx[1] = {-4.0};
	double Gy[1] = {1.0};
	double h[1] = {8.0};
	double Gx_tar[1], Gy_tar[1], h_tar[1];

	/* Manually setup the linprog2d_data_t structure */
	linprog2d_reset(&prog, 1U);
	prog.Gx = Gx_tar, prog.Gy = Gy_tar, prog.h = h_tar;

	EXPECT_EQ(TRUE, linprog2d_condition_problem(&prog, 0.0, 1.0, Gx, Gy, h));

	EXPECT_EQ(-1.0, Gx_tar[0]);
	EXPECT_EQ(0.25, Gy_tar[0]);
	EXPECT_EQ(2.0, h_tar[0]); /* Only rescaling has happened, no shifting */
}

void test_linprog2d_condition_problem_offset_and_rescale() {
	linprog2d_data_t prog;
	double Gx[2] = {-4.0, -8.0};
	double Gy[2] = {4.0, -8.0};
	double h[2] = {8.0, -24.0};
	double Gx_tar[2], Gy_tar[2], h_tar[2];

	/* Manually setup the linprog2d_data_t structure */
	linprog2d_reset(&prog, 2U);
	prog.Gx = Gx_tar, prog.Gy = Gy_tar, prog.h = h_tar;

	EXPECT_EQ(TRUE, linprog2d_condition_problem(&prog, 0.0, 1.0, Gx, Gy, h));

	EXPECT_EQ(-1.0, Gx_tar[0]);
	EXPECT_EQ(-1.0, Gx_tar[1]);

	EXPECT_EQ(1.0, Gy_tar[0]);
	EXPECT_EQ(-1.0, Gy_tar[1]);

	EXPECT_EQ(0.0, h_tar[0]);
	EXPECT_EQ(0.0, h_tar[1]);

	EXPECT_EQ(0.5, prog.o.x);
	EXPECT_EQ(2.5, prog.o.y);
}

void test_linprog2d_categorize() {
	linprog2d_data_t prog;

	double Gx[7] = {1.0, -1.0, 0.0, 0.0, 0.5, 0.5, -0.25};
	double Gy[7] = {0.0, 0.0, -1.0, 1.0, 0.1, 5.0, -1.0};
	double h[7] = {2.0, -7.0, -8.0, 2.0, 2.0, 15.0, -11.0};
	unsigned int ceil[7], floor[7];

	/* Manually setup the linprog2d_data_t structure */
	linprog2d_reset(&prog, 7U);
	prog.Gx = Gx, prog.Gy = Gy, prog.h = h;
	prog.ceil = ceil, prog.floor = floor;

	/* There are no contradictory constraints in this example */
	EXPECT_EQ(TRUE, linprog2d_categorize_constraints(&prog));

	/* The left and right boundaries are determined by the first two
	   constraints */
	EXPECT_EQ(2.0, prog.x0);
	EXPECT_EQ(7.0, prog.x1);

	/* There are two ceil constraints and three floor constraints. Abort if
	   these conditions are not met */
	ASSERT_EQ(2U, prog.ceil_len);
	ASSERT_EQ(3U, prog.floor_len);

	EXPECT_EQ(2U, ceil[0]);
	EXPECT_EQ(6U, ceil[1]);

	EXPECT_EQ(3U, floor[0]);
	EXPECT_EQ(4U, floor[1]);
	EXPECT_EQ(5U, floor[2]);
}

void test_linprog2d_calculate_intersect() {
#define LP2D_CI linprog2d_calculate_intersect
	double x, y;

	EXPECT_EQ(TRUE, LP2D_CI(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, &x, &y));
	EXPECT_EQ(0.0, x);
	EXPECT_EQ(0.0, y);

	EXPECT_EQ(TRUE, LP2D_CI(1.0, 0.0, 1.0, 0.0, 1.0, 1.0, &x, &y));
	EXPECT_EQ(1.0, x);
	EXPECT_EQ(1.0, y);

	EXPECT_EQ(TRUE, LP2D_CI(1.0, 0.0, -1.0, 0.0, 1.0, 1.0, &x, &y));
	EXPECT_EQ(-1.0, x);
	EXPECT_EQ(1.0, y);

	EXPECT_EQ(TRUE, LP2D_CI(-4.0, 4.0, 8.0, -8.0, -8.0, -24.0, &x, &y));
	EXPECT_EQ(0.5, x);
	EXPECT_EQ(2.5, y);

	EXPECT_EQ(FALSE, LP2D_CI(1.0, 0.0, 0.0, 1.0, 0.0, 0.0, &x, &y));
#undef LP2D_CI
}

void test_linprog2d_calculate_yoffset_form() {
	double Gx[7] = {1.0, 4.0, 1.2, -8.0, 1.5, 9.0, 1.2};
	double Gy[7] = {2.0, 2.0, 3.5, 16.0, -7.8, -3.0, 2.0};
	double h[7] = {4.0, -1.0, 2.0, -8.0, 0.1, 4.0, 1.0};
	double dx[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double y0[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	unsigned int idcs[3] = {1, 3, 5};

	linprog2d_calculate_yoffset_form(idcs, 3, Gx, Gy, h, dx, y0);

	EXPECT_EQ(0.0, dx[0]);
	EXPECT_EQ(-2.0, dx[1]);
	EXPECT_EQ(0.0, dx[2]);
	EXPECT_EQ(0.5, dx[3]);
	EXPECT_EQ(0.0, dx[4]);
	EXPECT_EQ(3.0, dx[5]);
	EXPECT_EQ(0.0, dx[6]);

	EXPECT_EQ(0.0, y0[0]);
	EXPECT_EQ(-0.5, y0[1]);
	EXPECT_EQ(0.0, y0[2]);
	EXPECT_EQ(-0.5, y0[3]);
	EXPECT_EQ(0.0, y0[4]);
	EXPECT_NEAR(-1.333333333333333333333333333, y0[5]);
	EXPECT_EQ(0.0, y0[6]);
}

void test_linprog2d_eliminate_constraint() {
#define LP2D_EC linprog2d_eliminate_constraint
	/* Parallel constraints. Result only depends on the offset h. */

	double h[2] = {0.0, 1.0};
	double dx[2] = {0.0, 0.0};

	EXPECT_EQ(1U, LP2D_EC(h, dx, 0U, 1U, TRUE, TRUE, FALSE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 1U, 0U, TRUE, TRUE, FALSE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 0U, 1U, FALSE, TRUE, FALSE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 1U, 0U, FALSE, TRUE, FALSE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 0U, 1U, TRUE, TRUE, TRUE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 1U, 0U, TRUE, TRUE, TRUE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 0U, 1U, FALSE, TRUE, TRUE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 1U, 0U, FALSE, TRUE, TRUE));

	h[0] = 1.0;
	h[1] = 0.0;

	EXPECT_EQ(0U, LP2D_EC(h, dx, 0U, 1U, TRUE, TRUE, FALSE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 1U, 0U, TRUE, TRUE, FALSE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 0U, 1U, FALSE, TRUE, FALSE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 1U, 0U, FALSE, TRUE, FALSE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 0U, 1U, TRUE, TRUE, TRUE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 1U, 0U, TRUE, TRUE, TRUE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 0U, 1U, FALSE, TRUE, TRUE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 1U, 0U, FALSE, TRUE, TRUE));

	/* Non-parallel constraint. Result depends on the slope, this being a
	   ceiling or floor constraint and whether the optimum is on the left or
	   right of the intersection. */

	dx[0] = 1.0;
	dx[1] = -1.0;

	EXPECT_EQ(1U, LP2D_EC(h, dx, 0U, 1U, FALSE, FALSE, TRUE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 1U, 0U, FALSE, FALSE, TRUE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 0U, 1U, TRUE, FALSE, TRUE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 1U, 0U, TRUE, FALSE, TRUE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 0U, 1U, FALSE, FALSE, FALSE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 1U, 0U, FALSE, FALSE, FALSE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 0U, 1U, TRUE, FALSE, FALSE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 1U, 0U, TRUE, FALSE, FALSE));

	dx[0] = -1.0;
	dx[1] = 1.0;

	EXPECT_EQ(0U, LP2D_EC(h, dx, 0U, 1U, FALSE, FALSE, TRUE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 1U, 0U, FALSE, FALSE, TRUE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 0U, 1U, TRUE, FALSE, TRUE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 1U, 0U, TRUE, FALSE, TRUE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 0U, 1U, FALSE, FALSE, FALSE));
	EXPECT_EQ(1U, LP2D_EC(h, dx, 1U, 0U, FALSE, FALSE, FALSE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 0U, 1U, TRUE, FALSE, FALSE));
	EXPECT_EQ(0U, LP2D_EC(h, dx, 1U, 0U, TRUE, FALSE, FALSE));

#undef LP2D_EC
}

void test_linprog2d_calculate_intersects() {
	linprog2d_data_t prog;
	double Gx[9] = {1.0, -1.0, 0.0, 0.0, 0.5, 0.5, -0.25, 4.0, 2.0};
	double Gy[9] = {0.0, 0.0, -1.0, 1.0, 0.1, 5.0, -1.0, -1.0, 9.0};
	double h[9] = {2.0, -7.0, -8.0, 2.0, 2.0, 15.0, -11.0, 5.0, 8.0};
	double dx[9];
	double y0[9];
	double x_intersect[4];
	unsigned int ceil[9], floor[9], tmp[9];

	/* Manually setup the linprog2d_data_t structure */
	linprog2d_reset(&prog, 9U);
	prog.Gx = Gx, prog.Gy = Gy, prog.h = h, prog.dx = dx, prog.y0 = y0;
	prog.x_intersect = x_intersect, prog.ceil = ceil, prog.floor = floor;
	prog.tmp = tmp;

	/* We don't have parallel constraints in the example, so we can skip
	   conditioning (eliminate_constraints will return incorrect results for
	   parallel constraints without conditioning) */
	linprog2d_categorize_constraints(&prog);
	EXPECT_EQ(3U, prog.ceil_len);
	EXPECT_EQ(4U, prog.floor_len);

	linprog2d_calculate_yoffset_form(ceil, prog.ceil_len, Gx, Gy, h, dx, y0);
	linprog2d_calculate_yoffset_form(floor, prog.floor_len, Gx, Gy, h, dx, y0);

	prog.intersect_len = 0U;
	linprog2d_calculate_intersects(&prog, ceil, &prog.ceil_len, TRUE, FALSE, 0,
	                               FALSE);
	EXPECT_EQ(0U, prog.intersect_len);
	EXPECT_EQ(2U, prog.ceil_len);
	EXPECT_EQ(2U, ceil[0]);
	EXPECT_EQ(7U, ceil[1]);

	linprog2d_calculate_intersects(&prog, floor, &prog.floor_len, FALSE, FALSE,
	                               0, FALSE);
	EXPECT_EQ(1U, prog.intersect_len);
	EXPECT_EQ(3U, prog.floor_len);
	EXPECT_EQ(3U, floor[0]);
	EXPECT_EQ(4U, floor[1]);
	EXPECT_EQ(5U, floor[2]);

	EXPECT_EQ(3.6, prog.x_intersect[0]);
}

void test_linprog2d_track_min_max() {
	double dx[5] = {-1., -2., -8., -4., -8.};
	double y0[5] = {2., 4., 32., 8., 16.};
	unsigned int idcs[4] = {0, 1, 3, 4};
	struct linprog2d_extremum e;

	e = linprog2d_track_extrema(2.0, dx, y0, idcs, 4, TRUE);
	EXPECT_EQ(0.0, e.y);
	EXPECT_EQ(-8.0, e.min_dx);
	EXPECT_EQ(-1.0, e.max_dx);
	EXPECT_EQ(TRUE, e.valid);

	e = linprog2d_track_extrema(2.0, dx, y0, idcs, 4, FALSE);
	EXPECT_EQ(0.0, e.y);
	EXPECT_EQ(-8.0, e.min_dx);
	EXPECT_EQ(-1.0, e.max_dx);
	EXPECT_EQ(TRUE, e.valid);

	e = linprog2d_track_extrema(1.0, dx, y0, idcs, 4, TRUE);
	EXPECT_EQ(1.0, e.y);
	EXPECT_EQ(-1.0, e.min_dx);
	EXPECT_EQ(-1.0, e.max_dx);
	EXPECT_EQ(TRUE, e.valid);

	e = linprog2d_track_extrema(1.0, dx, y0, idcs, 4, FALSE);
	EXPECT_EQ(8.0, e.y);
	EXPECT_EQ(-8.0, e.min_dx);
	EXPECT_EQ(-8.0, e.max_dx);
	EXPECT_EQ(TRUE, e.valid);

	e = linprog2d_track_extrema(3.0, dx, y0, idcs, 4, TRUE);
	EXPECT_EQ(-8.0, e.y);
	EXPECT_EQ(-8.0, e.min_dx);
	EXPECT_EQ(-8.0, e.max_dx);
	EXPECT_EQ(TRUE, e.valid);

	e = linprog2d_track_extrema(3.0, dx, y0, idcs, 4, FALSE);
	EXPECT_EQ(-1.0, e.y);
	EXPECT_EQ(-1.0, e.min_dx);
	EXPECT_EQ(-1.0, e.max_dx);
	EXPECT_EQ(TRUE, e.valid);

	e = linprog2d_track_extrema(3.0, dx, y0, idcs, 0, FALSE);
	EXPECT_EQ(FALSE, e.valid);
}

/* Macro assembling a linprog2d_data instance on the stack */
#define MKPROG(C)                                                         \
	linprog2d_result_t res;                                               \
	linprog2d_data_t prog;                                                \
	double Gx[C], Gy[C], h[C], dx[C], y0[C], x_intersect[C];              \
	unsigned int ceil[C], floor[C], tmp[C];                               \
	prog.Gx = Gx, prog.Gy = Gy, prog.h = h, prog.dx = dx, prog.y0 = y0;   \
	prog.x_intersect = x_intersect, prog.ceil = ceil, prog.floor = floor; \
	prog.capacity = C;                                                    \
	prog.tmp = tmp;

void test_linprog2d_solve_vee() {
	/* Result has its minimum at (0, 0)

	    \xxxx^xxxx/
	     \xxx|xxx/
	      \xx|xx/
	       \x|x/
	        \|/
	  -------X------->
	        /|\
	       / | \
	      /  |  \                       */

	double Gx_src[2] = {1.0, -1.0};
	double Gy_src[2] = {1.0, 1.0};
	double h_src[2] = {0.0, 0.0};

	MKPROG(2U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 2U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_EQ(0.0, res.x1);
	EXPECT_EQ(0.0, res.y1);
}

void test_linprog2d_solve_vee_offset() {
	/* Result has its minimum at (1, 2)

	    \x^xxxxxxx/
	     \|xxxxxx/
	      |xxxxx/
	      |\xxx/
	      | \x/
	      |  X
	      | / \
	------|/-------------->
	      /     \                         */

	double Gx_src[2] = {1.0, -1.0};
	double Gy_src[2] = {1.0, 1.0};
	double h_src[2] = {3.0, 1.0};

	MKPROG(2U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 2U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_EQ(1.0, res.x1);
	EXPECT_EQ(2.0, res.y1);
}

void test_linprog2d_solve_vee_offset_parallel1() {
	/* Result has its minimum at (1, 2)

	 \  \x^xxxxxxx/  /
	  \  \|xxxxxx/  /
	   \  |xxxxx/  /
	    \ |\xxx/  /
	     \| \x/  /
	      \  X  /
	      |\/ \/
	------|/\-/----------->
	      /  /  \                         */

	double Gx_src[4] = {1.0, -1.0, -1.0, 1.0};
	double Gy_src[4] = {1.0, 1.0, 1.0, 1.0};
	double h_src[4] = {3.0, -1.0, 1.0, 0.0};

	MKPROG(4U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 4U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_EQ(1.0, res.x1);
	EXPECT_EQ(2.0, res.y1);
}

void test_linprog2d_solve_vee_offset_parallel2() {
	double Gx_src[4] = {1.0, -1.0, -1.0, 1.0};
	double Gy_src[4] = {1.0, 1.0, 1.0, 1.0};
	double h_src[4] = {3.0, 1.0, -1.0, 0.0};

	MKPROG(4U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 4U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_EQ(1.0, res.x1);
	EXPECT_EQ(2.0, res.y1);
}

void test_linprog2d_solve_vee_offset_parallel3() {
	double Gx_src[4] = {1.0, -1.0, 1.0, -1.0};
	double Gy_src[4] = {1.0, 1.0, 1.0, 1.0};
	double h_src[4] = {3.0, 1.0, 0.0, -1.0};

	MKPROG(4U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 4U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_EQ(1.0, res.x1);
	EXPECT_EQ(2.0, res.y1);
}

void test_linprog2d_solve_vee_offset_parallel4() {
	double Gx_src[4] = {1.0, 1.0, -1.0, -1.0};
	double Gy_src[4] = {1.0, 1.0, 1.0, 1.0};
	double h_src[4] = {3.0, 0.0, 1.0, -1.0};

	MKPROG(4U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 4U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_EQ(1.0, res.x1);
	EXPECT_EQ(2.0, res.y1);
}

void test_linprog2d_solve_vee_offset_rotated() {
	/* Result has its minimum at (1, 2)

	xxxx\ ^       /
	xxxxx\|      /
	xxxxxx|     /
	xxxxxx|\   /
	xxxxxx|x\ /
	xxxxxx|xxX
	xxxxxx|x/ \
	------|/-------------->
	xxxxxx/     \                         */

	double Gx_src[2] = {-1.0, -1.0};
	double Gy_src[2] = {1.0, -1.0};
	double h_src[2] = {1.0, -3.0};

	MKPROG(2U)

	res = linprog2d_solve(&prog, -1.0, 0.0, Gx_src, Gy_src, h_src, 2U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_NEAR(1.0, res.x1);
	EXPECT_NEAR(2.0, res.y1);
}

void test_linprog2d_single_floor_horz_unbounded() {
	/* Result is unbounded.

	xxxxxx^xxxxxxxxxxxxxxx
	xxxxxx|xxxxxxxxxxxxxxx
	xxxxxx|xxxxxxxxxxxxxxx
	xxxxxx|xxxxxxxxxxxxxxx
	X-X-X-X-X-X-X-X-X-X-X-
	      |
	      |
	------|--------------->
	      |                               */

	double Gx_src[1] = {0.0};
	double Gy_src[1] = {1.0};
	double h_src[1] = {1.0};

	MKPROG(1U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 1U);
	EXPECT_EQ(LP2D_UNBOUNDED, res.status);
}

void test_linprog2d_single_floor_horz_edge() {
	/* Result is on a line.

	   |xx^xxxxxxx|
	   |xx|xxxxxxx|
	   |xx|xxxxxxx|
	   |xx|xxxxxxx|
	---|XXXXXXXXXX|-------
	   |  |       |
	   |  |       |
	---|--|-------|------->
	   |  |       |                       */

	double Gx_src[3] = {0.0, 1.0, -1.0};
	double Gy_src[3] = {1.0, 0.0, 0.0};
	double h_src[3] = {1.0, -2.0, -3.0};

	MKPROG(3U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 3U);
	EXPECT_EQ(LP2D_EDGE, res.status);
	EXPECT_EQ(-2.0, res.x1);
	EXPECT_EQ(1.0, res.y1);
	EXPECT_EQ(3.0, res.x2);
	EXPECT_EQ(1.0, res.y2);
}

void test_linprog2d_single_floor_ceil_parallel1() {
	/* Result is on a line.

	      ^
	----------------------
	xxxxxx|xxxxxxxxxxxxxxx
	xxxxxx|xxxxxxxxxxxxxxx
	X-X-X-X-X-X-X-X-X-X-X-
	      |
	      |
	------|--------------->
	      |                               */

	double Gx_src[2] = {0.0, 0.0};
	double Gy_src[2] = {1.0, -1.0};
	double h_src[2] = {1.0, -3.0};

	MKPROG(2U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 2U);
	EXPECT_EQ(LP2D_UNBOUNDED, res.status);
}

void test_linprog2d_single_floor_ceil_parallel2() {
	/* Result is on a line.

	      ^
	----------------------
	xxxxxx|xxxxxxxxxxxxxxx
	xxxxxx|xxxxxxxxxxxxxxx
	X-X-X-X-X-X-X-X-X-X-X-
	      |
	      |
	------|--------------->
	      |                               */

	double Gx_src[2] = {0.0, 0.0};
	double Gy_src[2] = {1.0, -1.0};
	double h_src[2] = {1.0, 3.0};

	MKPROG(2U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 2U);
	EXPECT_EQ(LP2D_INFEASIBLE, res.status);
}

void test_linprog2d_floor_ceil_intersect_edge1() {
	/* Result is on a line.

	\xxxxx^xxxxx\
	 \xxxx|xxxxxx\
	  \xxx|xxxxxxx\
	   \xx|xxxxxxxx\
	----\X|XXXXXXXXX\-----
	     \|          \
	      |           \
	------|\-----------\-->
	      | \           \                 */

	double Gx_src[3] = {0.0, -1.0, 1.0};
	double Gy_src[3] = {1.0, -1.0, 1.0};
	double h_src[3] = {1.0, -5.0, -5.0};

	MKPROG(3U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 3U);
	EXPECT_EQ(LP2D_EDGE, res.status);
	EXPECT_EQ(-6.0, res.x1);
	EXPECT_EQ(1.0, res.y1);
	EXPECT_EQ(4.0, res.x2);
	EXPECT_EQ(1.0, res.y2);
}

void test_linprog2d_floor_ceil_intersect_edge2() {
	/* Result is on a line.
	          /\
	      ^  /xx\
	      | /xxxx\
	      |/xxxxxx\
	      /xxxxxxxx\
	-----/|XXXXXXXXX\-----
	    / |          \
	   /  |           \
	--/---|------------\-->
	 /    |             \                 */

	double Gx_src[3] = {0.0, 1.0, -1.0};
	double Gy_src[3] = {1.0, -1.0, -1.0};
	double h_src[3] = {1.0, -5.0, -5.0};

	MKPROG(3U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 3U);
	EXPECT_EQ(LP2D_EDGE, res.status);
	EXPECT_EQ(-4.0, res.x1);
	EXPECT_EQ(1.0, res.y1);
	EXPECT_EQ(4.0, res.x2);
	EXPECT_EQ(1.0, res.y2);
}

void test_linprog2d_floor_ceil_intersect_edge3() {
	/* Result is on a line.

	      ^  /xxxxxxxxxxx/
	      | /xxxxxxxxxxx/
	      |/xxxxxxxxxxx/
	      /xxxxxxxxxxx/
	---- /X-X-X-X-X-X/----
	    / |         /
	   /  |        /
	--/---|-------/------->
	 /    |      /                        */

	double Gx_src[3] = {0.0, 1.0, -1.0};
	double Gy_src[3] = {1.0, -1.0, 1.0};
	double h_src[3] = {1.0, -5.0, -5.0};

	MKPROG(3U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 3U);
	EXPECT_EQ(LP2D_EDGE, res.status);
	EXPECT_EQ(-4.0, res.x1);
	EXPECT_EQ(1.0, res.y1);
	EXPECT_EQ(6.0, res.x2);
	EXPECT_EQ(1.0, res.y2);
}

void test_linprog2d_floor_floor_intersect_edge() {
	/* Result is on a line.

	\xxxxx^xxxxxxxxxxx/
	 \xxxx|xxxxxxxxxx/
	  \xxx|xxxxxxxxx/
	   \xx|xxxxxxxx/
	----\X|XXXXXXX/-------
	     \|      /
	      |     /
	------|\---/---------->
	      | \ /                           */

	double Gx_src[3] = {0.0, 1.0, -1.0};
	double Gy_src[3] = {1.0, 1.0, 1.0};
	double h_src[3] = {1.0, -5.0, 0.0};

	MKPROG(3U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 3U);
	EXPECT_EQ(LP2D_EDGE, res.status);
	EXPECT_EQ(-6.0, res.x1);
	EXPECT_EQ(1.0, res.y1);
	EXPECT_EQ(1.0, res.x2);
	EXPECT_EQ(1.0, res.y2);
}

void test_linprog2d_vert_infeasible() {
	/* Result is infeasible.

	<|    ^    |>
	----------------------
	<|    |    |>
	<|^ ^ |^ ^ |> ^ ^ ^ ^
	-|----|----|----------
	<|    |    |>
	<|    |    |>
	-|----|----|---------->
	<|    |    |>                         */

	double Gx_src[4] = {0.0, 0.0, 1.0, -1.0};
	double Gy_src[4] = {1.0, -1.0, 0.0, 0.0};
	double h_src[4] = {1.0, -3.0, 5.0, 5.0};

	MKPROG(4U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 4U);
	EXPECT_EQ(LP2D_INFEASIBLE, res.status);
}

void test_linprog2d_hatches() {
	double Gx_src[16] = {  1,  -1,   1,  -1,   1,  -1,   1,  -1,   1,  -1,   1,  -1,   1,  -1,   1,  -1};
	double Gy_src[16] = {  1,   1,   1,   1,   1,   1,   1,   1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1};
	double h_src[16] =  {-20, -20, -15, -15, -10, -10,  -5,  -5, -20, -20, -15, -15, -10, -10,  -5,  -5};

	MKPROG(16U)

	res = linprog2d_solve(&prog, 0.0, 1.0, Gx_src, Gy_src, h_src, 16U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_NEAR(0.0, res.x1);
	EXPECT_NEAR(-5.0, res.y1);
}

void test_linprog2d_nr_example() {
	/* Example from Numerical Recipes 3rd ed. pp. 529; see p. 534 for fig. */

	double Gx_src[3] = {-2.0, 1.0, -1.0};
	double Gy_src[3] = {-1.0, 1.0, -3.0};
	double h_src[3] = {-70.0, 40.0, -90.0};

	MKPROG(3U)

	res = linprog2d_solve(&prog, -40.0, -60.0, Gx_src, Gy_src, h_src, 3U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_NEAR(24.0, res.x1);
	EXPECT_NEAR(22.0, res.y1);
}

void test_linprog2d_barnfm10e_example() {
	/* Example from random lecture found on the internet
	   http://www3.govst.edu/kriordan/files/ssc/math161/ppt/barnfm10e_ppt_5_2.ppt
	 */

	double Gx_src[5] = {1.0, 0.0, -1.0, -8.0, -4.0};
	double Gy_src[5] = {0.0, 1.0, 0.0, -8.0, -12.0};
	double h_src[5] = {0.0, 0.0, -15.0, -160.0, -180.0};

	MKPROG(5U)

	res = linprog2d_solve(&prog, -5.0, -10.0, Gx_src, Gy_src, h_src, 5U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_NEAR(7.5, res.x1);
	EXPECT_NEAR(12.5, res.y1);
}

void test_linprog2d_solve_simple_nr_example() {
	/* Example from Numerical Recipes 3rd ed. pp. 529; see p. 534 for fig. */

	const double Gx_src[3] = {-2.0, 1.0, -1.0};
	const double Gy_src[3] = {-1.0, 1.0, -3.0};
	const double h_src[3] = {-70.0, 40.0, -90.0};

	linprog2d_result_t res =
	    linprog2d_solve_simple(-40.0, -60.0, Gx_src, Gy_src, h_src, 3U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_NEAR(24.0, res.x1);
	EXPECT_NEAR(22.0, res.y1);
}

void test_linprog2d_solve_simple_barnfm10e_example() {
	/* Example from random lecture found on the internet
	   http://www3.govst.edu/kriordan/files/ssc/math161/ppt/barnfm10e_ppt_5_2.ppt
	 */

	const double Gx_src[5] = {1.0, 0.0, -1.0, -8.0, -4.0};
	const double Gy_src[5] = {0.0, 1.0, 0.0, -8.0, -12.0};
	const double h_src[5] = {0.0, 0.0, -15.0, -160.0, -180.0};

	linprog2d_result_t res =
	    linprog2d_solve_simple(-5.0, -10.0, Gx_src, Gy_src, h_src, 5U);
	EXPECT_EQ(LP2D_POINT, res.status);
	EXPECT_NEAR(7.5, res.x1);
	EXPECT_NEAR(12.5, res.y1);
}

void test_linprog2d_solve_simple_fail() {
	/* Try to allocate 240 GiB of memory. Congrats if this does not fail on your
	   machine. */
	linprog2d_result_t res =
	    linprog2d_solve_simple(0.0, 1.0, NULL, NULL, NULL, 0xFFFFFFFFUL);
	EXPECT_EQ(LP2D_ERROR, res.status);
}

/******************************************************************************
 * Main program                                                               *
 ******************************************************************************/

int main() {
	RUN(test_feq);
	RUN(test_memalign64);
	RUN(test_sort);
	RUN(test_partition);
	RUN(test_kth_smallest);
	RUN(test_median);
	RUN(test_linprog2d_normalization_coeff);
#ifndef LINPROG2D_NO_ALLOC
	RUN(test_linprog2d_create_and_capacity);
#endif
	RUN(test_linprog2d_problem_too_large);
	RUN(test_linprog2d_condition_problem_rotation);
	RUN(test_linprog2d_condition_problem_eliminate_invalid);
	RUN(test_linprog2d_condition_problem_offset1);
	RUN(test_linprog2d_condition_problem_offset2);
	RUN(test_linprog2d_condition_problem_offset_and_rescale_single);
	RUN(test_linprog2d_condition_problem_offset_and_rescale);
	RUN(test_linprog2d_categorize);
	RUN(test_linprog2d_calculate_intersect);
	RUN(test_linprog2d_calculate_yoffset_form);
	RUN(test_linprog2d_eliminate_constraint);
	RUN(test_linprog2d_calculate_intersects);
	RUN(test_linprog2d_track_min_max);
	RUN(test_linprog2d_solve_vee);
	RUN(test_linprog2d_solve_vee_offset);
	RUN(test_linprog2d_solve_vee_offset_parallel1);
	RUN(test_linprog2d_solve_vee_offset_parallel2);
	RUN(test_linprog2d_solve_vee_offset_parallel3);
	RUN(test_linprog2d_solve_vee_offset_parallel4);
	RUN(test_linprog2d_solve_vee_offset_rotated);
	RUN(test_linprog2d_single_floor_horz_unbounded);
	RUN(test_linprog2d_single_floor_horz_edge);
	RUN(test_linprog2d_single_floor_ceil_parallel1);
	RUN(test_linprog2d_single_floor_ceil_parallel2);
	RUN(test_linprog2d_floor_ceil_intersect_edge1);
	RUN(test_linprog2d_floor_ceil_intersect_edge2);
	RUN(test_linprog2d_floor_ceil_intersect_edge3);
	RUN(test_linprog2d_floor_floor_intersect_edge);
	RUN(test_linprog2d_vert_infeasible);
	RUN(test_linprog2d_hatches);
	RUN(test_linprog2d_nr_example);
	RUN(test_linprog2d_barnfm10e_example);
#ifndef LINPROG2D_NO_ALLOC
	RUN(test_linprog2d_solve_simple_nr_example);
	RUN(test_linprog2d_solve_simple_barnfm10e_example);
	RUN(test_linprog2d_solve_simple_fail);
#endif

	fprintf(stderr, ANSI_GRAY "=====" ANSI_RESET "\n");
	if (n_failed) {
		fprintf(stderr, ANSI_RED "[ERR]" ANSI_RESET);
	} else {
		fprintf(stderr, ANSI_GREEN "[OK!]" ANSI_RESET);
	}
	fprintf(stderr, " Successful tests: %d; Failed tests: %d\n", n_success,
	        n_failed);
	return n_failed ? 1 : 0;
}

