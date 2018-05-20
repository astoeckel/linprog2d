#  linprog2d --- Two-dimensional linear programming solver
#  Copyright (C) 2018 Andreas Stöckel
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

.PHONY: all test clean

CCFLAGS := --std=c89 -Wall -Wextra -pedantic-errors

all: linprog2d.o test_linprog2d test_linprog2d_c

test_linprog2d: linprog2d.c linprog2d.h
	$(CC) $(CCFLAGS) -DLINPROG_2D_TEST -o test_linprog2d linprog2d.c -lm

test: test_linprog2d
	./test_linprog2d

clean:
	rm -Rf *.o test_linprog2d

