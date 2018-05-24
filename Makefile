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

.PHONY: all test clean cov

CCFLAGS := -g -fPIC --std=c89 -Wall -Wextra -pedantic-errors

all: build/liblinprog2d.a build/liblinprog2d.so build/test/test_linprog2d

build/linprog2d.o: linprog2d.c linprog2d.h
	mkdir -p build
	$(CC) $(CCFLAGS) -c linprog2d.c -o build/linprog2d.o

build/liblinprog2d.a: build/linprog2d.o
	ar rcs build/liblinprog2d.a build/linprog2d.o

build/liblinprog2d.so: build/linprog2d.o
	gcc -shared -o build/liblinprog2d.so build/linprog2d.o -lm

build/test/test_linprog2d: test/test_linprog2d.c linprog2d.c linprog2d.h
	mkdir -p build/test
	$(CC) $(CCFLAGS) -o build/test/test_linprog2d test/test_linprog2d.c -lm

build/test/test_linprog2d_cov: test/test_linprog2d.c linprog2d.c linprog2d.h
	mkdir -p build/test
	$(CC) $(CCFLAGS) -O0 -fprofile-arcs -ftest-coverage -o build/test/test_linprog2d_cov test/test_linprog2d.c -lm

test: build/test/test_linprog2d
	./build/test/test_linprog2d

cov: build/test/test_linprog2d_cov
	./build/test/test_linprog2d_cov
	gcovr -e test/test_linprog2d.c -r . --html --html-details -o test_linprog2d_coverage.html

clean:
	rm -Rf *.gcda *.gcno *.gcov *.vgcore build/linprog2d.o build/liblinprog2d.a build/liblinprog2d.so build/test/test_linprog2d build/test/test_linprog2d_cov test_linprog2d_coverage*.html

