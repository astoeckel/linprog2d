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

.PHONY: all test clean cov wasm dist

CCFLAGS := -g -L build/ -I . -fPIC --std=c89 -Wall -Wextra -pedantic-errors

all: build/liblinprog2d.a build/liblinprog2d.so \
     build/example/linprog2d_simple \
     build/test/test_linprog2d

build/linprog2d.o: linprog2d.c linprog2d.h
	mkdir -p build
	$(CC) $(CCFLAGS) -c linprog2d.c -o build/linprog2d.o

build/liblinprog2d.a: build/linprog2d.o
	ar rcs build/liblinprog2d.a build/linprog2d.o

build/liblinprog2d.so: build/linprog2d.o
	gcc -shared -o build/liblinprog2d.so build/linprog2d.o -lm

build/linprog2d.wasm: linprog2d.c linprog2d.h
	emcc -DLINPROG2D_REDUCED_INTERFACE -Oz -s WASM=1 -s SIDE_MODULE=1 linprog2d.c -o build/linprog2d.wasm

build/linprog2d.wasm.b64: build/linprog2d.wasm
	base64 -w0 build/linprog2d.wasm > build/linprog2d.wasm.b64

build/linprog2d.js: build/linprog2d.wasm.b64 linprog2d.in.js
	perl -pe 's/<_WASM_CODE_HERE_>/`cat build\/linprog2d.wasm.b64`/ge' linprog2d.in.js > build/linprog2d.js

build/linprog2d.min.js: build/linprog2d.js
	minify build/linprog2d.js > build/linprog2d.min.js # npm i babel-minify

build/example/linprog2d_simple: build/liblinprog2d.a examples/linprog2d_simple.c
	mkdir -p build/examples
	$(CC) $(CCFLAGS) -static -o build/examples/linprog2d_simple examples/linprog2d_simple.c -llinprog2d -lm

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

wasm: build/linprog2d.js build/linprog2d.min.js

dist: wasm
	cp build/linprog2d.min.js dist/

clean:
	rm -Rf \
		*.gcda *.gcno *.gcov *.vgcore \
		build/linprog2d.o \
		build/liblinprog2d.a \
		build/liblinprog2d.so \
		build/linprog2d.js \
		build/linprog2d.min.js \
		build/linprog2d.wasm.b64 \
		build/linprog2d.wasm \
		build/test/test_linprog2d \
		build/test/test_linprog2d_cov \
		test_linprog2d_coverage*.html

