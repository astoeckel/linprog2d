# linprog2d -- Linear Programming Solver for Two Dimensional Problems

`liblinprog2d` is a small library for solving two-dimensional linear programming problems in linear time (with respect to the number of constraints). It comes with bindings for Python and a JavaScript/WebAssembly version.

The library is written in pure C89/C90 and has very little dependencies aside from `sqrt` and `malloc`, `free` (the latter can be deactivated).

## How to use



## How to build

`liblinprog2d` is written in the C++-compatible subset of C89/C90. So all you need is a standards-compliant C or C++ compiler. Since `liblinprog2d` is just a single source file, building should be fairly straight forward on any platform; just compile `linprog2d.c` or `test/test_linprog2d.c` if you want to execute the unit tests.

For your convenience, this repository comes with a Makefile that should work on most Unix-like platforms. To build and test `liblinprog2d` just execute the following on your command line:
```sh
git clone https://github.com/astoeckel/linprog2d
cd linprog2d
make && make test
```

Building the JavaScript/WASM version is a little bit more complicated. You’ll need to install `emsdk`, as well as the `npm` package `babel-minify`. Make sure that the `emsdk` environment is active, then execute the following in the `linprog2d` directory:
```sh
make wasm
```
To run the unit-tests in the JavaScript/WebAssembly version you'll have to install `nodejs`. Execute the following command in the `linprog2d` directory
```sh
emcc test/test_linprog2d.c && node ./a.out.js
```

## License

```
linprog2d --- Two-dimensional linear programming solver
Copyright (C) 2018 Andreas Stöckel

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
```

