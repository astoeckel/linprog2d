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

"""
Python wrapper for linprog2d and numpy. Loads liblinprog2d.so using the ctypes
foreign function interface.
"""

import ctypes
import numpy as np

linprog2d_solve_simple = None

def init():
    global linprog2d_solve_simple

    # Abort if the library has already been loaded
    if not linprog2d_solve_simple is None:
        return

    # Import linprog2d_solve_simple using the foreign function interface
    from ctypes import cdll, c_double, c_uint32, POINTER, Structure, pointer
    liblinprog2d = cdll.LoadLibrary("liblinprog2d.so")
    c_linprog2d_solve_simple = liblinprog2d.linprog2d_solve_simple

    # Define datatypes
    class linprog2d_result(Structure):
        _fields_ = [
            ('x1', c_double),
            ('y1', c_double),
            ('x2', c_double),
            ('y2', c_double),
            ('status', c_uint32)
        ]
    c_double_p = POINTER(c_double)
    c_linprog2d_result_p = POINTER(linprog2d_result)

    # Define a wrapper function
    def linprog2d_solve_simple_(cx, cy, Gx, Gy, h):
        result = linprog2d_result(1.0, 2.0, 3.0, 4.0, 5)
        c_linprog2d_solve_simple(
            c_linprog2d_result_p(result),
            c_double(cx),
            c_double(cy),
            Gx.ctypes.data_as(c_double_p),
            Gy.ctypes.data_as(c_double_p),
            h.ctypes.data_as(c_double_p), c_uint32(h.size))
        return result

    # Store the wrapper in the global variable
    linprog2d_solve_simple = linprog2d_solve_simple_

def solve(cx, cy, Gx, Gy, h):
    """
    Solves the given two-dimensional linear programming problem.

    Parameters
    ==========

    cx: x-component of the gradient vector
    cy: y-component of the gradient vector
    Gx: is an n-element vector containing the x-components of the constraint
        normal vectors, where n is the number of constraints.
    Gy: is an n-element vector containing the y-components of the constraint
        normal vectors, where n is the number of constraints.
    h: is an n-element vector containing the offsets of the constraint vectors.
    """

    # Load the library
    init()

    # Convert the incoming parameters to numpy arrays of the correct layout
    Gx = np.atleast_1d(Gx).astype(np.float64, order='C', copy=False)
    Gy = np.atleast_1d(Gy).astype(np.float64, order='C', copy=False)
    h = np.atleast_1d(h).astype(np.float64, order='C', copy=False)

    # Make sure the arrays have the correct layout
    assert (Gx.size == Gy.size == h.size)

    # Call linprog2d_solve_simple
    return linprog2d_solve_simple(cx, cy, Gx, Gy, h)

# Status codes
STATUS_ERROR = 0
STATUS_INFEASIBLE = 1
STATUS_UNBOUNDED = 2
STATUS_EDGE = 3
STATUS_POINT = 4

