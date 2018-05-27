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
 * @file linprog2d.js
 *
 * JavaScript implementation of 2D linear programming. Uses the C library in
 * this repository compiled to WebAssembly.
 *
 * @author Andreas Stöckel
 */

this.linprog2d = (function(global) {
	'use strict';

	/* Base64 encoded WASM code gets injected here. Run "make wasm" to build the
	   final JavaScript file. */
	const WASM_CODE = '<_WASM_CODE_HERE_>';

	/**
	 * Decodes a base64 string into a Uint8Array.
	 */
	function _decode_base64(s) {
		/* Determine the length of the encoded string */
		let len = Math.ceil(s.length / 4 * 3) | 0;
		if (s.charAt(s.length - 1) == '=') { len--; };
		if (s.charAt(s.length - 2) == '=') { len--; };

		/* Instantiate the byte array containing the data */
		const data = new Uint8Array(new ArrayBuffer(len));

		/* Build a reverse lookup table from characters to alphabet index */
		const A = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
		const M = {'=': 0}; for (let i = 0; i < A.length; i++) { M[A.charAt(i)] = i; }

		/* Iterate over the string and convert quadtruples to triples */
		let j = 0;
		for (let i = 0; i < s.length;) {
			/* Assemble a 32-bit integer containing the decoded data */
			const c1 = M[s.charAt(i++)], c2 = M[s.charAt(i++)];
			const c3 = M[s.charAt(i++)], c4 = M[s.charAt(i++)];
			const d = (c1 << 18) | (c2 << 12) | (c3 << 6) | c4;

			/* Append the three bytes from the 32-bit integer to the data */
			if (j < len) { data[j++] = (d >> 16) & 0xFF; }
			if (j < len) { data[j++] = (d >> 8) & 0xFF; }
			if (j < len) { data[j++] = d & 0xFF; }
		}
		return data;
	}

	let _module; /* WASM module instance. */
	let _memory; /* Global memory space used in the module. */
	let _init = null;

	/**
	 * Loads and initialises the WASM module. Returns a promise which, for
	 * convenience, provides a reference at the solve() function.
	 */
	function init() {
		if (!_init) {
			/* Instantiate the global memory object. */
			_memory = new WebAssembly.Memory({ initial: 256, maximum: 256 });

			/* Compile and instantiate the WASM code. */
			_init = WebAssembly.instantiate(_decode_base64(WASM_CODE), {
				'env': {
					'table': new WebAssembly.Table({
						'initial': 8,
						'element': 'anyfunc'
					}),
					'tableBase': 0,
					'memory': _memory,
					'memoryBase': 0,
					'abort': (i) => { throw null; }
				},
				'global': {
					'Infinity': Infinity
				}
			}).then(module => {
				/* Store the module in the global variable */
				_module = module.instance.exports;

				/* Run post-initialisation code */
				_module.__post_instantiate();

				/* Return the "solve" function */
				return solve;
			});
		}
		return _init;
	}

	/**
	 * Solves a 2D linear programming problem of the form
	 *
	 * minimize c.x * x + c.y * y
	 * w.r.t.   Gx[i] * x + Gy[i] * y >= h[i] for all i
	 *
	 * where cx, cy are doubles and Gx, Gy, h are double arrays of the same
	 * length. This function must only be called after init() has completed.
	 * This function returns an object containing the following data
	 * {
	 *     'status': <one of the LP2D_* constants defined below>,
	 *     'x1': <solution point x-coordinate>,
	 *     'y1': <solution point y-coordinate>,
	 *     'x2': <solution second edge point x-coordinate>,
	 *     'y2': <solution second edge point y-coordinate>
	 * }
	 */
	function solve(cx, cy, Gx, Gy, h) {
		/* Make sure the input is valid */
		if ((Gx.length != Gy.length) || (Gy.length != h.length)) {
			throw 'Invalid input';
		}

		/* Compute all target pointers */
		const n = Gx.length, offs = 0x80000; /* skip the first page */
		const Gx_ptr = offs, Gy_ptr = Gx_ptr + n * 8;
		const h_ptr = Gy_ptr + n * 8, res_ptr = h_ptr + n * 8;
		const base_ptr = res_ptr + 5 * 8;

		/* Copy the input to the WebAssembly memory */
		const mu32 = new Uint32Array(_memory.buffer);
		const mf64 = new Float64Array(_memory.buffer);
		for (let i = 0; i < n; i++) {
			mf64[Gx_ptr / 8 + i] = Gx[i];
			mf64[Gy_ptr / 8 + i] = Gy[i];
			mf64[h_ptr / 8 + i] = h[i];
		}

		/* Instantiate the linprog2d instance */
		const prog = _module._linprog2d_init(n, base_ptr);

		/* Solve the problem */
		_module._linprog2d_solve(res_ptr, prog, cx, cy, Gx_ptr, Gy_ptr, h_ptr, n);

		/* Read the result */
		return {
			'x1': mf64[res_ptr / 8 + 0],
			'y1': mf64[res_ptr / 8 + 1],
			'x2': mf64[res_ptr / 8 + 2],
			'y2': mf64[res_ptr / 8 + 3],
			'status': mu32[res_ptr / 4 + 8]
		}
	}

	return {
		'init': init,
		'solve': solve,
		'ERROR': 0,
		'INFEASIBLE': 1,
		'UNBOUNDED': 2,
		'EDGE': 3,
		'POINT': 4
	};
})(this);
