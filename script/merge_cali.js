#!/usr/bin/env k8

"use strict";

function* k8_readline(fn) {
	let buf = new Bytes();
	let file = new File(fn);
	while (file.readline(buf) >= 0) {
		yield buf.toString();
	}
	file.close();
	buf.destroy();
}

function main(args) {
	if (args.length < 3) {
		print("Usage: merge_cali.js <r1>,<r2>,... <f1.eval> <f2.eval> ...");
		return;
	}
	let frac = args[0].split(",");
	if (frac.length != args.length - 1)
		throw Error("number of fractions don't match the number of files");
	let step = -1.0, n_bin = -1, a = [];
	for (let i = 1; i < args.length; ++i) {
		let f = parseFloat(frac[i-1]);
		if (f <= 0.0) continue;
		if (f > 1.0) f = 1.0;
		for (const line of k8_readline(args[i])) {
			let t = line.split("\t");
			if (t[0] == "ST") {
				const s = parseFloat(t[1]);
				if (step <= 0.0) step = s;
				else if (s < step - 1e-6 || s > step + 1e-6)
					throw("different steps");
			} else if (t[0] == "NB") {
				const n = parseInt(t[1]);
				if (n_bin <= 0) {
					n_bin = n;
					for (let j = 0; j < n_bin; ++j)
						a[j] = { tot:0, pos:0 };
				} else if (n_bin != n) {
					throw("different number of bins");
				}
			} else if (t[0] == "BN") {
				const k = parseInt(t[1]);
				const tot = parseInt(t[2]);
				const pos = parseInt(t[3]);
				a[k].tot += Math.floor(tot * f + .499);
				a[k].pos += Math.floor(pos * f + .499);
			}
		}
	}
	print("ST", step);
	print("NB", n_bin);
	for (let j = n_bin - 1; j >= 0; --j)
		print("BN", j, a[j].tot, a[j].pos);
}

main(arguments);
