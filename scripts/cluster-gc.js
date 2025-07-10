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

function two_mean(a) {
	a.sort(function(p,q) { return p.y / p.x - q.y / q.x });
	let X = 0, Y = 0, X2 = 0, Y2 = 0, XY = 0;
	for (let i = 0; i < a.length; ++i)
		X += a[i].x, Y += a[i].y, X2 += a[i].x * a[i].x, Y2 += a[i].y * a[i].y, XY += a[i].x * a[i].y;
	const avg = Y / X;
	let min_i = -1, min_d, min_u, min_v, min_x, min_y;
	let x = 0, y = 0, x2 = 0, y2 = 0, xy = 0; // y2 is not necessary but might help precision
	for (let i = 0; i < a.length - 1; ++i) {
		if (i >= 1) {
			let u = y / x, v = (Y - y) / (X - x);
			const d = (y2 + u * u * x2 - 2 * u * xy) + ((Y2 - y2) + v * v * (X2 - x2) - 2 * v * (XY - xy));
			if (min_i < 0 || d < min_d)
				min_i = i, min_d = d, min_u = u, min_v = v, min_x = x, min_y = y;
		}
		x += a[i].x, y += a[i].y, x2 += a[i].x * a[i].x, y2 += a[i].y * a[i].y, xy += a[i].x * a[i].y;
	}
	let d1 = 0, d2 = 0;
	for (let i = 0; i < a.length; ++i)
		d1 += (a[i].y - a[i].x * avg) * (a[i].y - a[i].x * avg);
	for (let i = 0; i < min_i; ++i)
		d2 += (a[i].y - a[i].x * min_u) * (a[i].y - a[i].x * min_u);
	for (let i = min_i; i < a.length; ++i)
		d2 += (a[i].y - a[i].x * min_v) * (a[i].y - a[i].x * min_v);
	d1 = Math.sqrt(d1 / X2);
	d2 = Math.sqrt(d2 / X2);
	return [avg, min_i, min_u, min_v, d1, d2];
}

function main(args) {
	if (args.length < 1) {
		print("Usage: cluster-gc.js <minigff-intron.out>");
		return;
	}
	let a = [];
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		const gc = parseInt(t[7]), at = parseInt(t[6]);
		a.push({ x: gc+at, y: gc, line: line });
	}
	const [avg, k, u1, u2, d1, d2] = two_mean(a);
	print("UA", a.length, avg.toFixed(6));
	print("U1", k, u1.toFixed(6));
	print("U2", a.length - k, u2.toFixed(6));
	print("DF", d2.toFixed(6), d1.toFixed(6));
}

main(arguments);
