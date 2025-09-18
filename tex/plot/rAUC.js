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
	let max_fpr = 0.1, min_sen = 0.5;
	for (let i = 0; i < args.length; ++i) {
		const fn = args[i];
		let a = [];
		for (const line of k8_readline(fn)) {
			let t = line.split("\t");
			if (t[0] != "BN") continue;
			a.push({ x:parseFloat(t[9]), y:parseFloat(t[8])});
		}
		let x0 = 0, y0 = 0;
		let area = 0.0;
		for (let j = 0; j < a.length; ++j) {
			const x1 = a[j].x, y1 = a[j].y;
			let da = 0.0;
			if (y1 >= min_sen && x0 <= max_fpr) {
				const dx = x1 - x0;
				const dy = y1 - y0;
				if (y0 >= min_sen && x1 <= max_fpr) {
					da = dx * ((y0 - min_sen) + (y1 - min_sen)) / 2;
				} else if (y0 < min_sen) {
					const dx1 = dx * ((y1 - min_sen) / dy);
					da = dx1 * (y1 - min_sen) / 2;
				} else if (x1 > max_fpr) {
					const dy1 = dy * ((max_fpr - x0) / dx);
					da = dy1 * (max_fpr - x0) / 2 + (max_fpr - x0) * (y0 - min_sen);
				}
			}
			area += da;
			//print(da, x1, y1, x0, y0);
			x0 = x1, y0 = y1;
		}
		print((area / ((1 - min_sen) * max_fpr)).toFixed(3), fn);
	}
}

main(arguments);
