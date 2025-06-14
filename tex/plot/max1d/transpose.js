#!/usr/bin/env k8

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
	let a = [];
	for (const line of k8_readline(args[0]))
		a.push(line.split("\t"));
	for (let i = 0; i < a[0].length; ++i) {
		let t = [];
		for (let j = 0; j < a.length; ++j)
			t.push(a[j][i]);
		print(t.join("\t"));
	}
}

main(arguments);
