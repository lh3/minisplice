#!/usr/bin/env k8

"use strict";

Array.prototype.delete_at = function(i) {
	for (let j = i; j < this.length - 1; ++j)
		this[j] = this[j + 1];
	--this.length;
}

function* getopt(argv, ostr, longopts) {
	if (argv.length == 0) return;
	let pos = 0, cur = 0;
	while (cur < argv.length) {
		let lopt = "", opt = "?", arg = "";
		while (cur < argv.length) { // skip non-option arguments
			if (argv[cur][0] == "-" && argv[cur].length > 1) {
				if (argv[cur] == "--") cur = argv.length;
				break;
			} else ++cur;
		}
		if (cur == argv.length) break;
		let a = argv[cur];
		if (a[0] == "-" && a[1] == "-") { // a long option
			pos = -1;
			let c = 0, k = -1, tmp = "", o;
			const pos_eq = a.indexOf("=");
			if (pos_eq > 0) {
				o = a.substring(2, pos_eq);
				arg = a.substring(pos_eq + 1);
			} else o = a.substring(2);
			for (let i = 0; i < longopts.length; ++i) {
				let y = longopts[i];
				if (y[y.length - 1] == "=") y = y.substring(0, y.length - 1);
				if (o.length <= y.length && o == y.substring(0, o.length)) {
					k = i, tmp = y;
					++c; // c is the number of matches
					if (o == y) { // exact match
						c = 1;
						break;
					}
				}
			}
			if (c == 1) { // find a unique match
				lopt = tmp;
				if (pos_eq < 0 && longopts[k][longopts[k].length-1] == "=" && cur + 1 < argv.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				}
			}
		} else { // a short option
			if (pos == 0) pos = 1;
			opt = a[pos++];
			let k = ostr.indexOf(opt);
			if (k < 0) {
				opt = "?";
			} else if (k + 1 < ostr.length && ostr[k+1] == ":") { // requiring an argument
				if (pos >= a.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				} else arg = a.substring(pos);
				pos = -1;
			}
		}
		if (pos < 0 || pos >= argv[cur].length) {
			argv.delete_at(cur);
			pos = 0;
		}
		if (lopt != "") yield { opt: `--${lopt}`, arg: arg };
		else if (opt != "?") yield { opt: `-${opt}`, arg: arg };
		else yield { opt: "?", arg: "" };
	}
}

function* k8_readline(fn) {
	let buf = new Bytes();
	let file = new File(fn);
	while (file.readline(buf) >= 0) {
		yield buf.toString();
	}
	file.close();
	buf.destroy();
}

function gff2bed(args)
{
	let fn_ucsc_fai = null, is_short = false, keep_gff = false, print_junc = false, output_gene = false, ens_canon_only = false;
	for (const o of getopt(args, "u:sgjGe")) {
		if (o.opt == 'u') fn_ucsc_fai = o.arg;
		else if (o.opt == 's') is_short = true;
		else if (o.opt == 'g') keep_gff = true;
		else if (o.opt == 'j') print_junc = true;
		else if (o.opt == 'G') output_gene = true;
		else if (o.opt == 'e') ens_canon_only = true;
	}

	if (args.length == 0) {
		print("Usage: gff2bed.js [options] <in.gff>");
		print("Options:");
		print("  -j       output junction BED");
		print("  -s       print names in the short form");
		print("  -u FILE  hg38.fa.fai for chr name conversion");
		print("  -e       only show transcript tagged with 'Ensembl_canonical'");
		print("  -g       output GFF (used with -u)");
		return;
	}

	let ens2ucsc = {};
	if (fn_ucsc_fai != null) {
		for (const line of k8_readline(fn_ucsc_fai)) {
			let t = line.split("\t");
			var s = t[0];
			if (/_(random|alt|decoy)$/.test(s)) {
				s = s.replace(/_(random|alt|decoy)$/, '');
				s = s.replace(/^chr\S+_/, '');
			} else {
				s = s.replace(/^chrUn_/, '');
			}
			s = s.replace(/v(\d+)/, ".$1");
			if (s != t[0]) ens2ucsc[s] = t[0];
		}
	}

	const colors = {
		'protein_coding':'0,128,255',
		'mRNA':'0,128,255',
		'lincRNA':'0,192,0',
		'snRNA':'0,192,0',
		'miRNA':'0,192,0',
		'misc_RNA':'0,192,0'
	};

	function print_bed12(exons, cds_st, cds_en, is_short, print_junc)
	{
		if (exons.length == 0) return;
		var name = is_short? exons[0][7] + "|" + exons[0][5] : exons[0].slice(4, 7).join("|");
		var a = exons.sort(function(a,b) {return a[1]-b[1]});
		if (print_junc) {
		for (var i = 1; i < a.length; ++i)
			print(a[i][0], a[i-1][2], a[i][1], name, 1000, a[i][3]);
			return;
		}
		var sizes = [], starts = [], st, en;
		st = a[0][1];
		en = a[a.length - 1][2];
		if (cds_st == 1<<30) cds_st = st;
		if (cds_en == 0) cds_en = en;
		if (cds_st < st || cds_en > en)
			throw Error("inconsistent thick start or end for transcript " + a[0][4]);
		for (var i = 0; i < a.length; ++i) {
			sizes.push(a[i][2] - a[i][1]);
			starts.push(a[i][1] - st);
		}
		let color = colors[a[0][5]];
		if (color == null) color = '196,196,196';
		print(a[0][0], st, en, name, 1000, a[0][3], cds_st, cds_en, color, a.length, sizes.join(",") + ",", starts.join(",") + ",");
	}

	const re_gtf  = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|transcript_name|tag) "([^"]+)";/g;
	const re_gff3 = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|transcript_name)=([^;]+)/g;
	const re_gtf_gene  = /\b(gene_id|gene_type|gene_name) "([^;]+)";/g;
	const re_gff3_gene = /\b(gene_id|gene_type|source_gene|gene_biotype|gene_name)=([^;]+);/g;

	let exons = [], cds_st = 1<<30, cds_en = 0, last_id = null;
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		if (keep_gff) {
			if (t[0].charAt(0) != '#' && ens2ucsc[t[0]] != null)
				t[0] = ens2ucsc[t[0]];
			print(t.join("\t"));
			continue;
		}
		if (t[0].charAt(0) == '#') continue;
		if (output_gene) {
			let id = null, src = null, biotype = null, type = "", name = "N/A";
			if (t[2] != "gene") continue;
			while ((m = re_gtf_gene.exec(t[8])) != null) {
				if (m[1] == "gene_id") id = m[2];
				else if (m[1] == "gene_type") type = m[2];
				else if (m[1] == "gene_name") name = m[2];
			}
			while ((m = re_gff3_gene.exec(t[8])) != null) {
				if (m[1] == "gene_id") id = m[2];
				else if (m[1] == "source_gene") src = m[2];
				else if (m[1] == "gene_type") type = m[2];
				else if (m[1] == "gene_biotype") biotype = m[2];
				else if (m[1] == "gene_name") name = m[2];
			}
			if (src != null) id = src;
			if (type == "" && biotype != null) type = biotype;
			print(t[0], parseInt(t[3]) - 1, t[4], [id, type, name].join("|"), 1000, t[6]);
			continue;
		}
		if (t[2] != "CDS" && t[2] != "exon") continue;
		t[3] = parseInt(t[3]) - 1;
		t[4] = parseInt(t[4]);
		var id = null, type = "", name = "N/A", biotype = "", m, tname = "N/A", ens_canonical = false;
		while ((m = re_gtf.exec(t[8])) != null) {
			if (m[1] == "transcript_id") id = m[2];
			else if (m[1] == "transcript_type") type = m[2];
			else if (m[1] == "transcript_biotype" || m[1] == "gbkey") biotype = m[2];
			else if (m[1] == "gene_name" || m[1] == "gene_id") name = m[2];
			else if (m[1] == "transcript_name") tname = m[2];
			else if (m[1] == "tag" && m[2] == "Ensembl_canonical") ens_canonical = true;
		}
		while ((m = re_gff3.exec(t[8])) != null) {
			if (m[1] == "transcript_id") id = m[2];
			else if (m[1] == "transcript_type") type = m[2];
			else if (m[1] == "transcript_biotype" || m[1] == "gbkey") biotype = m[2];
			else if (m[1] == "gene_name" || m[1] == "gene_id") name = m[2];
			else if (m[1] == "transcript_name") tname = m[2];
		}
		if (ens_canon_only && !ens_canonical) continue;
		if (type == "" && biotype != "") type = biotype;
		if (id == null) throw Error("No transcript_id");
		if (id != last_id) {
			print_bed12(exons, cds_st, cds_en, is_short, print_junc);
			exons = [], cds_st = 1<<30, cds_en = 0;
			last_id = id;
		}
		if (t[2] == "CDS") {
			cds_st = cds_st < t[3]? cds_st : t[3];
			cds_en = cds_en > t[4]? cds_en : t[4];
		} else if (t[2] == "exon") {
			if (fn_ucsc_fai != null) {
				if (ens2ucsc[t[0]] != null)
					t[0] = ens2ucsc[t[0]];
				else if (/^[A-Z]+\d+\.\d+$/.test(t[0]))
					t[0] = t[0].replace(/([A-Z]+\d+)\.(\d+)/, "chrUn_$1v$2");
			}
			exons.push([t[0], t[3], t[4], t[6], id, type, name, tname]);
		}
	}
	if (last_id != null)
		print_bed12(exons, cds_st, cds_en, is_short, print_junc);
}

gff2bed(arguments);
