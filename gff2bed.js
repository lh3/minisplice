#!/usr/bin/env k8

"use strict";

/***************
 * k8 routines *
 ***************/

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

/*****************
 * main routines *
 *****************/

const g2b_colors = {
	'protein_coding':'0,128,255',
	'mRNA':'0,128,255',
	'lincRNA':'0,192,0',
	'snRNA':'0,192,0',
	'miRNA':'0,192,0',
	'misc_RNA':'0,192,0'
};

class Transcript {
	constructor(tid, type, ctg, strand, gid, gname) {
		this.tid = tid, this.type = type, this.ctg = ctg, this.strand = strand, this.gid = gid, this.gname = gname;
		this.cds_st = -1, this.cds_en = -1;
		this.exon = [];
	}
	row() {
		const disp_name = [this.tid, this.type, this.gname].join("|");
		const a = this.exon.sort(function(a,b) {return a.st - b.st});
		const st = a[0].st;
		const en = a[a.length - 1].en;
		const cds_st = this.cds_st >= 0? this.cds_st : st;
		const cds_en = this.cds_en >= 0? this.cds_en : en;
		if (cds_st < st || cds_en > en)
			throw Error(`inconsistent thick start or end for transcript ${this.tid}`);
		let sizes = [], starts = [];
		for (let i = 0; i < a.length; ++i) {
			if (i > 0 && a[i].st <a[i-1].en)
				throw Error(`overlapping exon for transcript ${this.tid}`);
			sizes.push(a[i].en - a[i].st);
			starts.push(a[i].st - st);
		}
		let color = g2b_colors[this.type];
		if (color == null) color = "196,196,196";
		return [this.ctg, st, en, disp_name, 1000, this.strand, cds_st, cds_en, color, a.length, sizes.join(",")+",", starts.join(",")+","];
	}
}

function g2b_print(gene, longest)
{
	if (gene.length == 0) return;
	if (longest) {
		let max = 0, max_i = -1;
		for (let i = 0; i < gene.length; ++i) {
			let len = 0, t = gene[i];
			for (let j = 0; j < t.exon.length; ++j)
				len += t.exon[j].en - t.exon[j].st;
			if (max < len) max = len, max_i = i;
		}
		print(gene[max_i].row().join("\t"));
	} else {
		for (let i = 0; i < gene.length; ++i)
			print(gene[i].row().join("\t"));
	}
}

function g2b_main(args)
{
	let fn_ucsc_fai = null, ens_canon_only = false, longest = false, protein_only = false;
	for (const o of getopt(args, "u:elp")) {
		if (o.opt == '-u') fn_ucsc_fai = o.arg;
		else if (o.opt == '-l') longest = true;
		else if (o.opt == '-e') ens_canon_only = true;
		else if (o.opt == '-p') protein_only = true;
	}

	if (args.length == 0) {
		print("Usage: gff2bed.js [options] <in.gff>");
		print("Options:");
		print("  -l         choose the longest transcript for each gene");
		print("  -e         only include transcripts tagged with 'Ensembl_canonical'");
		print("  -p         only include 'protein_coding' transcripts");
		print("  -u FILE    hg38.fa.fai for chr name conversion");
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

	const re_gff = /\b(transcript_id|transcript_type|transcript_biotype|gene_name|gene_id|gbkey|tag)( "([^"]+)"|=([^;]+));/g;
	const re_gff_gene = /\b(gene_id|gene_type|source_gene|gene_biotype|gene_name)( "([^;]+)"|=([^;]+));/g;

	let gene = [];
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		if (t[0].charAt(0) == '#') continue;
		if (t[2] != "CDS" && t[2] != "exon") continue;

		// update contig name if requested
		let ctg = t[0];
		if (fn_ucsc_fai != null) {
			if (ens2ucsc[ctg] != null)
				tr.ctg = ens2ucsc[ctg];
			else if (/^[A-Z]+\d+\.\d+$/.test(ctg))
				ctg = ctg.replace(/([A-Z]+\d+)\.(\d+)/, "chrUn_$1v$2");
		}

		// collect meta annotation
		const st = parseInt(t[3]) - 1;
		const en = parseInt(t[4]);
		let m, tid = null, gid = null, type = "", gname = null, biotype = "", ens_canonical = false;
		while ((m = re_gff.exec(t[8])) != null) {
			const key = m[1], val = m[3]? m[3] : m[4];
			if (key == "transcript_id") tid = val;
			else if (key == "transcript_type") type = val;
			else if (key == "transcript_biotype" || key == "gbkey") biotype = val;
			else if (key == "gene_name") gname = val;
			else if (key == "gene_id") gid = val;
			else if (key == "tag" && val == "Ensembl_canonical") ens_canonical = true;
		}
		if (ens_canon_only && !ens_canonical) continue;
		if (tid == null) throw Error("No transcript_id");
		if (gname == null) gname = gid? gid : "*"; // infer gene name
		if (gid == null) gid = gname; // if gene_id missing, use gene name to identify a gene
		if (type == "" && biotype != "") type = biotype; // infer transcript type
		if (protein_only && type != "" && type != "protein_coding") continue;

		// finish transcript or gene
		if (gene.length == 0 || gene[gene.length - 1].tid != tid) { // changing transcript
			if (gene.length > 0 && gene[gene.length - 1].gid != gid) { // changing gene
				g2b_print(gene, longest);
				gene = [];
			}
			gene.push(new Transcript(tid, type, ctg, t[6], gid, gname));
		}

		// update the last transcript
		let tr = gene[gene.length - 1];
		if (t[2] == "CDS") {
			tr.cds_st = tr.cds_st >= 0 && tr.cds_st < st? tr.cds_st : st;
			tr.cds_en = tr.cds_en >= 0 && tr.cds_en > en? tr.cds_en : en;
		} else if (t[2] == "exon") {
			tr.exon.push({ st:st, en:en });
		}
	}
	g2b_print(gene, longest);
}

g2b_main(arguments);
