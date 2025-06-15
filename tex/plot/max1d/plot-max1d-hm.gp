set t po eps co enh "Helvetica,20"

set size 1.23,0.8

set tic scale 0
set xran [-0.5:31.5]
set yran [-0.5:15.5]

set cbrange [-0.6:0.6]

# reverse of 11-class RdBu
set palette defined (\
	0 "#053061", \
	1 "#2166ac", \
	2 "#4393c3", \
	3 "#92c5de", \
	4 "#d1e5f0", \
	5 "#f7f7f7", \
	6 "#fddbc7", \
	7 "#f4a582", \
	8 "#d6604d", \
	9 "#b2182b", \
	10 "#67001f" \
);

# adapted from DeepSeek solution

set view map
set style fill solid noborder  # Solid colored boxes

set ylab "Feature"
set ytic 1 font "Helvetica,18"
set cblabel "Activation difference"

set out "max1d-hm-D.eps"
set xtics ("0" 15.5, "+50" 10.5, "+100" 5.5, "-50" 20.5, "-100" 25.5, "+150" -0.5, "-150" 31.5)
set xlab "Approximate position on donor sequences (bp); from exon to intron"
plot '<grep ^MA vi2b-7k-D.max1d | cut -f4- | ./transpose.js -' matrix u 1:2:(0.5):(0.5):($3/100) not with boxxy fc palette

set out "max1d-hm-A.eps"
set xtics ("0" 15.5, "-50" 10.5, "-100" 5.5, "+50" 20.5, "+100" 25.5, "-150" -0.5, "+150" 31.5)
set xlab "Approximate position on acceptor sequences (bp); from intron to exon"
plot '<grep ^MA vi2b-7k-A.max1d | cut -f4- | ./transpose.js -' matrix u 1:2:(0.5):(0.5):($3/100) not with boxxy fc palette
