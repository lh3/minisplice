set t po eps co so enh "Helvetica,24.5"

#set size .9,1

#set style line 1 lt 1 lc rgb "#c083e3"
#set style line 2 lt 1 lc rgb "#92cdba"
#set style line 3 lt 1 lc rgb "#b5d9f1"
#
#set style line 1 lt 1 lc rgb "#33A02C"
#set style line 2 lt 1 lc rgb "#FF7F00"
#set style line 3 lt 1 lc rgb "#1F78B4"
#set style line 4 lt 1 lc rgb "#E31A1C"

set border 3

set auto x
set style data histogram
set style histogram cluster gap 1
set style fill solid border -1 #noborder
set boxwidth .9
set xtic rotate by -30 scale 0

set yran [0:*]

set out "dr-hs-hist.eps"
set xlab "No-splice-score protein alignment identity (%)"
set ylab "# zebrafish proteins"
plot "dr-hs-data.txt" u ($9):xtic(1) not ls 2

set out "dr-hs-exon.eps"
set ylab "Fraction of zebrafish junctions aligned"
plot "dr-hs-data.txt" u ($3/$2):xtic(1) t "Without splice score" ls 2, \
	"" u ($5/$2) t "With splice score" ls 3, \
	"" u ($7/$2) t "Spaln3" ls 4

set out "dr-hs-jerr.eps"
set ylab "Junction error rate (%)"
set key top left
plot "dr-hs-data.txt" u (100-$4*100):xtic(1) t "Without splice score" ls 2, \
	"" u (100-$6*100) t "With splice score" ls 3, \
	"" u (100-$8*100) t "Spaln3" ls 4

set out "dRNA004-jerr.eps"
set xlab "No-splice-score dRNA-seq alignment identity (%)"
plot "dRNA004-data.txt" u (100-$3*100):xtic(1) t "Without splice score" ls 2, \
	"" u (100-$5*100) t "With splice score" ls 3
