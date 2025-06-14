set t po eps co "Helvetica,19" lw 1

set xran [0:.1]
set key bot right
set grid

set yran [0.5:1]

set size 0.75,1

set style line 1 lt 1 lc rgb '#A6CEE3' pt 5 # light blue
set style line 2 lt 1 lc rgb '#1F78B4' pt 4 # dark blue
set style line 3 lt 1 lc rgb '#B2DF8A' pt 7 # light green
set style line 4 lt 1 lc rgb '#33A02C' pt 6 # dark green
set style line 5 lt 1 lc rgb '#FB9A99' pt 9 # light red
set style line 6 lt 1 lc rgb '#E31A1C' pt 8 # dark red
set style line 7 lt 1 lc rgb '#FDBF6F' pt 11 # light orange
set style line 8 lt 1 lc rgb '#FF7F00' pt 10 # dark orange
set style line 9 lt 1 lc rgb '#cab2d6' pt 13 # light purple
set style line 10 lt 1 lc rgb '#6a3d9a' pt 12 # dark purple

set out "roc-mod-self.eps"
set title "Trained on odd chr; tested on even chr (model org)"
set xlab "Splice site false positive rate: FP/(FP+TN)"
set ylab "Splice site sensitivity: TP/(TP+FN)"
plot "<grep ^BN human.self.eval" u 10:9 t "m-human (206k,10.4)" w lp ls 1 lw 2, \
	"<grep ^BN mouse.self.eval" u 10:9 t "m-mouse (209k,9.7)" w lp ls 2 lw 2, \
	"<grep ^BN zebrafish.self.eval" u 10:9 t "f-zebrafish (242k,9.6)" w lp ls 8 lw 2, \
	"<grep ^BN fruitfly.self.eval" u 10:9 t "id-Drosophila (60k,4.3)" w lp ls 3 lw 2, \
	"<grep ^BN mosquito.self.eval" u 10:9 t "id-Anopheles (57k,4.5)" w lp ls 4 lw 2, \
	"<grep ^BN honeybee.self.eval" u 10:9 t "iy-honeybee (74k,7.5)" w lp ls 5 lw 2

set out "roc-v-self.eps"
set title "Trained on odd chr; tested on even chr (vertebrate)"
plot "<grep ^BN bAnaAcu-n.self.eval" u 10:9 t "bAnaAcu ncbi (192k,11.3)" w lp ls 1 lw 2, \
	"<grep ^BN bAnaAcu-d.self.eval" u 10:9 t "bAnaAcu dtol (170k,10.7)" w lp ls 2 lw 2, \
	"<grep ^BN fPunPun-n.self.eval" u 10:9 t "fPunPun ncbi (243k,11.1)" w lp ls 7 lw 2, \
	"<grep ^BN fPunPun-d.self.eval" u 10:9 t "fPunPun dtol (224k,10.8)" w lp ls 8 lw 2

set out "roc-i-self.eps"
set title "Trained on odd chr; tested on even chr (insect)"
plot "<grep ^BN icTenMoli-n.self.eval" u 10:9 t "icTenMoli ncbi (84k,5.9)" w lp ls 1 lw 2, \
	"<grep ^BN icTenMoli-d.self.eval" u 10:9 t "icTenMoli dtol (111k,4.1)" w lp ls 2 lw 2, \
	"<grep ^BN idStoCalc-n.self.eval" u 10:9 t "idStoCalc ncbi (70k,4.6)" w lp ls 7 lw 2, \
	"<grep ^BN idStoCalc-d.self.eval" u 10:9 t "idStoCalc dtol (70k,4.4)" w lp ls 8 lw 2, \
	"<grep ^BN ihPlaCitr-n.self.eval" u 10:9 t "ihPlaCitr ncbi (128k,7.7)" w lp ls 3 lw 2, \
	"<grep ^BN ihPlaCitr-d.self.eval" u 10:9 t "ihPlaCitr dtol (128k,6.9)" w lp ls 4 lw 2, \
	"<grep ^BN ilOstNubi-n.self.eval" u 10:9 t "ilOstNubi ncbi (106k,7.6)" w lp ls 5 lw 2, \
	"<grep ^BN ilOstNubi-d.self.eval" u 10:9 t "ilOstNubi dtol (104k,6.2)" w lp ls 6 lw 2, \
	"<grep ^BN iyVesCrab-n.self.eval" u 10:9 t "iyVesCrab ncbi (78k,7.4)" w lp ls 9 lw 2, \
	"<grep ^BN iyVesCrab-d.self.eval" u 10:9 t "iyVesCrab dtol (76k,6.2)" w lp ls 10 lw 2

set out "roc-v-vi2.eps"
set title "Trained on self / vert+ins (vi2); tested on vert"
plot "<grep ^BN human.self.eval" u 10:9 t "human self" w lp ls 1 lw 2, \
	"<grep ^BN human.vi2-7k.eval" u 10:9 t "human vi2" w lp ls 2 lw 2, \
	"<grep ^BN bAnaAcu-n.self.eval" u 10:9 t "bAnaAcu self" w lp ls 7 lw 2, \
	"<grep ^BN bAnaAcu-n.vi2-7k.eval" u 10:9 t "bAnaAcu vi2" w lp ls 8 lw 2, \
	"<grep ^BN aDenEbr-n.self.eval" u 10:9 t "aDenEbr self" w lp ls 3 lw 2, \
	"<grep ^BN aDenEbr-n.vi2-7k.eval" u 10:9 t "aDenEbr vi2" w lp ls 4 lw 2, \
	"<grep ^BN fPunPun-n.self.eval" u 10:9 t "fPunPun self" w lp ls 5 lw 2, \
	"<grep ^BN fPunPun-n.vi2-7k.eval" u 10:9 t "fPunPun vi2" w lp ls 6 lw 2, \
	"<grep ^BN sMobHyp-n.self.eval" u 10:9 t "sMobHyp self" w lp ls 9 lw 2, \
	"<grep ^BN sMobHyp-n.vi2-7k.eval" u 10:9 t "sMobHyp vi2" w lp ls 10 lw 2

set out "roc-v-v2-vi2.eps"
set title "Trained on vert (v2) / vert+ins (vi2); tested on vert"
plot "<grep ^BN human.v2-7k.eval" u 10:9 t "human v2" w lp ls 1 lw 2, \
	"<grep ^BN human.vi2-7k.eval" u 10:9 t "human vi2" w lp ls 2 lw 2, \
	"<grep ^BN bAnaAcu-n.v2-7k.eval" u 10:9 t "bAnaAcu v2" w lp ls 7 lw 2, \
	"<grep ^BN bAnaAcu-n.vi2-7k.eval" u 10:9 t "bAnaAcu vi2" w lp ls 8 lw 2, \
	"<grep ^BN aDenEbr-n.v2-7k.eval" u 10:9 t "aDenEbr v2" w lp ls 3 lw 2, \
	"<grep ^BN aDenEbr-n.vi2-7k.eval" u 10:9 t "aDenEbr vi2" w lp ls 4 lw 2, \
	"<grep ^BN fPunPun-n.v2-7k.eval" u 10:9 t "fPunPun v2" w lp ls 5 lw 2, \
	"<grep ^BN fPunPun-n.vi2-7k.eval" u 10:9 t "fPunPun vi2" w lp ls 6 lw 2, \
	"<grep ^BN sMobHyp-n.v2-7k.eval" u 10:9 t "sMobHyp v2" w lp ls 9 lw 2, \
	"<grep ^BN sMobHyp-n.vi2-7k.eval" u 10:9 t "sMobHyp vi2" w lp ls 10 lw 2

set out "roc-i-vi2.eps"
set title "Trained on self / vert+ins (vi2); tested on ins"
plot "<grep ^BN icTenMoli-n.self.eval" u 10:9 t "*icTenMoli self" w lp ls 1 lw 2, \
	"<grep ^BN icTenMoli-n.vi2-7k.eval" u 10:9 t "*icTenMoli vi2" w lp ls 2 lw 2, \
	"<grep ^BN idStoCalc-n.self.eval" u 10:9 t "idStoCalc self" w lp ls 7 lw 2, \
	"<grep ^BN idStoCalc-n.vi2-7k.eval" u 10:9 t "idStoCalc vi2" w lp ls 8 lw 2, \
	"<grep ^BN ihPlaCitr-n.self.eval" u 10:9 t "*ihPlaCitr self" w lp ls 3 lw 2, \
	"<grep ^BN ihPlaCitr-n.vi2-7k.eval" u 10:9 t "*ihPlaCitr vi2" w lp ls 4 lw 2, \
	"<grep ^BN ilOstNubi-n.self.eval" u 10:9 t "ilOstNubi self" w lp ls 5 lw 2, \
	"<grep ^BN ilOstNubi-n.vi2-7k.eval" u 10:9 t "ilOstNubi vi2" w lp ls 6 lw 2, \
	"<grep ^BN iyVesCrab-n.self.eval" u 10:9 t "iyVesCrab self" w lp ls 9 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2-7k.eval" u 10:9 t "iyVesCrab vi2" w lp ls 10 lw 2

set out "roc-i-i2-vi2.eps"
set title "Trained on ins (i2) / vert+ins (vi2); tested on ins"
plot "<grep ^BN icTenMoli-n.i2-7k.eval" u 10:9 t "*icTenMoli i2" w lp ls 1 lw 2, \
	"<grep ^BN icTenMoli-n.vi2-7k.eval" u 10:9 t "*icTenMoli vi2" w lp ls 2 lw 2, \
	"<grep ^BN idStoCalc-n.i2-7k.eval" u 10:9 t "idStoCalc i2" w lp ls 7 lw 2, \
	"<grep ^BN idStoCalc-n.vi2-7k.eval" u 10:9 t "idStoCalc vi2" w lp ls 8 lw 2, \
	"<grep ^BN ihPlaCitr-n.i2-7k.eval" u 10:9 t "*ihPlaCitr i2" w lp ls 3 lw 2, \
	"<grep ^BN ihPlaCitr-n.vi2-7k.eval" u 10:9 t "*ihPlaCitr vi2" w lp ls 4 lw 2, \
	"<grep ^BN ilOstNubi-n.i2-7k.eval" u 10:9 t "ilOstNubi i2" w lp ls 5 lw 2, \
	"<grep ^BN ilOstNubi-n.vi2-7k.eval" u 10:9 t "ilOstNubi vi2" w lp ls 6 lw 2, \
	"<grep ^BN iyVesCrab-n.i2-7k.eval" u 10:9 t "iyVesCrab i2" w lp ls 9 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2-7k.eval" u 10:9 t "iyVesCrab vi2" w lp ls 10 lw 2

set out "roc-vi2-102.eps"
set title "202bp vs 102bp window"
plot "<grep ^BN human.vi2-7k.eval" u 10:9 t "human 202" w lp ls 1 lw 2, \
	"<grep ^BN human.vi2a-7k.eval" u 10:9 t "human 102" w lp ls 2 lw 2, \
	"<grep ^BN aDenEbr-n.vi2-7k.eval" u 10:9 t "aDenEbr 202" w lp ls 7 lw 2, \
	"<grep ^BN aDenEbr-n.vi2a-7k.eval" u 10:9 t "aDenEbr 102" w lp ls 8 lw 2, \
	"<grep ^BN sMobHyp-n.vi2-7k.eval" u 10:9 t "sMobHyp 202" w lp ls 3 lw 2, \
	"<grep ^BN sMobHyp-n.vi2a-7k.eval" u 10:9 t "sMobHyp 102" w lp ls 4 lw 2, \
	"<grep ^BN idStoCalc-n.vi2-7k.eval" u 10:9 t "idStoCalc 202" w lp ls 5 lw 2, \
	"<grep ^BN idStoCalc-n.vi2a-7k.eval" u 10:9 t "idStoCalc 102" w lp ls 6 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2-7k.eval" u 10:9 t "iyVesCrab 202" w lp ls 9 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2a-7k.eval" u 10:9 t "iyVesCrab 102" w lp ls 10 lw 2

set out "roc-vi2-302.eps"
set title "202bp vs 302bp window"
plot "<grep ^BN human.vi2-7k.eval" u 10:9 t "human 202" w lp ls 1 lw 2, \
	"<grep ^BN human.vi2b-7k.eval" u 10:9 t "human 302" w lp ls 2 lw 2, \
	"<grep ^BN aDenEbr-n.vi2-7k.eval" u 10:9 t "aDenEbr 202" w lp ls 7 lw 2, \
	"<grep ^BN aDenEbr-n.vi2b-7k.eval" u 10:9 t "aDenEbr 302" w lp ls 8 lw 2, \
	"<grep ^BN sMobHyp-n.vi2-7k.eval" u 10:9 t "sMobHyp 202" w lp ls 3 lw 2, \
	"<grep ^BN sMobHyp-n.vi2b-7k.eval" u 10:9 t "sMobHyp 302" w lp ls 4 lw 2, \
	"<grep ^BN idStoCalc-n.vi2-7k.eval" u 10:9 t "idStoCalc 202" w lp ls 5 lw 2, \
	"<grep ^BN idStoCalc-n.vi2b-7k.eval" u 10:9 t "idStoCalc 302" w lp ls 6 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2-7k.eval" u 10:9 t "iyVesCrab 202" w lp ls 9 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2b-7k.eval" u 10:9 t "iyVesCrab 302" w lp ls 10 lw 2
