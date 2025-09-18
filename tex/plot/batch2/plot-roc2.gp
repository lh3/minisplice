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
plot "<grep ^BN human.self.eval" u 10:9 t "human self (0.82)" w lp ls 1 lw 2, \
	"<grep ^BN human.vi2-7k.eval" u 10:9 t "human vi2 (0.79)" w lp ls 2 lw 2, \
	"<grep ^BN bAnaAcu-n.self.eval" u 10:9 t "bAnaAcu self (0.86)" w lp ls 7 lw 2, \
	"<grep ^BN bAnaAcu-n.vi2-7k.eval" u 10:9 t "bAnaAcu vi2 (0.86)" w lp ls 8 lw 2, \
	"<grep ^BN aDenEbr-n.self.eval" u 10:9 t "aDenEbr self (0.91)" w lp ls 3 lw 2, \
	"<grep ^BN aDenEbr-n.vi2-7k.eval" u 10:9 t "aDenEbr vi2 (0.88)" w lp ls 4 lw 2, \
	"<grep ^BN fPunPun-n.self.eval" u 10:9 t "fPunPun self (0.91)" w lp ls 5 lw 2, \
	"<grep ^BN fPunPun-n.vi2-7k.eval" u 10:9 t "fPunPun vi2 (0.86)" w lp ls 6 lw 2, \
	"<grep ^BN sMobHyp-n.self.eval" u 10:9 t "sMobHyp self (0.91)" w lp ls 9 lw 2, \
	"<grep ^BN sMobHyp-n.vi2-7k.eval" u 10:9 t "sMobHyp vi2 (0.88)" w lp ls 10 lw 2

set out "roc-v-v2-vi2.eps"
set title "Trained on vert (v2) / vert+ins (vi2); tested on vert"
plot "<grep ^BN human.v2-7k.eval" u 10:9 t "human v2 (0.81)" w lp ls 1 lw 2, \
	"<grep ^BN human.vi2-7k.eval" u 10:9 t "human vi2 (0.79)" w lp ls 2 lw 2, \
	"<grep ^BN bAnaAcu-n.v2-7k.eval" u 10:9 t "bAnaAcu v2 (0.88)" w lp ls 7 lw 2, \
	"<grep ^BN bAnaAcu-n.vi2-7k.eval" u 10:9 t "bAnaAcu vi2 (0.86)" w lp ls 8 lw 2, \
	"<grep ^BN aDenEbr-n.v2-7k.eval" u 10:9 t "aDenEbr v2 (0.88)" w lp ls 3 lw 2, \
	"<grep ^BN aDenEbr-n.vi2-7k.eval" u 10:9 t "aDenEbr vi2 (0.88)" w lp ls 4 lw 2, \
	"<grep ^BN fPunPun-n.v2-7k.eval" u 10:9 t "fPunPun v2 (0.87)" w lp ls 5 lw 2, \
	"<grep ^BN fPunPun-n.vi2-7k.eval" u 10:9 t "fPunPun vi2 (0.86)" w lp ls 6 lw 2, \
	"<grep ^BN sMobHyp-n.v2-7k.eval" u 10:9 t "sMobHyp v2 (0.89)" w lp ls 9 lw 2, \
	"<grep ^BN sMobHyp-n.vi2-7k.eval" u 10:9 t "sMobHyp vi2 (0.88)" w lp ls 10 lw 2

set out "roc-i-vi2.eps"
set title "Trained on self / vert+ins (vi2); tested on ins"
plot "<grep ^BN icTenMoli-n.self.eval" u 10:9 t "*icTenMoli self (0.82)" w lp ls 1 lw 2, \
	"<grep ^BN icTenMoli-n.vi2-7k.eval" u 10:9 t "*icTenMoli vi2 (0.77)" w lp ls 2 lw 2, \
	"<grep ^BN idStoCalc-n.self.eval" u 10:9 t "idStoCalc self (0.84)" w lp ls 7 lw 2, \
	"<grep ^BN idStoCalc-n.vi2-7k.eval" u 10:9 t "idStoCalc vi2 (0.76)" w lp ls 8 lw 2, \
	"<grep ^BN ihPlaCitr-n.self.eval" u 10:9 t "*ihPlaCitr self (0.89)" w lp ls 3 lw 2, \
	"<grep ^BN ihPlaCitr-n.vi2-7k.eval" u 10:9 t "*ihPlaCitr vi2 (0.83)" w lp ls 4 lw 2, \
	"<grep ^BN ilOstNubi-n.self.eval" u 10:9 t "ilOstNubi self (0.89)" w lp ls 5 lw 2, \
	"<grep ^BN ilOstNubi-n.vi2-7k.eval" u 10:9 t "ilOstNubi vi2 (0.86)" w lp ls 6 lw 2, \
	"<grep ^BN iyVesCrab-n.self.eval" u 10:9 t "iyVesCrab self (0.84)" w lp ls 9 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2-7k.eval" u 10:9 t "iyVesCrab vi2 (0.78)" w lp ls 10 lw 2

set out "roc-i-i2-vi2.eps"
set title "Trained on ins (i2) / vert+ins (vi2); tested on ins"
plot "<grep ^BN icTenMoli-n.i2-7k.eval" u 10:9 t "*icTenMoli i2 (0.81)" w lp ls 1 lw 2, \
	"<grep ^BN icTenMoli-n.vi2-7k.eval" u 10:9 t "*icTenMoli vi2 (0.77)" w lp ls 2 lw 2, \
	"<grep ^BN idStoCalc-n.i2-7k.eval" u 10:9 t "idStoCalc i2 (0.77)" w lp ls 7 lw 2, \
	"<grep ^BN idStoCalc-n.vi2-7k.eval" u 10:9 t "idStoCalc vi2 (0.76)" w lp ls 8 lw 2, \
	"<grep ^BN ihPlaCitr-n.i2-7k.eval" u 10:9 t "*ihPlaCitr i2 (0.85)" w lp ls 3 lw 2, \
	"<grep ^BN ihPlaCitr-n.vi2-7k.eval" u 10:9 t "*ihPlaCitr vi2 (0.83)" w lp ls 4 lw 2, \
	"<grep ^BN ilOstNubi-n.i2-7k.eval" u 10:9 t "ilOstNubi i2 (0.87)" w lp ls 5 lw 2, \
	"<grep ^BN ilOstNubi-n.vi2-7k.eval" u 10:9 t "ilOstNubi vi2 (0.86)" w lp ls 6 lw 2, \
	"<grep ^BN iyVesCrab-n.i2-7k.eval" u 10:9 t "iyVesCrab i2 (0.81)" w lp ls 9 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2-7k.eval" u 10:9 t "iyVesCrab vi2 (0.78)" w lp ls 10 lw 2

set out "roc-vi2-102.eps"
set title "202bp vs 102bp window"
plot "<grep ^BN human.vi2-7k.eval" u 10:9 t "human 202 (0.79)" w lp ls 1 lw 2, \
	"<grep ^BN human.vi2a-7k.eval" u 10:9 t "human 102 (0.75)" w lp ls 2 lw 2, \
	"<grep ^BN aDenEbr-n.vi2-7k.eval" u 10:9 t "aDenEbr 202 (0.88)" w lp ls 7 lw 2, \
	"<grep ^BN aDenEbr-n.vi2a-7k.eval" u 10:9 t "aDenEbr 102 (0.84)" w lp ls 8 lw 2, \
	"<grep ^BN sMobHyp-n.vi2-7k.eval" u 10:9 t "sMobHyp 202 (0.88)" w lp ls 3 lw 2, \
	"<grep ^BN sMobHyp-n.vi2a-7k.eval" u 10:9 t "sMobHyp 102 (0.84)" w lp ls 4 lw 2, \
	"<grep ^BN idStoCalc-n.vi2-7k.eval" u 10:9 t "idStoCalc 202 (0.76)" w lp ls 5 lw 2, \
	"<grep ^BN idStoCalc-n.vi2a-7k.eval" u 10:9 t "idStoCalc 102 (0.69)" w lp ls 6 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2-7k.eval" u 10:9 t "iyVesCrab 202 (0.78)" w lp ls 9 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2a-7k.eval" u 10:9 t "iyVesCrab 102 (0.71)" w lp ls 10 lw 2

set out "roc-vi2-302.eps"
set title "202bp vs 302bp window"
plot "<grep ^BN human.vi2-7k.eval" u 10:9 t "human 202 (0.79)" w lp ls 1 lw 2, \
	"<grep ^BN human.vi2b-7k.eval" u 10:9 t "human 302 (0.81)" w lp ls 2 lw 2, \
	"<grep ^BN aDenEbr-n.vi2-7k.eval" u 10:9 t "aDenEbr 202 (0.88)" w lp ls 7 lw 2, \
	"<grep ^BN aDenEbr-n.vi2b-7k.eval" u 10:9 t "aDenEbr 302 (0.88)" w lp ls 8 lw 2, \
	"<grep ^BN sMobHyp-n.vi2-7k.eval" u 10:9 t "sMobHyp 202 (0.88)" w lp ls 3 lw 2, \
	"<grep ^BN sMobHyp-n.vi2b-7k.eval" u 10:9 t "sMobHyp 302 (0.90)" w lp ls 4 lw 2, \
	"<grep ^BN idStoCalc-n.vi2-7k.eval" u 10:9 t "idStoCalc 202 (0.76)" w lp ls 5 lw 2, \
	"<grep ^BN idStoCalc-n.vi2b-7k.eval" u 10:9 t "idStoCalc 302 (0.76)" w lp ls 6 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2-7k.eval" u 10:9 t "iyVesCrab 202 (0.78)" w lp ls 9 lw 2, \
	"<grep ^BN iyVesCrab-n.vi2b-7k.eval" u 10:9 t "iyVesCrab 302 (0.80)" w lp ls 10 lw 2
