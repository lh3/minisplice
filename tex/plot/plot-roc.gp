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
set style line 9 lt 1 lc rgb '#cab2d6' pt 13
set style line 10 lt 1 lc rgb '#6a3d9a' pt 12

set out "human-to-basic.eps"
set title "Trained on different GenCode annotation"
set xlab "Donor false positive rate: FP/(FP+TN)"
set ylab "Donor sensitivity: TP/(TP+FN)"
plot "<grep ^BN human/human.plong.to-basic-D.eval" u 10:9 t "Longest protein (0.86)" w lp ls 2 lw 2, \
	"<grep ^BN human/human.basic.to-basic-D.eval" u 10:9 t "Basic (0.87)" w lp ls 8 lw 2, \
	"<grep ^BN human/human.cmphs.to-basic-D.eval" u 10:9 t "Comphrehensive (0.87)" w lp ls 4 lw 2

set out "human-from-plong.eps"
set title "Tested on different GenCode annotation"
plot "<grep ^BN human/human.plong.from-plong-D.eval" u 10:9 t "Longest protein (0.94)" w lp ls 2 lw 2, \
	"<grep ^BN human/human.basic.from-plong-D.eval" u 10:9 t "Basic (0.86)" w lp ls 8 lw 2, \
	"<grep ^BN human/human.cmphs.from-plong-D.eval" u 10:9 t "Comphrehensive (0.81)" w lp ls 4 lw 2

#set out "self-vi1.eps"
#set title "Intra- vs cross-species training"
#set xlab "Splice site false positive rate: FP/(FP+TN)"
#set ylab "Splice site sensitivity: TP/(TP+FN)"
#plot "<grep ^BN self/human.self.eval" u 10:9 t "Human intra" w lp ls 1 lw 2, \
#	"<grep ^BN vi1/human.bi1.eval" u 10:9 t "Human cross" w lp ls 2 lw 2, \
#	"<grep ^BN self/zebrafinch.self.eval" u 10:9 t "Zebrafinch intra" w lp ls 7 lw 2, \
#	"<grep ^BN vi1/zebrafinch.bi1.eval" u 10:9 t "Zebrafinch cross" w lp ls 8 lw 2, \
#	"<grep ^BN self/fugu.self.eval" u 10:9 t "Fugu intra" w lp ls 3 lw 2, \
#	"<grep ^BN vi1/fugu.bi1.eval" u 10:9 t "Fugu cross" w lp ls 4 lw 2, \
#	"<grep ^BN self/honeybee.self.eval" u 10:9 t "Honeybee intra" w lp ls 5 lw 2, \
#	"<grep ^BN vi1/honeybee.bi1.eval" u 10:9 t "Honeybee cross" w lp ls 6 lw 2

#set out "vi1-vi2.eps"
#set title "Training with 202bp vs 102bp window"
#set xlab "Splice site false positive rate: FP/(FP+TN)"
#set ylab "Splice site sensitivity: TP/(TP+FN)"
#plot "<grep ^BN vi1/human.bi1.eval" u 10:9 t "Human 202bp" w lp ls 1 lw 2, \
#	"<grep ^BN vi2/human.bi2.eval" u 10:9 t "Human 102bp" w lp ls 2 lw 2, \
#	"<grep ^BN vi1/chicken.bi1.eval" u 10:9 t "Chicken 202bp" w lp ls 7 lw 2, \
#	"<grep ^BN vi2/chicken.bi2.eval" u 10:9 t "Chicken 102bp" w lp ls 8 lw 2, \
#	"<grep ^BN vi1/fugu.bi1.eval" u 10:9 t "Fugu 202bp" w lp ls 3 lw 2, \
#	"<grep ^BN vi2/fugu.bi2.eval" u 10:9 t "Fugu 102bp" w lp ls 4 lw 2, \
#	"<grep ^BN vi1/mosquito.bi1.eval" u 10:9 t "Mosquito 202bp" w lp ls 5 lw 2, \
#	"<grep ^BN vi2/mosquito.bi2.eval" u 10:9 t "Mosquito 102bp" w lp ls 6 lw 2

#set out "vi1-vi5.eps"
#set title "Large (139k parameters) vs small model (35k)"
#set xlab "Splice site false positive rate: FP/(FP+TN)"
#set ylab "Splice site sensitivity: TP/(TP+FN)"
#plot "<grep ^BN vi1/human.bi1.eval" u 10:9 t "Human large" w lp ls 1 lw 2, \
#	"<grep ^BN vi5/human.vi5.eval" u 10:9 t "Human small" w lp ls 2 lw 2, \
#	"<grep ^BN vi1/chicken.bi1.eval" u 10:9 t "Chicken large" w lp ls 7 lw 2, \
#	"<grep ^BN vi5/chicken.vi5.eval" u 10:9 t "Chicken small" w lp ls 8 lw 2, \
#	"<grep ^BN vi1/fugu.bi1.eval" u 10:9 t "Fugu large" w lp ls 3 lw 2, \
#	"<grep ^BN vi5/fugu.vi5.eval" u 10:9 t "Fugu small" w lp ls 4 lw 2, \
#	"<grep ^BN vi1/mosquito.bi1.eval" u 10:9 t "Mosquito large" w lp ls 5 lw 2, \
#	"<grep ^BN vi5/mosquito.vi5.eval" u 10:9 t "Mosquito small" w lp ls 6 lw 2

#set out "vi1-vi6.eps"
#set title "Training with 202bp vs 302bp window"
#set xlab "Splice site false positive rate: FP/(FP+TN)"
#set ylab "Splice site sensitivity: TP/(TP+FN)"
#plot "<grep ^BN vi1/human.bi1.eval" u 10:9 t "Human 202bp" w lp ls 1 lw 2, \
#	"<grep ^BN vi6/human.vi6.eval" u 10:9 t "Human 302bp" w lp ls 2 lw 2, \
#	"<grep ^BN vi1/chicken.bi1.eval" u 10:9 t "Chicken 202bp" w lp ls 7 lw 2, \
#	"<grep ^BN vi6/chicken.vi6.eval" u 10:9 t "Chicken 302bp" w lp ls 8 lw 2, \
#	"<grep ^BN vi1/fugu.bi1.eval" u 10:9 t "Fugu 202bp" w lp ls 3 lw 2, \
#	"<grep ^BN vi6/fugu.vi6.eval" u 10:9 t "Fugu 302bp" w lp ls 4 lw 2, \
#	"<grep ^BN vi1/mosquito.bi1.eval" u 10:9 t "Mosquito 202bp" w lp ls 5 lw 2, \
#	"<grep ^BN vi6/mosquito.vi6.eval" u 10:9 t "Mosquito 302bp" w lp ls 6 lw 2

#set out "self-D.eps"
#set title "Trained on odd chr; tested on even chr (intra-species)"
#set xlab "Donor false positive rate: FP/(FP+TN)"
#set ylab "Donor sensitivity: TP/(TP+FN)"
#plot "<grep ^BN self/human.selfD.eval" u 10:9 t "Human" w lp ls 2 lw 2, \
#	"<grep ^BN self/mouse.selfD.eval" u 10:9 t "Mouse" w lp ls 1 lw 2, \
#	"<grep ^BN self/chicken.selfD.eval" u 10:9 t "Chicken" w lp ls 8 lw 2, \
#	"<grep ^BN self/zebrafish.selfD.eval" u 10:9 t "Zebrafish" w lp ls 4 lw 2, \
#	"<grep ^BN self/fugu.selfD.eval" u 10:9 t "Fugu" w lp ls 3 lw 2, \
#	"<grep ^BN self/fruitfly.selfD.eval" u 10:9 t "Fruitfly" w lp ls 6 lw 2, \
#	"<grep ^BN self/mosquito.selfD.eval" u 10:9 t "Mosquito" w lp ls 5 lw 2

#set out "self-A.eps"
#set xlab "Acceptor false positive rate: FP/(FP+TN)"
#set ylab "Acceptor sensitivity: TP/(TP+FN)"
#plot "<grep ^BN self/human.selfA.eval" u 10:9 t "Human" w lp ls 2 lw 2, \
#	"<grep ^BN self/mouse.selfA.eval" u 10:9 t "Mouse" w lp ls 1 lw 2, \
#	"<grep ^BN self/chicken.selfA.eval" u 10:9 t "Chicken" w lp ls 8 lw 2, \
#	"<grep ^BN self/zebrafish.selfA.eval" u 10:9 t "Zebrafish" w lp ls 4 lw 2, \
#	"<grep ^BN self/fugu.selfA.eval" u 10:9 t "Fugu" w lp ls 3 lw 2, \
#	"<grep ^BN self/fruitfly.selfA.eval" u 10:9 t "Fruitfly" w lp ls 6 lw 2, \
#	"<grep ^BN self/mosquito.selfA.eval" u 10:9 t "Mosquito" w lp ls 5 lw 2

set out "self-BD.eps"
set title "Donor-acceptor joint training vs donor-only"
set xlab "Donor false positive rate: FP/(FP+TN)"
set ylab "Donor sensitivity: TP/(TP+FN)"
plot "<grep ^BN self/human.selfD.eval" u 10:9 t "Human donor (0.86)" w lp ls 1 lw 2, \
	"<grep ^BN self/human.selfSD.eval" u 10:9 t "Human joint (0.85)" w lp ls 2 lw 2, \
	"<grep ^BN self/chicken.selfD.eval" u 10:9 t "Chicken donor (0.81)" w lp ls 7 lw 2, \
	"<grep ^BN self/chicken.selfSD.eval" u 10:9 t "Chicken joint (0.81)" w lp ls 8 lw 2, \
	"<grep ^BN self/zebrafish.selfD.eval" u 10:9 t "Zebrafish donor (0.90)" w lp ls 3 lw 2, \
	"<grep ^BN self/zebrafish.selfSD.eval" u 10:9 t "Zebrafish joint (0.90)" w lp ls 4 lw 2, \
	"<grep ^BN self/mosquito.selfD.eval" u 10:9 t "Mosquito donor (0.87)" w lp ls 5 lw 2, \
	"<grep ^BN self/mosquito.selfSD.eval" u 10:9 t "Mosquito joint (0.87)" w lp ls 6 lw 2

set out "self-BA.eps"
set title "Donor-acceptor joint training vs acceptor-only"
set xlab "Acceptor false positive rate: FP/(FP+TN)"
set ylab "Acceptor sensitivity: TP/(TP+FN)"
plot "<grep ^BN self/human.selfA.eval" u 10:9 t "Human acceptor (0.83)" w lp ls 1 lw 2, \
	"<grep ^BN self/human.selfSA.eval" u 10:9 t "Human joint (0.82)" w lp ls 2 lw 2, \
	"<grep ^BN self/chicken.selfA.eval" u 10:9 t "Chicken acceptor (0.78)" w lp ls 7 lw 2, \
	"<grep ^BN self/chicken.selfSA.eval" u 10:9 t "Chicken joint (0.77)" w lp ls 8 lw 2, \
	"<grep ^BN self/zebrafish.selfA.eval" u 10:9 t "Zebrafish acceptor (0.89)" w lp ls 3 lw 2, \
	"<grep ^BN self/zebrafish.selfSA.eval" u 10:9 t "Zebrafish joint (0.88)" w lp ls 4 lw 2, \
	"<grep ^BN self/mosquito.selfA.eval" u 10:9 t "Mosquito acceptor (0.88)" w lp ls 5 lw 2, \
	"<grep ^BN self/mosquito.selfSA.eval" u 10:9 t "Mosquito joint (0.87)" w lp ls 6 lw 2

set out "cross-human.eps"
set title "Trained on different species; tested on human"
set xlab "Splice site false positive rate: FP/(FP+TN)"
set ylab "Splice site sensitivity: TP/(TP+FN)"
plot "<grep ^BN cross/human-human.eval" u 10:9 t "Human (0.83)" w lp ls 2 lw 2, \
	"<grep ^BN cross/human-mouse.eval" u 10:9 t "Mouse (0.83)" w lp ls 1 lw 2, \
	"<grep ^BN cross/human-chicken.eval" u 10:9 t "Chicken (0.79)" w lp ls 8 lw 2, \
	"<grep ^BN cross/human-zebrafish.eval" u 10:9 t "Zebrafish (0.66)" w lp ls 4 lw 2, \
	"<grep ^BN cross/human-fruitfly.eval" u 10:9 t "Fruitfly (0.51)" w lp ls 6 lw 2, \
	"<grep ^BN cross/human-mosquito.eval" u 10:9 t "Mosquito (0.54)" w lp ls 5 lw 2

set out "cross-fish.eps"
set title "Cross-species: zebrafish vs fugu"
set xlab "Splice site false positive rate: FP/(FP+TN)"
set ylab "Splice site sensitivity: TP/(TP+FN)"
plot "<grep ^BN self/zebrafish.self.eval" u 10:9 t "Zebrafish-to-Zebrafish (0.89)" w lp ls 1 lw 2, \
	"<grep ^BN cross/zebrafish-fugu.eval" u 10:9 t "Fugu-to-Zebrafish (0.87)" w lp ls 2 lw 2, \
	"<grep ^BN self/fugu.self.eval" u 10:9 t "Fugu-to-Fugu (0.67)" w lp ls 7 lw 2, \
	"<grep ^BN cross/fugu-zebrafish.eval" u 10:9 t "Zebrafish-to-Fugu (0.62)" w lp ls 8 lw 2

set out "cross-insect.eps"
set title "Cross-species: fruitfly vs mosquito"
plot "<grep ^BN self/fruitfly.self.eval" u 10:9 t "Fruitfly-to-fruitfly" w lp ls 1 lw 2, \
	"<grep ^BN cross/fruitfly-mosquito.eval" u 10:9 t "Mosquito-to-fruitfly" w lp ls 2 lw 2, \
	"<grep ^BN self/mosquito.self.eval" u 10:9 t "Mosquito-to-mosquito" w lp ls 7 lw 2, \
	"<grep ^BN cross/mosquito-fruitfly.eval" u 10:9 t "Fruitfly-to-mosquito" w lp ls 8 lw 2
