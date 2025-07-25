Release 0.4-r118 (15 July 2025)
-------------------------------

Notable changes:

 * Bugfix: minisplice supports contigs longer than 2Gb.

 * Change: don't output scores -6 or lower. This reduces splice score file
   sizes by 30% at almost no cost of miniprot accuracy.

(0.4: 15 July 2025, r118)



Release 0.3-r110 (5 July 2025)
------------------------------

Notable changes:

 * Improvement: automatically load the .cali file based on the name of the
   model file.

 * Bugfix: fixed a potential off-by-one error but apparently this has not
   happened.

(0.3: 5 July 2025, r110)



Release 0.2-r104 (15 June 2025)
-------------------------------

Notable changes:

 * Change: added one more max-pooling layer and removed dropout. This reduced
   the default model size to 7k.

 * Improvement: support NEON for faster matrix multiplication

 * New feature: output values a the last max1d layer and activation differences
   between positive and negative samples.

(0.2: 15 June 2025, r104)
