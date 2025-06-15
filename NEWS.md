Release 0.2-r104 (15 June 2025)
---------------------------

Notable changes:

 * Change: added one more max-pooling layer and removed dropout. This reduced
   the default model size to 7k.

 * Improvement: support NEON for faster matrix multiplication

 * New feature: output values a the last max1d layer and activation differences
   between positive and negative samples.

(0.2: 15 June 2025, r104)
