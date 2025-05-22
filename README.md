## Getting Started
```sh
# compile
git clone https://github.com/lh3/minisplice
cd minisplice && make

# download vertebrate-insect pre-trained model and calibration data
wget https://zenodo.org/records/15492781/files/vi1-35k.kan
wget https://zenodo.org/records/15492781/files/vi1-35k.kan.cali

# compute the splice score for GT and AG sites; see below for model training
./minisplice predict -t16 -c vi1-35k.kan.cali vi1-35k.kan genome.fa.gz > score.tsv

# use splice scores (miniprot r272+ recommended)
miniprot -Iut16 --gff -j2 --spsc=score.tsv genome.fa.gz proteins.faa > align.gff

# use pre-calculated human or Drosophila scores
wget https://zenodo.org/records/15492781/files/human-GRCh38.spsc.tsv.gz
miniprot -Iut16 --gff -j2 --spsc=human-GRCh38.spsc.tsv.gz hg38.fa proteins.faa
```

## Introduction

**What:** minisplice is a command-line tool to estimate the odds-ratio score of
canonical donor (GT) and acceptor (AG) splice sites. It is intended to be used
with [miniprot][mp] (r272+) for improving alignment accuracy especially for
distant homologs. Pre-trained models and pre-computed splice scores can be found
[at Zenodo][zn].

**Why:** protein-to-genome aligners like miniprot and GeneWise are effectively
gene-finders that trace open-reading frames and model splice signals. For
distant homologs, the splice model plays an important role in resolving ambiguous
alignment around splice junctions. Miniprot uses a simplistic model with 4-5
parameters based on human data. While this model is reasonably robust in
practice, it losses power in comparison to more sophisticated solutions such as
position weight matrix.

**How:** minisplice trains a 1D convolutional neural network (1D-CNN) on
sequences around annotated and random GT- or -AG sites from the genome,
calibrates the model output to empirical probability, and computes the
odds-ratio score of each GT- or -AG in the genome. Training can be applied to
multiple distantly related species. For example, **one** model trained from mouse,
chicken, zebrafish, *Drosophila* (fruit fly) and *Anopheles* (mosquito) apparently
works well for vertebrtes and insect.

## Usage

### Prediction

If your target genome is a vertebrate or insect, you can use pre-trained model:
```sh
./minisplice predict -t16 -c vi1-35k.kan.cali vi1-35k.kan genome.fa.gz > score.tsv
```
where `vi1-35k.kan` encodes the model and `vi1-35k.kan.cali` provides calibration data
which is used to translate the model output to empirical probability. The
output looks like:
```txt
2L  1005295  +   A   2
2L  1005339  -   D   13
2L  1005339  -   A   0
2L  1005410  -   A   -5
2L  1005415  -   A   -5
```
Each line gives contig/chromosome name, offset, strand, D (for donor) or
A (for acceptor) and the score, which is 2log2-scaled odds ratio of the
estimated probability of the site being real over the genome-wide fraction of
annotated splice sites (the null model). Miniprot can take this file as input
with option `--spsc`. Note that this option was added in miniprot-0.14 (r265),
but versions before r271 may lead to an assertion failure.

### Training

The following command lines show how to train and calibrate a model for one genome:
```sh
# convert gene annotation in GTF/GFF3 to BED12
script/gff2bed.js -pl anno.gtf.gz | gzip > anno-long.bed.gz  # longest protein-coding only
script/gff2bed.js anno.gtf.gz | gzip > anno-all.bed.gz       # all annotation

# generate training data
./minisplice gentrain anno-long.bed.gz genome-odd.fa.gz | gzip > train.txt.gz

# model training; 8 or 16 threads are recommended
./minisplice train -t16 -o model.kan train.txt.gz

# calibration (computing empirical odds-ratio scores)
./minisplice predict -t16 -b anno-all.bed.gz model.kan genome-even.fa.gz > model.cali

# prediction
./minisplice predict -t16 -c model.cali model.kan genome.fa.gz | gzip > score.tsv.gz
```
You can train on odd chromosomes and calibrate on even chromosomes. To train a
model from multiple species:
```sh
./minisplice gentrain anno1-long.bed.gz genome1-odd.fa.gz | gzip > train1.txt.gz
./minisplice gentrain anno2-long.bed.gz genome2-odd.fa.gz | gzip > train2.txt.gz
cat train*.txt.gz | ./minisplice train -t16 -o model.kan -
./minisplice predict -t16 -b anno1-all.bed.gz model.kan genome1-even.fa.gz > cali1.txt
./minisplice predict -t16 -b anno1-all.bed.gz model.kan genome2-even.fa.gz > cali2.txt
script/merge-cali.js 1,1 cali1.txt cali2.txt > model.cali
```
Usually one genome provides enough training data. To save training time, you
can subsample training data from each genome before combining them. If you
subsample training data at different rates, it is recommended to provides the
rates on the `merge-cali.js` command line.

[mp]: https://github.com/lh3/miniprot
[mm]: https://github.com/lh3/minimap2
[zn]: https://zenodo.org/records/15446314
