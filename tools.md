---
layout: workshop_main_3day
permalink: /EPI_2023
title: Epigenomic Data Analysis
header1: Workshop Pages for Students
header2: Epigenomic Analysis 2023
image: /site_images/CBW_Epigenome-data_icon.jpg
length: 3 days
---
## conda
Install Miniconda by following the instructoion at [Miniconda official site](https://docs.conda.io/en/main/miniconda.html)
## bedtools
```
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
mv bedtools.static.binary bedtools
chmod +x bedtools
```
## bwa
```
wget https://newcontinuum.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
tar -jxvf bwa-0.7.17.tar.bz2
cd bwa-0.7.17/
make
```
## FastQC
```
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
```
## samtools
```
wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2
tar -jxvf samtools-1.17.tar.bz2
cd samtools-1.17
make
sudo make install
```
## picard
```
wget https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar
```
## deeptools
```
pip install deeptools
```
## MACS2
```
pip install MACS2
```
## Homer
```
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl ./configureHomer.pl -install
perl ./configureHomer.pl -install hg19
perl ./configureHomer.pl -install hg38
```
## UCSC tools
```
rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/ ./
```
## Biamark
```
wget https://github.com/FelixKrueger/Bismark/archive/refs/tags/v0.24.2.tar.gz
tar -zxf v0.24.2.tar.gz
```
