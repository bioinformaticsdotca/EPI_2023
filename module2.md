---
output: html_document
permalink: /EPI_2023_module2
title: Epigenomics Analysis Module 2 Lab
header1: Workshop Pages for Students
header2: Epigenomic Analysis 2023
image: /site_images/CBW_Epigenome-data_icon.jpg
home: https://bioinformaticsdotca.github.io/EPI_2023
---

# Introduction to ChIP-Seq Analysis
by Martin Hirst and Edmund Su


## Table of Contents
1. Breakdown of markdown file
2. Module 2. Peak Calling
    - Dedup BAM file
    - Running Peak Caller (Narrow + Broad)
    - Running Peak call for ATAC
    - Filter Blacklist regions
    - Visualization coverage and peak files
    - Quality control via Fraction of reads in peaks and known enriched areas
3. Server side resources
## Breakdown of markdown file
- At the start of each step, the intention will declare.
- this is then follow by a code block

**Code:**
```
Like this! This is the main code to run for the step.

Additionally the code block will include a header to indicate what environment to the run code for example:
###Shell###
pwd

###R###
getwd()
```

- explaining for commands will be broken down following code block
- `pwd` & `getwd()`
    - see your work directory
- sprinkled throughout will also include comments

> [!NOTE]
> Important points and considerations will also be raised as so.

## Module 2. Peak Calling
### Step 1: Dedup BAM file
- We first remove duplicates here b/c `MACS2`(peak caller) identifies duplicates via genomic coordinates (excessive removal of a lot reads).
- dedup now and instruct our peak caller to "keep duplicates" after.

**Code:**
```
###Shell###
treatment=MCF10A_H3K27ac_chr19
treatment_bam=~/workspace/module123/alignments/${treatment}.sorted.dup_marked.bam
treatment_dedup=~/workspace/module123/alignments/${treatment}.sorted.dup_marked.dedup.bam
samtools view -@4 ${treatment_bam} -bh -q10 -F1028 -o ${treatment_dedup}

treatment=MCF10A_H3K27me3_chr19
treatment_bam=~/workspace/module123/alignments/${treatment}.sorted.dup_marked.bam
treatment_dedup=~/workspace/module123/alignments/${treatment}.sorted.dup_marked.dedup.bam
samtools view -@4 ${treatment_bam} -bh -q10 -F1028 -o ${treatment_dedup}

treatment=MCF10A_H3K4me3_chr19
treatment_bam=~/workspace/module123/alignments/${treatment}.sorted.dup_marked.bam
treatment_dedup=~/workspace/module123/alignments/${treatment}.sorted.dup_marked.dedup.bam
samtools view -@4 ${treatment_bam} -bh -q10 -F1028 -o ${treatment_dedup}

treatment=MCF10A_ATAC_chr19
treatment_bam=~/workspace/module123/alignments/${treatment}.sorted.dup_marked.bam
treatment_dedup=~/workspace/module123/alignments/${treatment}.sorted.dup_marked.dedup.bam
samtools view -@4 ${treatment_bam} -bh -q10 -F1028 -o ${treatment_dedup}

input=MCF10A_input_chr19
input_bam=~/workspace/module123/alignments/${input}.sorted.dup_marked.bam
input_dedup=~/workspace/module123/alignments/${input}.sorted.dup_marked.dedup.bam
samtools view -@4 ${input_bam} -bh -q10 -F1028 -o ${input_dedup}
```
- Run time <1 min per command 
- Silently completes; no log needed
- `samtools view -@4 ${input_bam} -bh -q10 -F1028 -o ${input_dedup}`
  - we subset our file based on the following criteria
  - `-bh` we would like the output to be in `BAM` format and contain the header
  - `-q10` reads with mapping qual `>10`
  - `-F1028` any reads that do not have the flags `UNMAP` and `DUP`
    - whats the difference between `-f` and `-F`?
### Step 2A : Run Peak Caller (narrow)
- MACS has two modes for narrow marks and broad marks.
- Refer to this [link](https://www.encodeproject.org/chip-seq/histone/#:~:text=must%20pass%20routine%20metadata%20audits%20in%20order%20to%20be%20released.-,Target%2Dspecific%20Standards,-For%20narrow%2Dpeak%20histone%20experiments%2C%20each%20replicate%C2%A0should%20have%2020) to reference which mark are narrow or broad

**Code:**
```
###Shell###
mkdir -p ~/workspace/module123/peaks
name=MCF10A_H3K27ac
treatment=~/workspace/module123/alignments/${name}_chr19.sorted.dup_marked.dedup.bam
input=~/workspace/module123/alignments/MCF10A_input_chr19.sorted.dup_marked.dedup.bam


macs2 callpeak -t ${treatment} -c ${input} -f BAMPE -g 58617616 -n ${name} --keep-dup all --outdir ~/workspace/module123/peaks/ --bdg 1> ~/workspace/module123/peaks/${name}.out.log 2> ~/workspace/module123/peaks/${name}.err.log

name=MCF10A_H3K4me3
treatment=~/workspace/module123/alignments/${name}_chr19.sorted.dup_marked.dedup.bam
input=~/workspace/module123/alignments/MCF10A_input_chr19.sorted.dup_marked.dedup.bam


macs2 callpeak -t ${treatment} -c ${input} -f BAMPE -g 58617616 -n ${name} --keep-dup all --outdir ~/workspace/module123/peaks/ --bdg 1> ~/workspace/module123/peaks/${name}.out.log 2> ~/workspace/module123/peaks/${name}.err.log
```
- `macs2 callpeak -t ${treatment} -c ${input} -f BAMPE -g 58617616 -n ${name} --keep-dup all --outdir ~/workspace/module123/peaks/ --bdg 1> ~/workspace/module123/peaks/${name}.out.log 2> ~/workspace/module123/peaks/${name}.err.log`
    - `macs2 callpeak` General purpose peak calling mode
    - `-t ${treatment}` Treatment file, can provide more than one to "pool"
    - `-c ${input}` Control File/Input
    - `-f BAMPE` Instructs MACS2 on what kind of file to expect. Single/Paired-end bed/bam
    - `-g hs` Sets appropriate genome size for background sampling. Typically would set `hs=human` `mm=mouse` but in our case we use the HG38 size of chr19 
    - `-n ${name}` name or prefix to use
    - `-q 0.05` FDR q value default
    - `--outdir ~/workspace/module123/peaks/` where to output files otherwise stores in current working directory
    - `--bdg` outputs pileup into bedgraph (a `BED` file where the fourth column is pileup/fragment count/coverage)
    - `1> ~/workspace/module123/peaks/${name}.out.log` output log
    - `2>  ~/workspace/module123/peaks/${name}.err.log` error log
    - let's inspect the peaks file
        - chromosome name
        - peak start
        - peak stop
        - peak name
        - int(-10*log10pvalue)
        - strand
        - Fold change at peak summit
        - -log10 P-value at Peak
        - -log10 Q-value at Peak
        - Summit position relative to peak

### Step 2B : Run Peak Caller (broad)
**Code:**
```
###Shell###
mkdir -p ~/workspace/module123/peaks
name=MCF10A_H3K27me3
treatment=~/workspace/module123/alignments/${name}_chr19.sorted.dup_marked.dedup.bam
input=~/workspace/module123/alignments/MCF10A_input_chr19.sorted.dup_marked.dedup.bam

macs2 callpeak -t ${treatment} -c ${input} -f BAMPE -g 58617616 -n ${name} --keep-dup all --outdir ~/workspace/module123/peaks/ --bdg --broad 1> ~/workspace/module123/peaks/${name}.out.log 2> ~/workspace/module123/peaks/${name}.err.log
```
- `macs2 callpeak -t ${treatment} -c ${input} -f BAMPE -g 58617616 -n ${name} --keep-dup all --outdir ~/workspace/module123/peaks/ --bdg --broad 1> ~/workspace/module123/peaks/${name}.out.log 2> ~/workspace/module123/peaks/${name}.err.log`
    - `--broad` for broad marks - stitches small peaks together
    - let's note the differences between `broadPeak` vs `gappedPeak`
    - `gapped` has the additional columns starting from `strand` where most of the differences are visualization support for the narrrow peaks within the broadPeak
        - thickStart 
        - thickEnd 
        - itemRgb 
        - blockCount 
        - blockSizes 
        - blockStarts 
        - signalValue 
        - pValue 
        - qValue
### Step 3A : Run Peak Caller for ATAC - Make fragment file
- Convert our `BAM` to `BED` to easily manipulate coordinates
- perform shift to account for Tn5 binding as a dimer and inserts two adapters separated by 9 bp
> [!NOTE]
> This step can also be done for chipseq

**Code:**
```
name=MCF10A_ATAC
dedup=~/workspace/module123/alignments/${name}_chr19.sorted.dup_marked.dedup.bam
nsort=~/workspace/module123/alignments/${name}_chr19.nsorted.dup_marked.dedup.bam
tags=~/workspace/module123/peaks/${name}_chr19.frags.bed
samtools sort -@ 8 -n ${dedup} -o ${nsort}
bedtools bamtobed -bedpe -mate1 -i ${nsort} > ${tags}
```

- `samtools sort -@ 8 -n ${dedup} -o ${nsort}`
    - specifying the `-n` sorts according to read name rather than coordinates by default
- `bedtools bamtobed -bedpe -mate1 -i ${nsort} > ${tags}`
    - `bedtools bamtobed` converts our `BAM` to `BED`/coordinates for our reads
    - `-bedpe` instructs bedtools to report read pairs instead of each read individually
    - unpaired reads will emit a warning
    - `-mate1` instructs bedtools to report read1 information first followed by read2
    - columns : `chrR1,startR1,endR1,chrR2,startR2,endR2,readName,mapQ,strandR1,strandR2`
    - e.g. `chr19 4534039 4534190 chr19 4534110 4534248 SRR20814384.28 60 + -`
    - bedtools will check if file is name sorted or coordinate sorted
### Step 3B : Run Peak Caller for ATAC - Tn5 Shift
- perform shift to account for Tn5 binding as a dimer and inserts two adapters separated by 9 bp
> [!NOTE]
> Tn5 shift can be skipped if one is not interested in footprinting.

**Code:**
```
###Shell###
name=MCF10A_ATAC

tags=~/workspace/module123/peaks/${name}_chr19.frags.bed
tn5_tags=~/workspace/module123/peaks/${name}_chr19.frags.tn5.bed

<!-- {% raw %} -->

cat ${tags} \
| awk 'BEGIN{{OFS="\t"}}{{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}}'\
| awk 'BEGIN{{OFS = "\t"}} {{if ($6 == "+") {{$2 = $2 + 4}} else if ($6 == "-") {{$3 = $3 - 5}} if ($2 >= $3) {{ if ($6 == "+") {{$2 = $3 - 1}} else {{$3 = $2 + 1}} }} print $0}}'\
> ${tn5_tags}
<!-- {% endraw %} -->

```
- `cat ${tags} | 
awk 'BEGIN{{OFS="\t"}}{{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}}' |
awk 'BEGIN {{OFS = "\t"}} {{if ($6 == "+") {{$2 = $2 + 4}} else if ($6 == "-") {{$3 = $3 - 5}} if ($2 >= $3) {{ if ($6 == "+") {{$2 = $3 - 1}} else {{$3 = $2 + 1}} }} print $0}}' > ${tn5_tags}`
  - `read tagFile | rearrange columns | shift Tag depending on strand`
  - see below for a more through breakdown of code
```
awk '
BEGIN{{OFS="\t"}}
{
    {
        printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",
        $1,$2,$3,$9,$4,$5,$6,$10
    }
}'
```
- `BEGIN{{OFS="\t"}}` output file as tab delimited
- `printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2...` print string where `%s` is substituted with the next specified element
  - equivalent to `print chrR1,startR1,endR1,"N","1000",strandR1,chrR2,startR2,endR2,strandR2`
```
awk '
BEGIN {{OFS = "\t"}} 
{
    {
        if ($6 == "+") 
            {{$2 = $2 + 4}} 
        else if ($6 == "-") 
            {{$3 = $3 - 5}} 
        if ($2 >= $3) 
        {
            { 
                if ($6 == "+") 
                {{$2 = $3 - 1}} 
                else 
                {{$3 = $2 + 1}} 
            }
        } print $0
    }
}' 
> ${tn5_tags}
```
- if fragment is on `+` shifted 4 bp to the right 
- if fragment is on `-` shift 5 bp to the left
    - while this is recommended, others have shifted `4/4` instead of `4/5`  
- after the fix is start coordinates is greater than

### Step 3C : Run Peak Caller for ATAC - Peak calling
**Code:**
```
###Shell###
name=MCF10A_ATAC
tn5_tags=~/workspace/module123/peaks/${name}_chr19.frags.tn5.bed

macs2 callpeak \
-t ${tn5_tags} \
-f BED \
-n ${name} \
-g 58617616 \
-p 0.01 \
--shift -100 \
--extsize 200 \
--nomodel \
--bdg \
--keep-dup all \
--outdir ~/workspace/module123/peaks/
```
- `macs2 callpeak -t /home/ubuntu/workspace/module123/peaks/MCF10A_ATAC_chr19.frags.tn5.bed -f BED -n {name} -g 58617616 -p 0.01 --nomodel --extsize 200 --shift -100 --bdg --keep-dup all `
    - `-f BED` we specify a `BED` format input instead of `BAM`
    - `-name ${name}` name of file
    - `-g 58617616` size of chr19
    - `-p 0.01` Pvalue cutoff
    - `--nomodel` off by default, normally calculates the `extsize` and `shift` parameters
    - `--extsize 200`, b/c Tn5's activity is on the `5'`, we extend the `5'` to get more representation
    - `--shift -100` should be `-1*(extsize/2)`, this focuses `MACS` on our `5'`
    - `--bdg` generate a pileup bedgraph
    - `--keep-dup all` retain "duplicates". As we've already filtered out duplicates, MACS2 will call duplicates via genomic coordinates
    - `--call-summits` 

### Step 4 : Blacklist removal
- [Problematic regions](https://www.nature.com/articles/s41598-019-45839-z) can obscure our results, thus filter any peaks that coincide with those regions.

**Code:**
```
###Shell###
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O ~/workspace/module123/resources/hg38_blacklist.bed.gz

gunzip ~/workspace/module123/resources/hg38_blacklist.bed.gz

blacklist=~/workspace/module123/resources/hg38_blacklist.bed

sample="MCF10A_H3K27ac_peaks"
bedtools intersect -v -a ~/workspace/module123/peaks/${sample}.narrowPeak -b ${blacklist} > ~/workspace/module123/peaks/${sample}.blacklistRemoved.narrowPeak
bedtools intersect -u -a ~/workspace/module123/peaks/${sample}.narrowPeak -b ${blacklist} > ~/workspace/module123/peaks/${sample}.blacklist.narrowPeak

sample="MCF10A_H3K4me3_peaks"
bedtools intersect -v -a ~/workspace/module123/peaks/${sample}.narrowPeak -b ${blacklist} > ~/workspace/module123/peaks/${sample}.blacklistRemoved.narrowPeak
bedtools intersect -u -a ~/workspace/module123/peaks/${sample}.narrowPeak -b ${blacklist} > ~/workspace/module123/peaks/${sample}.blacklist.narrowPeak

sample="MCF10A_H3K27me3_peaks"
bedtools intersect -v -a ~/workspace/module123/peaks/${sample}.broadPeak -b ${blacklist} > ~/workspace/module123/peaks/${sample}.blacklistRemoved.broadPeak
bedtools intersect -u -a ~/workspace/module123/peaks/${sample}.broadPeak -b ${blacklist} > ~/workspace/module123/peaks/${sample}.blacklist.broadPeak

sample="MCF10A_ATAC_peaks"
bedtools intersect -v -a ~/workspace/module123/peaks/${sample}.narrowPeak -b ${blacklist} > ~/workspace/module123/peaks/${sample}.blacklistRemoved.narrowPeak
bedtools intersect -u -a ~/workspace/module123/peaks/${sample}.narrowPeak -b ${blacklist} > ~/workspace/module123/peaks/${sample}.blacklist.narrowPeak
```
- `bedtools intersect -u -a ~/workspace/module123/peaks/${sample}.narrowPeak -b ${blacklist}`
    - `bedtools intersect` identify elements from `fileA` and `fileB` that overlap genomic coordiante wise
    - `-u` by default, bedtools will output every time a element intersection is detected i.e. if fileA_ele1 overlaps fileB_ele1 and fileB_ele2. The `-u` instead reports the element once regardless of how many overlaps
- `bedtools intersect -v -a ~/workspace/module123/peaks/${sample}.narrowPeak -b ${blacklist}`
    - `-v` reverse the behaviour, identify elements that do not overlap
- we'll return to this later when we can visualize the peaks

### Step 5A : Visualization of pileup tracks
- in the next steps, we convert our pipleup bedgraphs and bed peak files into a smaller managable formats.

**Code:**
```
###Shell###
mkdir -p ~/workspace/module123/{bigBed,bigWig}
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes -O ~/workspace/module123/resources/hg38.chrom.sizes

sample="MCF10A_H3K27ac"
chrom_sizes=~/workspace/module123/resources/hg38.chrom.sizes
input_bedgraph=~/workspace/module123/peaks/${sample}_treat_pileup.bdg
output_bigwig=~/workspace/module123/bigWig/${sample}_treat_pileup.bw

sort -k1,1 -k2,2n ${input_bedgraph} > ~/workspace/module123/bigWig/tmp
bedGraphToBigWig ~/workspace/module123/bigWig/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/module123/bigWig/tmp

input_bedgraph=~/workspace/module123/peaks/MCF10A_H3K27me3_control_lambda.bdg
output_bigwig=~/workspace/module123/bigWig/MCF10A_Input_control_pileup.bw

sort -k1,1 -k2,2n ${input_bedgraph} > ~/workspace/module123/bigWig/tmp
bedGraphToBigWig ~/workspace/module123/bigWig/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/module123/bigWig/tmp
```
- `sort -k1,1 -k2,2n ${input_bedgraph}` sorting the BED file
    - `-k1,1` sort first by chromosome alphabetically
    - `-k2,2n` sort secondarily by genomic coordinates
- `bedGraphToBigWig` convert bedGraph file to bigWig
### Step 5B : Visualization of pileup tracks continued
- we're converting the `bedgraph` file to a `bigWig` file.
- note `bedgraph` is an extension of bed(genomic coordinates) + a column with a numeric value
    - in our case pileup/coverage
- a bigWig is a binary version of a `wig` file
    - wig has a different format : https://useast.ensembl.org/info/website/upload/wig.html
    - the preferred convention for displaying data on a track  
  
**Code:**  

```
###Shell###
sample="MCF10A_H3K27me3"
chrom_sizes=~/workspace/module123/resources/hg38.chrom.sizes
input_bedgraph=~/workspace/module123/peaks/${sample}_treat_pileup.bdg
output_bigwig=~/workspace/module123/bigWig/${sample}_treat_pileup.bw

sort -k1,1 -k2,2n ${input_bedgraph} > ~/workspace/module123/bigWig/tmp
bedGraphToBigWig ~/workspace/module123/bigWig/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/module123/bigWig/tmp

sample="MCF10A_H3K4me3"
chrom_sizes=~/workspace/module123/resources/hg38.chrom.sizes
input_bedgraph=~/workspace/module123/peaks/${sample}_treat_pileup.bdg
output_bigwig=~/workspace/module123/bigWig/${sample}_treat_pileup.bw

sort -k1,1 -k2,2n ${input_bedgraph} > ~/workspace/module123/bigWig/tmp
bedGraphToBigWig ~/workspace/module123/bigWig/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/module123/bigWig/tmp

sample="MCF10A_ATAC"
chrom_sizes=~/workspace/module123/resources/hg38.chrom.sizes
input_bedgraph=~/workspace/module123/peaks/${sample}_treat_pileup.bdg
output_bigwig=~/workspace/module123/bigWig/${sample}_treat_pileup.bw

sort -k1,1 -k2,2n ${input_bedgraph} > ~/workspace/module123/bigWig/tmp
bedGraphToBigWig ~/workspace/module123/bigWig/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/module123/bigWig/tmp
```
### Step 5C : Visualization of peak tracks
- next we convert our `BED` files
**Code:**
```
###Shell###
mkdir -p ~/workspace/module123/{bigBed,bigWig}
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes -O ~/workspace/module123/resources/hg38.chrom.sizes

sample="MCF10A_H3K27ac"
chrom_sizes=~/workspace/module123/resources/hg38.chrom.sizes
input_bed=~/workspace/module123/peaks/${sample}_peaks.blacklistRemoved.narrowPeak
output_bigbed=~/workspace/module123/bigBed/${sample}.blacklistRemoved.bb


sort -k1,1 -k2,2n ${input_bed} | cut -f1-3 > ~/workspace/module123/bigBed/tmp
bedToBigBed ~/workspace/module123/bigBed/tmp ${chrom_sizes} ${output_bigbed}
rm ~/workspace/module123/bigBed/tmp
```
- `sort -k1,1 -k2,2n ${input_bedgraph}` sorting the BED file
    - `-k1,1` sort first by chromosome alphabetically
    - `-k2,2n` sort secondarily by genomic coordinates
- `bedGraphToBigWig` convert bedGraph file to bigBed
### Step 5D : Visualization of peak tracks continued

**Code:**
```
###Shell###
sample="MCF10A_ATAC"
chrom_sizes=~/workspace/module123/resources/hg38.chrom.sizes
input_bed=~/workspace/module123/peaks/${sample}_peaks.blacklistRemoved.narrowPeak
output_bigbed=~/workspace/module123/bigBed/${sample}.blacklistRemoved.bb


sort -k1,1 -k2,2n ${input_bed} | cut -f1-3 > ~/workspace/module123/bigBed/tmp
bedToBigBed ~/workspace/module123/bigBed/tmp ${chrom_sizes} ${output_bigbed}
rm ~/workspace/module123/bigBed/tmp

sample="MCF10A_H3K4me3"
chrom_sizes=~/workspace/module123/resources/hg38.chrom.sizes
input_bed=~/workspace/module123/peaks/${sample}_peaks.blacklistRemoved.narrowPeak
output_bigbed=~/workspace/module123/bigBed/${sample}.blacklistRemoved.bb


sort -k1,1 -k2,2n ${input_bed} | cut -f1-3 > ~/workspace/module123/bigBed/tmp
bedToBigBed ~/workspace/module123/bigBed/tmp ${chrom_sizes} ${output_bigbed}
rm ~/workspace/module123/bigBed/tmp

sample="MCF10A_H3K27me3"
chrom_sizes=~/workspace/module123/resources/hg38.chrom.sizes
input_bed=~/workspace/module123/peaks/${sample}_peaks.blacklistRemoved.broadPeak
output_bigbed=~/workspace/module123/bigBed/${sample}.blacklistRemoved.bb


sort -k1,1 -k2,2n ${input_bed} | cut -f1-3 > ~/workspace/module123/bigBed/tmp
bedToBigBed ~/workspace/module123/bigBed/tmp ${chrom_sizes} ${output_bigbed}
rm ~/workspace/module123/bigBed/tmp
```
### Step 5E : Visualization of peaks and tracks
- can either download your tracks and load them from files or via URL
- colours used:
  - H3K27ac (`~/workspace/module123/bigWig/MCF10A_H3K27ac_treat_pileup.bw`) : blue
  - H3K27me3 (`~/workspace/module123/bigWig/MCF10A_H3K27me3_treat_pileup.bw`): Brown
  - H3K4me3 (`~/workspace/module123/bigWig/MCF10A_H3K4me3_treat_pileup.bw`) : Green
  - ATAC (`~/workspace/module123/bigWig/MCF10A_ATAC_treat_pileup.bw`): LightBlue
  - Input/Control (`~/workspace/module123/bigWig/MCF10A_Input_control_pileup.bw`): Black
    
<img src="https://github.com/bioinformaticsdotca/EPI_2023/blob/module123/module123_images/BCL3.png?raw=true" alt="Region" width="750" />

- BCL3 start site is enriched with H3K4me3 and H3K27ac and accessible
  
<img src="https://github.com/bioinformaticsdotca/EPI_2023/blob/module123/module123_images/RDH8.png?raw=true" alt="Region" width="750" />

- Genes around RDH8 show enrichment of H3K27me3 and lack of accessibility and enrichment for H3K27ac and H3K4me3
  
<img src="https://github.com/bioinformaticsdotca/EPI_2023/blob/module123/module123_images/blacklist.png?raw=true" alt="Region" width="750" />

- Blacklist identified region hightlighted in red, represented by both H3K27me3 and ATAC

### Step 6A : Quality Control (Enrichment in key genomic areas)
- A way to determine the efficacy of your enrichment is to benchmark % of reads to called peaks (FRIP) or known regions
- in our toy example we'll be examining enrichment at promoters (TSS+/- 2kb) and encode defining enhancer regions
- For H3K27me3 we'd typically look at HOX regions however no HOX regions on chr19
**Code :** 
```
###Shell###
mkdir ~/workspace/module123/qc
wget https://www.bcgsc.ca/downloads/esu/touchdown/hg38v79_genes_tss_2000.bed -O ~/workspace/module123/resources/hg38v79_genes_tss_2000.bed

sort -k1,1 -k2,2n ~/workspace/module123/resources/hg38v79_genes_tss_2000.bed > tmp
mv tmp ~/workspace/module123/resources/hg38v79_genes_tss_2000.bed

wget https://www.bcgsc.ca/downloads/esu/touchdown/encode_enhancers_liftover.bed -O ~/workspace/module123/resources/encode_enhancers_liftover.bed

TSS=~/workspace/module123/resources/hg38v79_genes_tss_2000.bed
ENH=~/workspace/module123/resources/encode_enhancers_liftover.bed

sample="MCF10A_H3K27ac"
query_bam=~/workspace/module123/alignments/${sample}_chr19.sorted.dup_marked.bam

samtools view -@4 -q 10 -F 1028 $query_bam -c
samtools view -@4 -q 10 -F 1028 $query_bam -L ${ENH} -c
samtools view -@4 -q 10 -F 1028 $query_bam -L ~/workspace/module123/peaks/${sample}_peaks.blacklistRemoved.narrowPeak -c
```
- `samtools view -@4 -q 10 -F 1028 $query_bam -L ~/workspace/module123/peaks/${sample}_peaks.blacklistRemoved.narrowPeak -c`
    - `samtools view` parse through our BAM file
    - `-@4` resources to use
    - `-q 10` filter for reads with mapping quality(MAPQ) >10
    - `-F 1028` Remove reads that have the following flags:UNMAP,DUP
    - `$query_bam` Bam of interest
    - `-c ` count the number of reads that fulfil the criteria
    - `-L` perform actions on reads that fall within the specified genomic coordiantes (Bed File)
### Step 6B : Quality Control (Enrichment in key genomic areas) Continued:

**Code:**
```
###Shell###
TSS=~/workspace/module123/resources/hg38v79_genes_tss_2000.bed
ENH=~/workspace/module123/resources/encode_enhancers_liftover.bed

for histone in H3K27ac H3K27me3 H3K4me3 ATAC input;
do
    sample="MCF10A_${histone}"
    query_bam=~/workspace/module123/alignments/${sample}_chr19.sorted.dup_marked.bam
    samtools view -@4 -q 10 -F 1028 $query_bam -c > ~/workspace/module123/qc/col_${histone}
    samtools view -@4 -q 10 -F 1028 $query_bam -L ${TSS} -c >> ~/workspace/module123/qc/col_${histone}
    samtools view -@4 -q 10 -F 1028 $query_bam -L ${ENH} -c >> ~/workspace/module123/qc/col_${histone}

    if [[ "$histone" == "H3K27me3" ]]; then
        peaks=~/workspace/module123/peaks/${sample}_peaks.blacklistRemoved.broadPeak
        samtools view -@4 -q 10 -F 1028 $query_bam -L ~/workspace/module123/peaks/${sample}_peaks.blacklistRemoved.broadPeak -c >> ~/workspace/module123/qc/col_${histone}
    elif [[ "$histone" == "H3K27ac" || "$histone" == "H3K4me3" || "$histone" == "ATAC" ]]; then
        peaks=~/workspace/module123/peaks/${sample}_peaks.blacklistRemoved.narrowPeak
        samtools view -@4 -q 10 -F 1028 $query_bam -L ~/workspace/module123/peaks/${sample}_peaks.blacklistRemoved.narrowPeak -c >> ~/workspace/module123/qc/col_${histone}
    else 
        echo
    fi
done

paste ~/workspace/module123/qc/col_H3K27ac ~/workspace/module123/qc/col_H3K27me3 ~/workspace/module123/qc/col_H3K4me3 ~/workspace/module123/qc/col_ATAC ~/workspace/module123/qc/col_input > ~/workspace/module123/qc/integers.tsv
paste \
<(awk '{print $1/442829*100}' ~/workspace/module123/qc/col_H3K27ac) \
<(awk '{print $1/1610910*100}' ~/workspace/module123/qc/col_H3K27me3) \
<(awk '{print $1/1512352*100}' ~/workspace/module123/qc/col_H3K4me3) \
<(awk '{print $1/2751833*100}' ~/workspace/module123/qc/col_ATAC) \
<(awk '{print $1/1023265*100}' ~/workspace/module123/qc/col_input) > ~/workspace/module123/qc/percentages.tsv
```
`paste ~/workspace/module123/qc/col_H3K27ac ~/workspace/module123/qc/col_H3K27me3 ~/workspace/module123/qc/col_H3K4me3 ~/workspace/module123/qc/col_ATAC ~/workspace/module123/qc/col_input > ~/workspace/module123/qc/integers.tsv`
  - `paste` lets us aggregate our results where each is a column)
- 
- Reads enrichment in key region 

| |H3K27ac|H3K27me3|H3K4me3|ATAC|Input|
|------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
|Total|442829|1610910|1512352|2751833|1023265
|TSS|127577|298326|1087032|924326|194996
|Enhancer|242489|760089|834053|1721794|514115
|In Peaks|69584|783294|1239717|1082674|


 	 	 	 	 
- Reads % enrichment in key region
				
| |H3K27ac|H3K27me3|H3K4me3|ATAC|Input|
|------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
|TSS|28.8095|18.5191|71.8769|33.5895|19.0563
|Enhancer|54.7591|47.1838|55.1494|62.569|50.2426
|In Peaks|15.7135|48.6243|81.9728|39.3437|

- `H3K27ac` FRIP is poor as a result of the poorer sequencing depth.
- `H3K4me3` is performing well: high FRIP and enriched in promoter regions

## Server resources
### QC Resources
```
###TSS+/-2kb
mkdir ~/workspace/module123/qc
wget https://www.bcgsc.ca/downloads/esu/touchdown/hg38v79_genes_tss_2000.bed -O ~/workspace/module123/resources/hg38v79_genes_tss_2000.bed

sort -k1,1 -k2,2n ~/workspace/module123/resources/hg38v79_genes_tss_2000.bed > tmp
mv tmp ~/workspace/module123/resources/hg38v79_genes_tss_2000.bed

###Enhancer liftover
wget https://www.bcgsc.ca/downloads/esu/touchdown/encode_enhancers_liftover.bed -O ~/workspace/module123/resources/encode_enhancers_liftover.bed

###Blacklist
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O ~/workspace/module123/resources/hg38_blacklist.bed.gz

gunzip ~/workspace/module123/resources/hg38_blacklist.bed.gz
```
- `hg38v79_genes_tss_2000.bed`
    - Generated by downloading Ensemblv79 GTF convert to Bed +/-2kb of TSS. See https://www.biostars.org/p/56280/
- `encode_enhancers_liftover.bed`
    - download various [ChroHMM state7](https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/) and merge


### Encode Bed
```
ls ~/CourseData/EPI_data/module123/encode_bed
basal.H3K27ac.peak_calls.bed
basal.H3K27me3.peak_calls.bed
basal.H3K4me1.peak_calls.bed
basal.H3K4me3.peak_calls.bed
lp.H3K27ac.peak_calls.bed
lp.H3K27me3.peak_calls.bed
lp.H3K4me1.peak_calls.bed
lp.H3K4me3.peak_calls.bed
luminal.H3K27ac.peak_calls.bed
luminal.H3K27me3.peak_calls.bed
luminal.H3K4me1.peak_calls.bed
luminal.H3K4me3.peak_calls.bed
stromal.H3K27ac.peak_calls.bed
stromal.H3K27me3.peak_calls.bed
stromal.H3K4me1.peak_calls.bed
stromal.H3K4me3.peak_calls.bed
```
- https://epigenomesportal.ca/tracks/CEEHRC/hg38/
- Breast Basal CEMT0035
- Breast Stromal CEMT0036
- Breast Luminal CEMT0037
- Breast Luminal Progenitor CEMT0038
### Encode BigWig
```
ls ~/CourseData/EPI_data/module123/encode_bigWig
basal.H3K27ac.signal_unstranded.bigWig
basal.H3K27me3.signal_unstranded.bigWig
basal.H3K4me1.signal_unstranded.bigWig
basal.H3K4me3.signal_unstranded.bigWig
lp.H3K27ac.signal_unstranded.bigWig
lp.H3K27me3.signal_unstranded.bigWig
lp.H3K4me1.signal_unstranded.bigWig
lp.H3K4me3.signal_unstranded.bigWig
luminal.H3K27ac.signal_unstranded.bigWig
luminal.H3K27me3.signal_unstranded.bigWig
luminal.H3K4me1.signal_unstranded.bigWig
luminal.H3K4me3.signal_unstranded.bigWig
stromal.H3K27ac.signal_unstranded.bigWig
stromal.H3K27me3.signal_unstranded.bigWig
stromal.H3K4me1.signal_unstranded.bigWig
stromal.H3K4me3.signal_unstranded.bigWig
```
- https://epigenomesportal.ca/tracks/CEEHRC/hg38/
- Breast Basal CEMT0035
- Breast Stromal CEMT0036
- Breast Luminal CEMT0037
- Breast Luminal Progenitor CEMT0038
### MCF10A Fastq
```
ls ~/CourseData/EPI_data/module123/fastq
MCF10A.ATAC.chr19.R1.fastq.gz
MCF10A.ATAC.chr19.R2.fastq.gz
MCF10A.H3K27ac.chr19.R1.fastq.gz
MCF10A.H3K27ac.chr19.R2.fastq.gz
MCF10A.H3K27me3.chr19.R1.fastq.gz
MCF10A.H3K27me3.chr19.R2.fastq.gz
MCF10A.H3K4me3.chr19.R1.fastq.gz
MCF10A.H3K4me3.chr19.R2.fastq.gz
MCF10A.Input.chr19.R1.fastq.gz
MCF10A.Input.chr19.R2.fastq.gz
```
- MCF10A histone marks and input come courtesy of Dr.Hirst
- ATACseq data originates from [GSM6431322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6431322)
### Triplicates
```
CourseData/EPI_data/module123/triplicates/triplicates.csv

CourseData/EPI_data/module123/triplicates/alignments:
MCF10A_H3K4me3_chr19.CondA.Rep1.bam      MCF10A_H3K4me3_chr19.CondB.Rep2.bam      MCF10A_input_chr19.CondA.Rep3.bam
MCF10A_H3K4me3_chr19.CondA.Rep1.bam.bai  MCF10A_H3K4me3_chr19.CondB.Rep2.bam.bai  MCF10A_input_chr19.CondA.Rep3.bam.bai
MCF10A_H3K4me3_chr19.CondA.Rep2.bam      MCF10A_H3K4me3_chr19.CondB.Rep3.bam      MCF10A_input_chr19.CondB.Rep1.bam
MCF10A_H3K4me3_chr19.CondA.Rep2.bam.bai  MCF10A_H3K4me3_chr19.CondB.Rep3.bam.bai  MCF10A_input_chr19.CondB.Rep1.bam.bai
MCF10A_H3K4me3_chr19.CondA.Rep3.bam      MCF10A_input_chr19.CondA.Rep1.bam        MCF10A_input_chr19.CondB.Rep2.bam
MCF10A_H3K4me3_chr19.CondA.Rep3.bam.bai  MCF10A_input_chr19.CondA.Rep1.bam.bai    MCF10A_input_chr19.CondB.Rep2.bam.bai
MCF10A_H3K4me3_chr19.CondB.Rep1.bam      MCF10A_input_chr19.CondA.Rep2.bam        MCF10A_input_chr19.CondB.Rep3.bam
MCF10A_H3K4me3_chr19.CondB.Rep1.bam.bai  MCF10A_input_chr19.CondA.Rep2.bam.bai    MCF10A_input_chr19.CondB.Rep3.bam.bai

CourseData/EPI_data/module123/triplicates/bigWig:
CondA.Rep1.bw  CondA.Rep2.bw  CondA.Rep3.bw  CondB.Rep1.bw  CondB.Rep2.bw  CondB.Rep3.bw

CourseData/EPI_data/module123/triplicates/peaks:
CondA.Rep1_peaks.narrowPeak  CondA.Rep3_peaks.narrowPeak  CondB.Rep2_peaks.narrowPeak
CondA.Rep2_peaks.narrowPeak  CondB.Rep1_peaks.narrowPeak  CondB.Rep3_peaks.narrowPeak
```
- triplicates were generated from MCF10A_H3K4me3 by choosing a list of exclusive peaks for condA and condB and randomly subsampling replicates accordingly
