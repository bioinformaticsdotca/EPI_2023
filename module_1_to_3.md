# Introduction to ChIP-Seq Analysis
by Martin Hirst and Edmund Su


## Table of Contents
1. Breakdown of markdown file
2. Module 1. Alignment
    - FastQC and Inspection
    - Reference Indexing
    - Alignment
    - Coordinate Sort
    - MarkDuplicates
    - Flagstats and Stats
    - Clean up
3. Module 2. Peak Calling
    - Dedup BAM file
    - Running Peak Caller
    - Filter Blacklist regions
    - Visualization
4. Module 3. Differential Analysis
    - Bedtools
    - Deeptools
    - Deeptools
    - edgeR
    - DiffBind
5. Server side resources

## Module 1. Alignment
### Step 1: Setup
- It's always good practice to organize your directories beforehand
- This avoids file conflicts and accidental data loss.
- We'll be making a subdirectory for each major step.
<Br>
**Code:**
```
mkdir -p ~/workspace/BWA_index
mkdir -p ~/workspace/alignments
mkdir -p ~/workspace/fastqc
mkdir -p ~/workspace/stats
mkdir -p ~/workspace/peaks
mkdir -p ~/workspace/bigBed
mkdir -p ~/workspace/bigWig
mkdir -p ~/workspace/resources
mkdir -p ~/workspace/diffBind
```
> [!NOTE]
> Note we can also simplfy this command.

### Step 2A: Retrieve a reference genome
- We'll need a reference genome to align our seqeuences to.
- A great resource is [UCSC's genome repository](https://hgdownload.soe.ucsc.edu/downloads.html) containing latest Hg38 and older hg19 human genome.
- For our tutorial, we'll be working `chr19`. This is to speed demonstration purposes.

**Code:** 
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr19.fa.gz -O ~/workspace/BWA_index/chr19.fa.gz
gunzip ~/workspace/BWA_index/chr19.fa.gz
```
- Runtime : <1 minute

- `wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr19.fa.gz -O ~/workspace/BWA_index/chr19.fa.gz`

    - `Wget` is a command to pull online files.
    - We're instructing it to pull `chr19.fa.gz` from the example address.
    - `-O ~/workspace/BWA_index/chr19.fa.gz` instructs where and what we would like to save our file as.
    - _What were to happen if we ran the command without `-O`?_
    - Not all systems carry/support `wget`, `Curl` is another useful alternative.

- `gunzip ~/workspace/BWA_index/chr19.fa.gz -c > ~/workspace/BWA_index/chr19.fa`
    - `gzip` is a data compression format useful for reducing storage impacts
    - `gunzip` is the opposite decompresses `gzipped` data
    - `-c` is a STDOUT redirect. Normallt running `gunzip file.gz` will turn `file.tsv.gz` into `file.tsv`. Running `-c` allows us to save the contents elsewhere. We do this we can compare the size difference.
    - compare `ls ~/workspace/BWA_index/chr19.fa.gz -ilth` vs `ls ~/workspace/BWA_index/chr19.fa -ilth` and note the size difference

> [!NOTE]
> Prior to any tool usage, it's good practice to explore the tool by invoking the help command. This may vary depending on tool. This is important as the same `-O` flag could have differing functions in different tools.<Br>
### Step 2B: Index the genome

- Indexing step creates a database to enable quick retrieve of sequences
- Indexing is only required once per genome.
- That being said, because there are different versions of Hg38/GRCh38 (such as those with or [without ALT contigs and EBV](https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19)) each version would need their own indices.
- Additionally index files are typically program/tool specific.
- For the lab, we're utilizing [BWA](https://bio-bwa.sourceforge.net/bwa.shtml) for alignment.

> [!NOTE]
> [A newer version of BWA-mem](https://github.com/bwa-mem2/bwa-mem2) does exist, however we'll be using the older version due to a bug that requires significant memory for index.

> [!NOTE]
> Consortias will also have available the reference genome and related resources. E.g. https://ewels.github.io/AWS-iGenomes/
```
mkdir -p ~/workspace/BWA_index
bwa index ~/workspace/BWA_index/chr19.fa > ~/workspace/BWA_index/index.log
```
- Runtime : <1 minute
- take a look inside your folder `ls ~/CourseData/BWA_index/`
### Step 3: FastQC and interpretation
- Checking the quality of your FASTQs is sanity check step (garbage in garbage out)
**Code:** 
```
###Shell###
fastq_input=~/CourseData/fastq/MCF10A.Input.chr19.R1.fastq.gz
fastq_h3k27ac=~/CourseData/fastq/MCF10A.H3K27ac.chr19.R1.fastq.gz
mkdir -p ~/workspace/fastqc
fastqc -t 8 ${fastq_input} -o ~/workspace/fastqc > ~/workspace/fastqc/input.log 2>&1
fastqc -t 8 ${fastq_h3k27ac} -o ~/workspace/fastqc > ~/workspace/fastqc/h3k27ac.log 2>&1
```
- `fastqc`
    - `-t 8` specifies the number of threads we want to program to utilize
    - `-o` specifies our output directory
- Investigating the files!
- TBA Add discussion of individual files. Highlight H3k27ac vs input
### Step 4: Alignment
- Mapping the read pairs to a position of the reference genome
**Code:** 
```
###Shell###
ref=~/workspace/BWA_index/chr19.fa
read1=~/CourseData/fastq/MCF10A.Input.chr19.R1.fastq.gz
read2=~/CourseData/fastq/MCF10A.Input.chr19.R2.fastq.gz
sample=MCF10A_input_chr19
bwa mem -M -t 4 ${ref} ${read1} ${read2} 2>~/workspace/alignments/${sample}.alignment.log | samtools view -hbS -o ~/workspace/alignments/${sample}.bam
```
- run time ~2 min
- The command can be broken down into the following pseudo code `ALIGNMENT | SAMtoBAM_conversion`. The `|` operate streams the results from the first command into the second.

- `bwa mem -M -t 4 ${ref} ${read1} ${read2} 2>~/workspace/alignments/alignment.log`
    - `bwa mem` while BWA has many alignment algorithms, `mem` is best suited for efficiently handling >70bp paired end reads.
    - `-M` alters flagging of chimeric reads. By default chimeric reads are flagged as `SUPPLEMENTARY` (partially mapping), the option instead turns them to `SECONDARY` (multi-mapping). Needed to GATK/PICARD support downstream
    - `-t 4` resource specification
    - `2>~/workspace/alignments/alignment.log` data outputted occurs in two streams `1` (the data) and `2` (debug messages, warnings and status updates). We redirect `2` to a log file.
    - What do we see in the logs?
> [!NOTE]
> Logs contain valuable information, helpful for debugging.
- `samtools view -hbS -o ~/workspace/alignments/${sample}.bam`
    - `samtools view` a tool to read our `SAM`,`BAM`,`CRAM` files
    - `-hbS` include header, output as `BAM`, input is `SAM`
    - `-o` specify what we want to save our file as
    - take a look our new `BAM` file. Note the 
- What other important functions of `samtools`?
    - CRAM conversion
    - samtools quickcheck
    - samtools view header
    - samtools view subset
### Step 5: Coordinate Sort
- Rearrange our alignments by coordinates

**Code:**
```
###Shell###
sample=MCF10A_input_chr19
samtools sort -@8 ~/workspace/alignments/${sample}.bam -o ~/workspace/alignments/${sample}.sorted.bam
```
- run time <1 min
- `samtools sort` sorts our reads by genome coordinates
- Observe the files before and after via `samtools view | head`
- How to sort by readname?
### Step 6: Duplicate Marking
- Identify and tag alignments that are duplicates

**Code:**
```
###Shell###
sample=MCF10A_input_chr19
java -jar /usr/local/picard/picard.jar MarkDuplicates I=~/workspace/alignments/${sample}.sorted.bam O=~/workspace/alignments/${sample}.sorted.dup_marked.bam M=~/workspace/alignments/${sample}.dup_marked.output.log ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT > ~/workspace/alignments/${sample}.dup_marked.error.log
``` 
- `java -jar /usr/local/picard/picard.jar MarkDuplicates I=~/workspace/alignments/${sample}.sorted.bam O=~/workspace/alignments/${sample}.sorted.dup_marked.bam M=~/workspace/alignments/${sample}.dup_marked.output.log ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT`
    - `java -jar /usr/local/picard/picard.jar` produced by BroadInstitute, [Picard tools are a toolset](https://broadinstitute.github.io/picard/) for HTS/NGS data
    - `MarkDuplicates` identifies reads/clusters duplicates arising from library construciton during PCR and cluster formation during sequencing
    - `I=` and `O=` are input and output respectively
    - `M=` saves our output log
    - `ASSUME_SORTED=TRUE` informs `PICARD`, our input is already coordinate sorted
    - `VALIDATION_STRINGENCY=LENIENT` informs `PICARD` to continue and notify problems in the `work.log` instead of failing

### Step 7: Stats
- Identify and tag alignments that are duplicates
**Code:**
```
###Shell###
sample=MCF10A_input_chr19
mkdir -p  ~/workspace/stats
samtools flagstat ~/workspace/alignments/${sample}.sorted.dup_marked.bam > ~/workspace/stats/${sample}.sorted.dup_marked.flagstat
samtools stats ~/workspace/alignments/${sample}.sorted.dup_marked.bam > ~/workspace/stats/${sample}.sorted.dup_marked.flagstat

```
### Step 8: Processing the remaining files
**Code:**
```
###Shell##
ref=~/workspace/BWA_index/chr19.fa
for histone in H3K27ac H3K27me3 H3K4me3 ATAC;
    do
    read1=CourseData/fastq/MCF10A.${histone}.chr19.R1.fastq.gz
    read2=CourseData/fastq/MCF10A.${histone}.chr19.R2.fastq.gz
    sample=MCF10A_${histone}_chr19
    echo "aligning ${histone}"
    bwa mem -M -t 4 ${ref} ${read1} ${read2} 2>workspace/alignments/${sample}.alignment.log | samtools view -hbS -o workspace/alignments/${sample}.bam
    echo "sorting ${histone}"
    samtools sort -@8 workspace/alignments/${sample}.bam -o workspace/alignments/${sample}.sorted.bam
    echo "dupmarking ${histone}"
    picard MarkDuplicates I=workspace/alignments/${sample}.sorted.bam O=workspace/alignments/${sample}.sorted.dup_marked.bam M=workspace/alignments/${sample}.dup_marked.output.log ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT > workspace/alignments/${sample}.dup_marked.error.log
    echo "calculating stats ${histone}"
    samtools flagstat workspace/alignments/${sample}.sorted.dup_marked.bam > workspace/stats/${sample}.sorted.dup_marked.flagstat
    samtools stats workspace/alignments/${sample}.sorted.dup_marked.bam > workspace/stats/${sample}.sorted.dup_marked.stat
    done
```

### Step 9: Clean up!
- Remove temporary files

**Code:**
```
###Shell###
ls ~/workspace/alignments/*.bam | grep -v sorted.dup_marked | xargs -I {} sh -c "echo rm {}; rm {}"
```
- `|` the pipe delimiter meaning we want to do `ACTION A | ACTION B | ACTION C` (like an assembly line)
- `ls ~/workspace/alignments/*.bam`
    - a look up of `BAM` files in our directory : `~/workspace/alignments`
    - The `*` in `*.bam` acts as a wild card for any string of variable length ending in the suffix `bam`
    - running the command on it's own produces
- `grep -v sorted.dup_marked`
    - `grep` a powerful tool that searches and matches text
    - `-v` return matches that do not much our pattern or regex
    - `sorted.dup_marked` is our pattern, returning 2/3 of the files : `MCF10A_input_chr19.bam`,`MCF10A_input_chr19.sorted.bam`,`MCF10A_input_chr19.sorted.dup_marked.bam`
- `xargs -I {} sh -c "echo rm {}; rm {}"`
    - `xargs` receives input and performs an action. Think a for-loop
    - `-I {}` the variable we want to use
    - `sh -c` interpret the string provided in bash shell
    - `echo rm {}` echos the command we want to perform
    - `rm {}` removes the file

## Module 2. Peak Calling
### Step 1: Dedup BAM file
**Code:**
```
###Shell###
treatment=MCF10A_H3K27ac_chr19
treatment_bam=~/workspace/alignments/${treatment}.sorted.dup_marked.bam
treatment_dedup=~/workspace/alignments/${treatment}.sorted.dup_marked.dedup.bam
samtools view -@4 ${treatment_bam} -bh -q10 -F1028 -o ${treatment_dedup}

treatment=MCF10A_H3K27me3_chr19
treatment_bam=~/workspace/alignments/${treatment}.sorted.dup_marked.bam
treatment_dedup=~/workspace/alignments/${treatment}.sorted.dup_marked.dedup.bam
samtools view -@4 ${treatment_bam} -bh -q10 -F1028 -o ${treatment_dedup}

treatment=MCF10A_H3K4me3_chr19
treatment_bam=~/workspace/alignments/${treatment}.sorted.dup_marked.bam
treatment_dedup=~/workspace/alignments/${treatment}.sorted.dup_marked.dedup.bam
samtools view -@4 ${treatment_bam} -bh -q10 -F1028 -o ${treatment_dedup}

input=MCF10A_input_chr19
input_bam=~/workspace/alignments/${input}.sorted.dup_marked.bam
input_dedup=~/workspace/alignments/${input}.sorted.dup_marked.dedup.bam
samtools view -@4 ${input_bam} -bh -q10 -F1028 -o ${input_dedup}
```
- Run time <1 min per command 
- Silently completes; no log needed
### Step 2A : Run Peak Caller (narrow)
**Code:**
```
###Shell###
mkdir -p ~/workspace/peaks
name=MCF10A_H3K27ac
treatment=~/CourseData/EPI_data/Module1/MCF10A_resources/${name}_chr19.sorted.dup_marked.dedup.bam
input=~/CourseData/EPI_data/Module1/MCF10A_resources/MCF10A_Input.dedup.bam


macs2 callpeak -t ${treatment} -c ${input} -f BAMPE -g 58617616 -n ${name} --keep-dup all --outdir ~/workspace/peaks/ --bdg 1> ~/workspace/peaks/${name}.out.log 2> ~/workspace/peaks/${name}.err.log

name=MCF10A_H3K4me3
treatment=~/CourseData/EPI_data/Module1/MCF10A_resources/${name}_chr19.sorted.dup_marked.dedup.bam
input=~/CourseData/EPI_data/Module1/MCF10A_resources/MCF10A_Input.dedup.bam


macs2 callpeak -t ${treatment} -c ${input} -f BAMPE -g 58617616 -n ${name} --keep-dup all --outdir ~/workspace/peaks/ --bdg 1> ~/workspace/peaks/${name}.out.log 2> ~/workspace/peaks/${name}.err.log
```
- `macs2 callpeak -t ${treatment} -c ${input} -f BAMPE -g 58617616 -n ${name} --keep-dup all --outdir ~/workspace/peaks/ --bdg 1> ~/workspace/peaks/${name}.out.log 2> ~/workspace/peaks/${name}.err.log`
    - `macs2 callpeak` General purpose peak calling mode
    - `-t ${treatment}` Treatment file, can provide more than one to "pool"
    - `-c ${input}` Control File/Input
    - `-f BAMPE` Instructs MACS2 on what kind of file to expect. Single/Paired-end bed/bam
    - `-g hs` Sets appropriate genome size for background sampling. Typically would set `hs=human` `mm=mouse` but in our case we use the HG38 size of chr19 
    - `-n ${name}` name or prefix to use
    - `-q 0.05` FDR q value default
    - `--outdir ~/workspace/peaks/` where to output files otherwise stores in current working directory
    - `--bdg` outputs pileup into bedgraphs
    - `1> ~/workspace/peaks/${name}.out.log` output log
    - `2>  ~/workspace/peaks/${name}.err.log` error log
    - let's inspect the peaks file
        - chromosome name
        - peak start
        - peak stop
        - peak name
        - int(-10*log10pvalue)
        - N/A
        - Fold change at peak summit
        - -log10 P-value at Peak
        - -log10 Q-value at Peak
        - Summit position relative to peak

### Step 2B : Run Peak Caller (broad)
**Code:** 
```
###Shell###
mkdir -p ~/workspace/peaks
name=MCF10A_H3K27me3
treatment=~/CourseData/EPI_data/Module1/MCF10A_resources/${name}_chr19.sorted.dup_marked.dedup.bamm
input=~/CourseData/EPI_data/Module1/MCF10A_resources/MCF10A_Input.dedup.bam

macs2 callpeak -t ${treatment} -c ${input} -f BAMPE -g 58617616 -n ${name} --keep-dup all --outdir ~/workspace/peaks/ --bdg --broad 1> ~/workspace/peaks/${name}.out.log 2> ~/workspace/peaks/${name}.err.log
```
- `macs2 callpeak -t ${treatment} -c ${input} -f BAMPE -g 58617616 -n ${name} --keep-dup all --outdir ~/workspace/peaks/ --bdg --broad 1> ~/workspace/peaks/${name}.out.log 2> ~/workspace/peaks/${name}.err.log`
    - `--broad` for broad marks - stitches small peaks together
    - let's note the differences between `broadPeak` vs `gappedPeak`
### Step 3 : Blacklist removal
```
###Shell###
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -o workspace/resources/hg38_blacklist.bed.gz

gunzip workspace/resources/hg38_blacklist.bed.gz

blacklist=~/CourseData/EPI_data/Module1/QC_resources/hg38_blacklist.bed

sample="MCF10A_H3K27ac_chr19_peaks"
bedtools intersect -v -a ~/workspace/peaks/${sample}.narrowPeak -b ${blacklist} > ~/workspace/peaks/${sample}.blacklistRemoved.narrowPeak
bedtools intersect -u -a ~/workspace/peaks/${sample}.narrowPeak -b ${blacklist} > ~/workspace/peaks/${sample}.blacklist.narrowPeak

sample="MCF10A_H3K4me3_chr19_peaks"
bedtools intersect -v -a ~/workspace/peaks/${sample}.narrowPeak -b ${blacklist} > ~/workspace/peaks/${sample}.blacklistRemoved.narrowPeak
bedtools intersect -u -a ~/workspace/peaks/${sample}.narrowPeak -b ${blacklist} > ~/workspace/peaks/${sample}.blacklist.narrowPeak

sample="MCF10A_H3K27me3_chr19_peaks"
bedtools intersect -v -a ~/workspace/peaks/${sample}.broadPeak -b ${blacklist} > ~/workspace/peaks/${sample}.blacklistRemoved.narrowPeak
bedtools intersect -u -a ~/workspace/peaks/${sample}.broadPeak -b ${blacklist} > ~/workspace/peaks/${sample}.blacklist.narrowPeak
```
- lets inspect the blacklist peaks

### Step 4A : Visualization of pileup tracks
**Code:**
```
###Shell###
mkdir -p ~/workspace/{bigBed,bigWig}
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes -O workspace/resources/hg38.chrom.sizes

sample="MCF10A_H3K27ac"
chrom_sizes=~/workspace/resources/hg38.chrom.sizes
input_bedgraph=~/workspace/peaks/${sample}_treat_pileup.bdg
output_bigwig=~/workspace/bigWig/${sample}_treat_pileup.bw

sort -k1,1 -k2,2n ${input_bedgraph} > ~/workspace/bigWig/tmp
bedGraphToBigWig ~/workspace/bigWig/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/bigWig/tmp

input_bedgraph=~/workspace/peaks/MCF10A_H3K27me3_control_lambda.bdg
output_bigwig=~/workspace/bigWig/MCF10A_Input_control_pileup.bw

sort -k1,1 -k2,2n ${input_bedgraph} > ~/workspace/bigWig/tmp
bedGraphToBigWig ~/workspace/bigWig/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/bigWig/tmp
```
- `sort -k1,1 -k2,2n ${input_bedgraph}` sorting the BED file
    - `-k1,1` sort first by chromosome alphabetically
    - `-k2,2n` sort secondarily by genomic coordinates
- `bedGraphToBigWig` convert bedGraph file to bigWig
### Step 4B : Visualization of pileup tracks continued
```
sample="MCF10A_H3K27me3"
chrom_sizes=~/workspace/resources/hg38.chrom.sizes
input_bedgraph=~/workspace/peaks/${sample}_treat_pileup.bdg
output_bigwig=~/workspace/bigWig/${sample}_treat_pileup.bw

sort -k1,1 -k2,2n ${input_bedgraph} > ~/workspace/bigWig/tmp
bedGraphToBigWig ~/workspace/bigWig/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/bigWig/tmp

sample="MCF10A_H3K4me3"
chrom_sizes=~/workspace/resources/hg38.chrom.sizes
input_bedgraph=~/workspace/peaks/${sample}_treat_pileup.bdg
output_bigwig=~/workspace/bigWig/${sample}_treat_pileup.bw

sort -k1,1 -k2,2n ${input_bedgraph} > ~/workspace/bigWig/tmp
bedGraphToBigWig ~/workspace/bigWig/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/bigWig/tmp
```
### Step 4C : Visualization of peak tracks
**Code:**
```
mkdir -p ~/workspace/{bigBed,bigWig}
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes -O workspace/resources/hg38.chrom.sizes

sample="MCF10A_H3K27ac"
chrom_sizes=~/workspace/resources/hg38.chrom.sizes
input_bed=~/workspace/peaks/${sample}_peaks.blacklistRemoved.narrowPeak
output_bigwig=~/workspace/bigBed/${sample}.blacklistRemoved.bb


sort -k1,1 -k2,2n ${input_bed} | cut -f1-3 > ~/workspace/bigBed/tmp
bedToBigBed ~/workspace/bigBed/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/bigWig/tmp
```
- `sort -k1,1 -k2,2n ${input_bedgraph}` sorting the BED file
    - `-k1,1` sort first by chromosome alphabetically
    - `-k2,2n` sort secondarily by genomic coordinates
- `bedGraphToBigWig` convert bedGraph file to bigWig
### Step 4D : Visualization of peak tracks continued
```
sample="MCF10A_H3K4me3"
chrom_sizes=~/workspace/resources/hg38.chrom.sizes
input_bed=~/workspace/peaks/${sample}_peaks.blacklistRemoved.narrowPeak
output_bigwig=~/workspace/bigBed/${sample}.blacklistRemoved.bb


sort -k1,1 -k2,2n ${input_bed} | cut -f1-3 > ~/workspace/bigBed/tmp
bedToBigBed ~/workspace/bigBed/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/bigWig/tmp

sample="MCF10A_H3K27me3"
chrom_sizes=~/workspace/resources/hg38.chrom.sizes
input_bed=~/workspace/peaks/${sample}_peaks.blacklistRemoved.braodPeak
output_bigwig=~/workspace/bigBed/${sample}.blacklistRemoved.bb


sort -k1,1 -k2,2n ${input_bed} | cut -f1-3 > ~/workspace/bigBed/tmp
bedToBigBed ~/workspace/bigBed/tmp ${chrom_sizes} ${output_bigwig}
rm ~/workspace/bigWig/tmp
```
### Step 4E : Visualization of peaks and tracks


### Step 5A : Quality Control (Enrichment in key genomic areas)
**Code :** 
```
###Shell###
mkdir workspace/qc
wget https://www.bcgsc.ca/downloads/esu/touchdown/hg38v79_genes_tss_2000.bed -O workspace/resources/hg38v79_genes_tss_2000.bed
wget https://www.bcgsc.ca/downloads/esu/touchdown/encode_enhancers_liftover.bed -O workspace/resources/encode_enhancers_liftover.bed

TSS=workspace/resources/hg38v79_genes_tss_2000.bed
ENH=workspace/resources/encode_enhancers_liftover.bed

sample="MCF10A_H3K27ac"
query_bam=~/CourseData/EPI_data/Module1/MCF10A_resources/${sample}.bam

samtools view -@4 -q 10 -F 1028 $query_bam -c
samtools view -@4 -q 10 -F 1028 $query_bam -L ~/CourseData/EPI_data/Module1/QC_resources/encode_enhancers_liftover.bed -c
samtools view -@4 -q 10 -F 1028 $query_bam -L ~/workspace/peaks/${sample}_peaks.blacklistRemoved.narrowPeak -c
```
- `samtools view -@4 -q 10 -F 1028 $query_bam -L ~/workspace/peaks/${sample}_peaks.blacklistRemoved.narrowPeak -c`
    - `samtools view` parse through our BAM file
    - `-@4` resources to use
    - `-q 10` filter for reads with mapping quality(MAPQ) >10
    - `-F 1028` Remove reads that have the following flags:UNMAP,DUP
    - `$query_bam` Bam of interest
    - `-c ` count the number of reads that fulfil the criteria
    - `-L` perform actions on reads that fall within the specified genomic coordiantes (Bed File)
### Step 5B : Quality Control (Enrichment in key genomic areas) Continued:
```
TSS=workspace/resources/hg38v79_genes_tss_2000.bed
ENH=workspace/resources/encode_enhancers_liftover.bed

for histone in H3K27ac H3K27me3 H3K4me3 Input;
do
    sample="MCF10A_${histone}"
    query_bam=workspace/alignments/${sample}_chr19.sorted.dup_marked.bam
    samtools view -@4 -q 10 -F 1028 $query_bam -c > workspace/qc/col_${histone}
    samtools view -@4 -q 10 -F 1028 $query_bam -L ${TSS} -c >> workspace/qc/col_${histone}
    samtools view -@4 -q 10 -F 1028 $query_bam -L ${ENH} -c >> workspace/qc/col_${histone}

    if [[ "$histone" == "H3K27me3" ]]; then
        peaks=workspace/peaks/${sample}_peaks.blacklistRemoved.broadPeak
        samtools view -@4 -q 10 -F 1028 $query_bam -L workspace/peaks/${sample}_peaks.blacklistRemoved.broadPeak -c >> workspace/qc/col_${histone}
    elif [[ "$histone" == "H3K27ac" || "$histone" == "H3K4me3" ]]; then
        peaks=workspace/peaks/${sample}_peaks.blacklistRemoved.narrowPeak
        samtools view -@4 -q 10 -F 1028 $query_bam -L workspace/peaks/${sample}_peaks.blacklistRemoved.narrowPeak -c >> workspace/qc/col_${histone}
    else 
        echo
    fi
done

paste workspace/qc/col_H3K27ac workspace/qc/col_H3K27me3 workspace/qc/col_H3K4me3 workspace/qc/col_Input > workspace/qc/integers.tsv
paste \
<(awk '{print $1/442829*100}' workspace/qc/col_H3K27ac) \
<(awk '{print $1/1610910*100}' workspace/qc/col_H3K27me3) \
<(awk '{print $1/1512352*100}' workspace/qc/col_H3K4me3) \
<(awk '{print $1/1023265*100}' workspace/qc/col_Input) > workspace/qc/percentages.tsv
```
- Reads enrichment in key region 

| |H3K27ac|H3K27me3|H3K4me3|Input|
|------------- | ------------- | ------------- | ------------- | ------------- |
|Total|442829|1610910|1512352|1023265
|TSS|127577|298326|1087032|194996
|Enhancer|242489|760089|834053|514115
|In Peaks|69584|783294|1239717

 	 	 	 	 
- Reads % enrichment in key region
				
| |H3K27ac|H3K27me3|H3K4me3|Input|
|------------- | ------------- | ------------- | ------------- | ------------- |
|TSS|28.8|18.52|71.88|19.06
|Enhancer|54.76|47.18|55.15|50.24
|In Peaks|15.71|48.62|81.97

## Module 3 - Differential Analysis
### Step1A: Using Bedtools to compare marks
**Code**
```
###Shell###
MCF10A_H3K27ac=workspace/peaks/MCF10A_H3K27ac_peaks.blacklistRemoved.narrowPeak
MCF10A_H3K27me3=workspace/peaks/MCF10A_H3K27me3_peaks.blacklistRemoved.broadPeak
MCF10A_H3K4me3=workspace/peaks/MCF10A_H3K4me3_peaks.blacklistRemoved.narrowPeak

bedtools intersect -u -a ${MCF10A_H3K27ac} -b ${MCF10A_H3K27me3} | wc -l
bedtools intersect -u -a ${MCF10A_H3K27ac} -b ${MCF10A_H3K4me3} | wc -l
```
- `bedtools intersect -u -a ${MCF10A_H3K27ac} -b ${MCF10A_H3K27me3} | wc -l`
  - results in an intersect of `6/988`
- `bedtools intersect -u -a ${MCF10A_H3K27ac} -b ${MCF10A_H3K4me3} | wc -l`
  - results in an intersect of `789/988`
- what do we know about the relationship of H3K27ac vs H3K27me3 vs H3K4me3?
  - H3K27ac and H3K27me3 tend to be antagonistic, hence the very small intersect
  - H3K27ac co-occurs with H3K4me3 at promoters, hence the larger intersect
### Step1B: Using Bedtools to compare samples
**Code**
```
###Shell###
basal_H3K27ac=CourseData/bed/69057.CEEHRC.CEMT0035.H3K27ac.peak_calls.bed
luminal_H3K27ac=CourseData/bed/69095.CEEHRC.CEMT0037.H3K27ac.peak_calls.bed
stromal_H3K27ac=CourseData/bed/69076.CEEHRC.CEMT0036.H3K27ac.peak_calls.bed
lp_H3K27ac=CourseData/bed/69114.CEEHRC.CEMT0038.H3K27ac.peak_calls.bed
MCF10A_H3K27ac=workspace/peaks/MCF10A_H3K27ac_peaks.blacklistRemoved.narrowPeak

bedtools intersect -u -a ${MCF10A_H3K27ac} -b ${basal_H3K27ac} | wc -l
bedtools intersect -u -a ${MCF10A_H3K27ac} -b ${luminal_H3K27ac} | wc -l
bedtools intersect -u -a ${MCF10A_H3K27ac} -b ${stromal_H3K27ac} | wc -l
bedtools intersect -u -a ${MCF10A_H3K27ac} -b ${lp_H3K27ac} | wc -l

paste \
<(ls CourseData/Bed/*H3K4me3* | xargs -I {} sh -c "bedtools intersect -u -a workspace/peaks/MCF10A_H3K4me3_peaks.blacklistRemoved.narrowPeak -b {} | wc -l") \
<(ls CourseData/Bed/*H3K27ac* |  xargs -I {} sh -c "bedtools intersect -u -a workspace/peaks/MCF10A_H3K27ac_peaks.blacklistRemoved.narrowPeak -b {} | wc -l") \
<(ls CourseData/Bed/*H3K27me3* |  xargs -I {} sh -c "bedtools intersect -u -a workspace/peaks/MCF10A_H3K27me3_peaks.blacklistRemoved.broadPeak -b {} | wc -l")
```
- Intersect numbers:

| |H3K4me3|H3K27ac|H3K27me3|
|------------- | ------------- | ------------- | ------------- |
|MCF10A|1406|988|4797|
|Intersecting Basal|1039|656|2420|
|Intersecting Stromal|978|717|2604|
|Intersecting Luminal|1024|766|2430|
|Intersecting Luminal Progenitor|1099|778|2496|

- MCF10A is luminal progenitor like, how is that relationship reflect in the epigenetic landscape?
  - higher amount of intersect in permissive marks of H3K4me3 and H3K27ac
### Step1D: Using Bedtools to compare binary conditions/models
```
###Shell###
condA_rep1=CourseData/triplicates/peaks/CondA.Rep1_peaks.narrowPeak  
condB_rep1=CourseData/triplicates/peaks/CondB.Rep1_peaks.narrowPeak
condA_rep2=CourseData/triplicates/peaks/CondA.Rep2_peaks.narrowPeak  
condB_rep2=CourseData/triplicates/peaks/CondB.Rep2_peaks.narrowPeak
condA_rep3=CourseData/triplicates/peaks/CondA.Rep3_peaks.narrowPeak
condB_rep3=CourseData/triplicates/peaks/CondB.Rep3_peaks.narrowPeak

bedtools intersect -u -a ${condA_rep1} -b ${condA_rep2} | wc -l
bedtools intersect -u -a ${condA_rep1} -b ${condB_rep2} | wc -l 
bedtools intersect -v -a ${condA_rep1} -b ${condA_rep2} | wc -l
bedtools intersect -v -a ${condA_rep1} -b ${condB_rep2} | wc -l

bedtools intersect -u -a ${condA_rep1} -b ${condA_rep2} ${condA_rep3} | wc -l

bedtools intersect -u -a ${condA_rep1} -b ${condA_rep2} ${condA_rep3} -f 0.5 -F 0.5 | wc -l

bedtools intersect -wao -a ${condA_rep1} -b ${condA_rep2} ${condA_rep3} | head
```

- `bedtools intersect -u -a ${condA_rep1} -b ${condA_rep2} | wc -l `
    - counting the number of condA_rep1 peaks that intersect condA_rep2
    - returns `1191`
- `bedtools intersect -u -a ${condA_rep1} -b ${condB_rep2} | wc -l `
    - counting the number of condA_rep1 peaks that intersect condB_rep2
    - returns `1093`
- `bedtools intersect -v -a ${condA_rep1} -b ${condA_rep2} | wc -l`
    - counting the number of condA_rep1 peaks that do not intersect condA_rep2
    - return `50`
- `bedtools intersect -v -a ${condA_rep1} -b ${condB_rep2} | wc -l`
    - counting the number of condA_rep1 peaks that do not intersect condB_rep2
    - returns `148`
- as expected our replicates of matching conditions have more in common
- `bedtools intersect -u -a ${condA_rep1} -b ${condA_rep2} ${condA_rep3} | wc -l`
   - counting the number of condA_rep1 peaks that intersect condA_rep2 or condA_rep3
- `bedtools intersect -wao -a ${condA_rep1} -b ${condA_rep2} ${condA_rep3} | head`
   - specify `-wao` returns the original line of `${condA_rep1}` and the element it intersects
   - additionally adds an identify column for the database and number of base pairs overlapping
- `bedtools intersect -u -a ${condA_rep1} -b ${condA_rep2} ${condA_rep3} -f 0.5 -F 0.5 | wc -l`
   - the flag `-f 0.5` adds the conditions that intersects are only counted when 50% overlap of A occurs
   - the flag `-F 0.5` adds the conditions that intersects are only counted when 50% overlap of B occurs
   - if we remove one of the flags, how does the number change?
   - what if we wanted integer threshold instead of percentage?
- `bedtools intersect -wao -a ${condA_rep1} -b ${condA_rep2} ${condA_rep3} | head`
   - specify `-wao` returns the original line of `${condA_rep1}` and the element it intersects
   - additionally adds an identify column for the database and number of base pairs overlapping


### Step1C: Other useful bedtool functions
```
###Shell###
bedtools closest
bedtools map
bedtools merge
bedtools cluster
```
### Step2: Differential peaks utilizing triplicates and DiffBind
**Code :**
```
###R###
library(DiffBind)
setwd("/Users/esu/Desktop/work/epiworkshop/")
read.csv("CourseData/triplicates/triplicates.csv")
MCF10A <- dba(sampleSheet="CourseData/triplicates/triplicates.csv")
MCF10A <- dba(sampleSheet=samples)
MCF10A <- dba.count(MCF10A, bUseSummarizeOverlaps=TRUE)
dba.plotPCA(MCF10A, attributes=DBA_CONDITION,label=DBA_ID)
plot(MCF10A)
MCF10A <- dba.contrast(MCF10A, categories=DBA_CONDITION)

MCF10A <- dba.analyze(MCF10A, method=DBA_EDGER)

analyzed_peaks <- dba.report(MCF10A, method=DBA_EDGER, fold=1)

dba.plotMA(MCF10A, bXY=TRUE , method=DBA_EDGER, fold=1)

write.table(analyzed_peaks, file="workspace/diffBind/differential_peaks.tsv", sep="\t", quote=F, row.names=F, col.names=F)
```
- `library(DiffBind)`
    - we load R package [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)
- `setwd("/Users/esu/Desktop/work/epiworkshop/")`
    - set our working directory
- `read.csv("CourseData/triplicates/triplicates.csv")`
    - read in our csv
    - let's inspect the columns
- `MCF10A <- dba(sampleSheet=samples)`
    - read our samplesheet into the `dba` object that will be saved as `MCF10A`
- `MCF10A <- dba.count(MCF10A, bUseSummarizeOverlaps=TRUE)`
    - count the number of fragments that intersect with peaks
    - `bUseSummarizeOverlaps=TRUE` indicates the counting module to be used. `SummarizeOverlaps` comes from [GenomicAlignments](https://www.rdocumentation.org/packages/GenomicAlignments/versions/1.8.4/topics/summarizeOverlaps-methods).
- `dba.plotPCA(MCF10A, attributes=DBA_CONDITION,label=DBA_ID)`
    - generate a principle component analysis using data from our object `MCF10A` where the annotations are `DBA_CONDITION` and the labelling is `DBA_ID`
- `plot(MCF10A)`
    - generate a heatmap with correlation and dendrogram
    - should note the correlation is based on score of overlap and not pearson and spearman, should recalculate
- `MCF10A <- dba.contrast(MCF10A, categories=DBA_CONDITION)`
   - declares what are the conditions for our differential groups
   - `categories=DBA_CONDITION` the category we want to compare
- `MCF10A <- dba.analyze(MCF10A, method=DBA_EDGER)`
   - perform an analysis based on the `contrast` we previously established.
   - `method=DBA_EDGER`, analysis engine is a library called `edgeR`
   - note for our specific example `deseq2` does not work. `Deseq2` has a built in check for variablity which our synthetic dataset is lacking

- `analyzed_peaks <- dba.report(MCF10A, method=DBA_EDGER, fold=1)`
    - report the peaks identified by `DBA_EDGER` to be significant and have an absolute fold change `>1`
- `dba.plotMA(MCF10A, bXY=TRUE , method=DBA_EDGER, fold=1)`
    - generates Scatter plot
    - `method=DBA_EDGER` fetch results based on our previous analysis using `edgeR`
    - `bXY=TRUE` produces a scatter plot, `FALSE` produces a MA plot
    - `fold=1` report differential positions that meet fold change threshold

- `write.table(analyzed_peaks, file="workspace/diffBind/differential_peaks.tsv", sep="\t", quote=F, row.names=F, col.names=T)`
   - save our differential peaks to a TSV
   - `sep="\t"` the seperator to be used
   - `col.names=T` include column names
   - `row.names=F` include row names
   - `quote=F` if we want to include quotations around values
### Step3: Differential peaks utilizing Fold change and significance 