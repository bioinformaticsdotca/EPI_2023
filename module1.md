# Introduction to ChIP-Seq Analysis
by Martin Hirst and Edmund Su


## Table of Contents
1. Breakdown of markdown file
2. Server side resources
3. Module 1. Alignment
    - Perform FastQC and Inspect results
    - Index a Reference and popular reference resources
    - Align Data using BWA
    - Coordinate Sort BWA alignment file
    - MarkDuplicates in an Alignment File
    - Generate Flagstats and Stats and how to interpret for QC
    - Process multiple files using the above steps
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

## Module 1. Alignment
### Step 1: Setup
- It's always good practice to organize your directories beforehand
- This avoids file conflicts and accidental data loss.
- We'll be making a subdirectory for each major step.

**Code:**
```
###Shell###
mkdir -p ~/workspace/module123/BWA_index
mkdir -p ~/workspace/module123/alignments
mkdir -p ~/workspace/module123/fastqc
mkdir -p ~/workspace/module123/stats
mkdir -p ~/workspace/module123/peaks
mkdir -p ~/workspace/module123/bigBed
mkdir -p ~/workspace/module123/bigWig
mkdir -p ~/workspace/module123/resources
mkdir -p ~/workspace/module123/diffBind
mkdir -p ~/workspace/module123/deeptools
mkdir -p ~/workspace/module123/qc
```
- `mkdir -p ~/workspace/module123/BWA_index`
  -  `mkdir` creates a directory
  - `-p` creates the "parent" if the initial parent directory did not exist
  - how could we simplfy the command?
  - `mkdir -p ~/workspace/module123/{BWA_index,alignments,fastqc,stats,peaks,bigBed,bigWig,resources,diffBind,deeptools,qc}`



### Step 2A: Retrieve a reference genome
- We'll need a reference genome to align our sequences to.
- A great resource is [UCSC's genome repository](https://hgdownload.soe.ucsc.edu/downloads.html) containing latest `Hg38` and older `hg19` human genome.
- For our tutorial, we'll be working `chr19` of `Hg38`/`Grch38`. This is to speed demonstration purposes.

**Code:** 
```
###Shell###
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr19.fa.gz -O ~/workspace/module123/BWA_index/chr19.fa.gz
gunzip ~/workspace/module123/BWA_index/chr19.fa.gz -c > ~/workspace/module123/BWA_index/chr19.fa
```
- Runtime : <1 minute

- `wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr19.fa.gz -O ~/workspace/module123/BWA_index/chr19.fa.gz`

    - `Wget` is a command to pull online files.
    - We're instructing it to pull `chr19.fa.gz` from the example address.
    - `-O ~/workspace/module123/BWA_index/chr19.fa.gz` instructs where and what we would like to save our file as.
    - _What were to happen if we ran the command without `-O`?_
    - Not all systems carry/support `wget`, `Curl` is another useful alternative.

- `gunzip ~/workspace/module123/BWA_index/chr19.fa.gz -c > ~/workspace/module123/BWA_index/chr19.fa`
    - `gzip` is a data compression format useful for reducing storage impacts
    - `gunzip` is the opposite decompresses `gzipped` data
    - `-c` is a STDOUT redirect. Normally running `gunzip file.gz` will turn `file.tsv.gz` into `file.tsv`. Running `-c` allows us to save the contents elsewhere. We do this we can compare the size difference.
    - compare `ls ~/workspace/module123/BWA_index/chr19.fa.gz -ilth` vs `ls ~/workspace/module123/BWA_index/chr19.fa -ilth` and _note the size difference_

> [!NOTE]
> Prior to any tool usage, it's good practice to explore the tool by invoking the help command. This may vary depending on tool. This is important as the same `-O` flag could have differing functions in different tools.
### Step 2B: Index the genome
- Ready the genome fasta file by creating databases allowing tools to quickly access different parts of the file
- Indexing is only required once per genome.
- That being said, because there are different versions of `Hg38`/`GRCh38` (such as those with or [without ALT contigs and EBV](https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19)) each version would need their own indices.
- Additionally index files are typically program/tool specific.
- Consortias will also have available the reference genome and related resources. E.g. https://ewels.github.io/AWS-iGenomes/

> [!NOTE]
> [A newer version of BWA-mem](https://github.com/bwa-mem2/bwa-mem2) does exist, however we'll be using the older version due to a bug that requires significant memory for indexing.

**Code:** 
```
###Shell###
mkdir -p ~/workspace/module123/BWA_index
bwa index ~/workspace/module123/BWA_index/chr19.fa > ~/workspace/module123/BWA_index/index.log
```
- Runtime : <1 minute
- `bwa index ~/workspace/module123/BWA_index/chr19.fa`
    - Indexing step creates a database to enable quick retrieve of sequences
    - take a look inside your folder `ls ~/workspace/module123/BWA_index/`
    - For the lab, we're utilizing [BWA](https://bio-bwa.sourceforge.net/bwa.shtml) for alignment.
### Step 3A: FastQC and interpretation
- Checking the quality of your FASTQs is sanity check step (garbage in garbage out)

**Code:** 
```
###Shell###
fastq_input=~/CourseData/EPI_data/module123/fastq/MCF10A.Input.chr19.R1.fastq.gz
fastq_h3k4me3=~/CourseData/EPI_data/module123/fastq/MCF10A.H3K4me3.chr19.R1.fastq.gz
mkdir -p ~/workspace/module123/fastqc
fastqc -t 8 ${fastq_input} -o ~/workspace/module123/fastqc > ~/workspace/module123/fastqc/input.log 2>&1
fastqc -t 8 ${fastq_h3k4me3} -o ~/workspace/module123/fastqc > ~/workspace/module123/fastqc/h3k4me3.log 2>&1
```
- `fastqc`
    - `-t 8` specifies the number of threads we want to program to utilize
    - `-o` specifies our output directory
### Step 3A: FastQC and interpretation
- Let's take a look at our FASTQC results and compare

**Code:**
```
###Browser###
http://main.uhn-hpc.ca/module123/fastqc/MCF10A.Input.chr19.R1_fastqc.html
http://main.uhn-hpc.ca/module123/fastqc/MCF10A.H3K4me3.chr19.R1_fastqc.html
```
> [!NOTE]
> The URLs link to the TA's instance, make sure to replace the domain with your own.
- Observe the summary report, note the warnings and errors.
- we can also look at reports at https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- `Per base sequence quality`
  - a representation of quality across all reads where the ideal is average,median,25th and 75th resides within the green good quality space
  - X-axis is the BP position in your read
  - Y-axis is the base quality score
  - yellow boxes represent 25th and 75th percentiles
  - whiskers represent 10th and 90th percentiles
  - read line represents median
  - blue line represents average
  - For illumina reads, the first couple of bases tend to be poor, quality normalizes for majority of the read and slowly worsens approaching end of read
  - possible remedies are to trim reads
- `Per tile sequence quality`
  - a representation of quality by tile [(a discrete subdivision with the lane of a flowcell)](https://gtk-teaching.github.io/NGS-intro/fig/hiseq-flow-cell.png) where heat is bad
  - would indicate an issue with the flowcell
- `Per sequence quality score`
  - a distribution of mean read quality score where the peak and body should be as far right as possible
- `Per base sequence content`
  - The % of AGCT across reads according to read position
  - generally expect an equal/even except for the first couple of base pairs (due to sequencing calibration)
  - Certain assays would change this for example : Bisulphite-Sequencing (unmethylated Cs become Us) or amplicon-sequencing (fixed ratios)
  - look at H3K4me3, why is it so high in GC? Promoters!
- `Per sequence GC content`
  - a distrubtion of reads' GC content, where warning or failures are issued based on deviation of 15% or 30%
  - A shift in distribution could mean sequencing bias or contamination
- `Per base N content`
  - B/C `N`s are low confidence calls by the sequencer, we expect very few instances
- `Sequence Length Distribution`
  - With Illumina sequencing, expect uniform distribution
- `Sequence duplication levels`
  - Using the first 200,000 to find duplicates: expect most reads to unique (left side), spikes otherwise indicate contamination or PCR artifacts
- `Overrepresented sequences`
  - a table breakdown of the sequences found above in `Sequence duplication levels`
- `Adapter Content`
  - detects if commonly used adapters remain in read sequence

### Step 4: Alignment
- Mapping the read pairs to a position of the reference genome

**Code:** 
```
###Shell###
ref=~/workspace/module123/BWA_index/chr19.fa
read1=~/CourseData/EPI_data/module123/fastq/MCF10A.Input.chr19.R1.fastq.gz
read2=~/CourseData/EPI_data/module123/fastq/MCF10A.Input.chr19.R2.fastq.gz
sample=MCF10A_input_chr19
bwa mem -M -t 4 ${ref} ${read1} ${read2} 2>~/workspace/module123/alignments/${sample}.alignment.log | samtools view -hbS -o ~/workspace/module123/alignments/${sample}.bam
```
- run time ~2 min
- The command can be broken down into the following pseudo code `ALIGNMENT | SAMtoBAM_conversion`. The `|` operate streams the results from the first command into the second.

- `bwa mem -M -t 4 ${ref} ${read1} ${read2} 2>~/workspace/module123/alignments/alignment.log`
    - - `|` the pipe delimiter passes the results to the next step `ACTION A | ACTION B | ACTION C` (like an assembly line)
    - `bwa mem` while BWA has many alignment algorithms, `mem` is best suited for efficiently handling >70bp paired end reads.
    - `-M` alters flagging of chimeric reads. By default chimeric reads are flagged as `SUPPLEMENTARY` (partially mapping), the option instead turns them to `SECONDARY` (multi-mapping). Needed to GATK/PICARD support downstream
    - `-t 4` resource specification
    - `2>~/workspace/module123/alignments/alignment.log` data outputted occurs in two streams `1` (the data) and `2` (debug messages, warnings and status updates). We redirect `2` to a log file.
    - What do we see in the logs?
> [!NOTE]
> Logs contain valuable information, helpful for debugging.
- `samtools view -hbS -o ~/workspace/module123/alignments/${sample}.bam`
    - `samtools view` a tool to read our `SAM`,`BAM`,`CRAM` files
    - `-hbS` include header, output as `BAM`, input is `SAM`
    - `-o` specify what we want to save our file as
    - take a look our new `BAM` file. Note the header and the body.
- What other important functions of `samtools`?
    - CRAM conversion
    - samtools quickcheck
    - samtools view header
    - samtools flags
### Step 5: Coordinate Sort
- Rearrange our alignments by coordinates

**Code:**
```
###Shell###
sample=MCF10A_input_chr19
samtools sort -@8 ~/workspace/module123/alignments/${sample}.bam -o ~/workspace/module123/alignments/${sample}.sorted.bam
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
java -jar /usr/local/picard.jar MarkDuplicates \
I=~/workspace/module123/alignments/${sample}.sorted.bam \
O=~/workspace/module123/alignments/${sample}.sorted.dup_marked.bam \
M=~/workspace/module123/alignments/${sample}.dup_marked.output.log \
ASSUME_SORTED=TRUE \
VALIDATION_STRINGENCY=LENIENT \
> ~/workspace/module123/alignments/${sample}.dup_marked.error.log
``` 
- `java -jar /usr/local/picard.jar MarkDuplicates I=~/workspace/module123/alignments/${sample}.sorted.bam O=~/workspace/module123/alignments/${sample}.sorted.dup_marked.bam M=~/workspace/module123/alignments/${sample}.dup_marked.output.log ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT`
    - `java -jar /usr/local/picard.jar` produced by BroadInstitute, [Picard tools are a toolset](https://broadinstitute.github.io/picard/) for HTS/NGS data
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
mkdir -p  ~/workspace/module123/stats
samtools flagstat ~/workspace/module123/alignments/${sample}.sorted.dup_marked.bam > ~/workspace/module123/stats/${sample}.sorted.dup_marked.flagstat
samtools stats ~/workspace/module123/alignments/${sample}.sorted.dup_marked.bam > ~/workspace/module123/stats/${sample}.sorted.dup_marked.stat
```
 - `samtools flagstat` tallies `FLAGS` and produces a summarized report
 - `samtools stats` collects various metrics on the BAM file include:
    - summary stats similar to `flagstats`
    - distribution of insert sizes
    - distribution of read lengths
    - distribution of Coverage
    - distribution of GC-depth
    - Includes detailed instructions on how to extract parts of interest

### Step 8: Processing the remaining files
- we use a unix for loop to perform all the steps previously mentioned but for the other histone marks and ATAC-seq
**Code:**
```
###Shell##
ref=~/workspace/module123/BWA_index/chr19.fa
for histone in H3K27ac H3K27me3 H3K4me3 ATAC;
    do
    read1=~/CourseData/EPI_data/module123/fastq/MCF10A.${histone}.chr19.R1.fastq.gz
    read2=~/CourseData/EPI_data/module123/fastq/MCF10A.${histone}.chr19.R2.fastq.gz
    sample=MCF10A_${histone}_chr19
    echo "aligning ${histone}"
    bwa mem -M -t 4 ${ref} ${read1} ${read2} 2> ~/workspace/module123/alignments/${sample}.alignment.log | samtools view -hbS -o ~/workspace/module123/alignments/${sample}.bam
    echo "sorting ${histone}"
    samtools sort -@8 ~/workspace/module123/alignments/${sample}.bam -o ~/workspace/module123/alignments/${sample}.sorted.bam
    echo "dupmarking ${histone}"
    java -jar /usr/local/picard.jar MarkDuplicates I=~/workspace/module123/alignments/${sample}.sorted.bam O=~/workspace/module123/alignments/${sample}.sorted.dup_marked.bam M=~/workspace/module123/alignments/${sample}.dup_marked.output.log ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT > ~/workspace/module123/alignments/${sample}.dup_marked.error.log
    echo "calculating stats ${histone}"
    samtools flagstat ~/workspace/module123/alignments/${sample}.sorted.dup_marked.bam > ~/workspace/module123/stats/${sample}.sorted.dup_marked.flagstat
    samtools stats ~/workspace/module123/alignments/${sample}.sorted.dup_marked.bam > ~/workspace/module123/stats/${sample}.sorted.dup_marked.stat
    done
```
- run time ~8 mins
- for loops are different in every language, best practice to review beforehand
- have feedback to indicate where in the process we are e.g. the echo lines
- not show here but also good to capture success and failures 

### Step 9: Clean up!
- Remove temporary files

**Code:**
```
###Shell###
ls ~/workspace/module123/alignments/*.bam | grep -v sorted.dup_marked | xargs -I {} sh -c "echo rm {}; rm {}"
```
- `ls ~/workspace/module123/alignments/*.bam`
    - a look up of `BAM` files in our directory : `~/workspace/module123/alignments`
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


## Server resources
### QC Resources
```
###TSS+/-2kb
mkdir workspace/module123/qc
wget https://www.bcgsc.ca/downloads/esu/touchdown/hg38v79_genes_tss_2000.bed -O workspace/module123/resources/hg38v79_genes_tss_2000.bed

sort -k1,1 -k2,2n workspace/module123/resources/hg38v79_genes_tss_2000.bed > tmp
mv tmp workspace/module123/resources/hg38v79_genes_tss_2000.bed

###Enhancer liftover
wget https://www.bcgsc.ca/downloads/esu/touchdown/encode_enhancers_liftover.bed -O workspace/module123/resources/encode_enhancers_liftover.bed

###Blacklist
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz -O ~/workspace/module123/resources/hg38_blacklist.bed.gz

gunzip ~/workspace/module123/resources/hg38_blacklist.bed.gz
```
- `hg38v79_genes_tss_2000.bed`
    - Generated by downloading Ensemblv79 GTF convert to Bed +/-2kb of TSS. See https://www.biostars.org/p/56280/
- `encode_enhancers_liftover.bed`
    - download various [ChroHMM state7](https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/) and merge
    - converted from hg19 to hg38 using UCSC liftover tool


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
- triplicates were generated from MCF10A_H3K4me3 by choosing a list of exclusive peaks for condA and condB and subsampling replicates accordingly