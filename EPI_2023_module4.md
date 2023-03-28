---
layout: tutorial_page
permalink: /EPI_2021_Module4_lab
title: EPI 2021 Module 4 Lab
header1: Workshop Pages for Students
header2: Introduction to WGBS and Analysis 
image: /site_images/CBW_Epigenome-data_icon.jpg
home: https://bioinformaticsdotca.github.io/EPI_2021
---

-----------------------

# Introduction to WGBS and Analysis 

*by Guillaume Bourque, PhD and Jose Hector Galvez*

## Contents: 

1. [Introduction](#introduction)

2. [Mapping Tutorial](#map_tutorial)

    2.1 [Getting Started](#getting_started)
  
    2.2 [Map Data using Bismark](#map_bismark)
    
    2.3 [Repeat Alignment for All Datasets](#repeat)
    
    2.4 [Load Data and Explor with IGV](#load_igv)
    
    2.5 [Generate Methylation Profiles](#methylation_profiles)
    
3. [Differential Methylation Analysis in MethylKit](#methylkit)

    3.1 [Load R and MethylKit](#load_r)
    
    3.2 [Import the Alignment Data into methylKit](#r_import)
    
    3.3 [Find Differentially Methylated Regions](#r_dmr)


---

<a name="introduction"></a>

## 1. Introduction 

#### Description of the lab:

This module will cover the basics of Whole Genome Bisulfite-Sequencing (WGBS) data analysis including data visualization in IGV.


#### Objectives: 

1) Learn how to align WGBS data using `Bismark`
2) Learn how to generate a methylation profile with `Bismark`
3) Learn how to open alignments and methylation profiles in the `IGV` genome browser
4) Learn how to perform a basic differential methylation analysis with `MethylKit`


#### Local software that we will use:
Before you begin, make sure you have the following programs ready in your local computer: 

* A connection to the `EPI_2021` AWS instance
* An internet browser 
* IGV
* R (or RStudio), *optional*


---

<a name="map_tutorial"></a> 

## 2. Mapping Tutorial

<a name="getting_started"></a>

### 2.1 Getting Started


##### Prepare Directory for the Lab

```bash
mkdir -p ~/workspace/module4
cd ~/workspace/module4

```

##### Save Genome Location

```bash
GENOME=~/CourseData/EPI_data/module4/Homo_sapiens.GRCh38.chr19
export GENOME

```

This will define a variable `$GENOME` that will simplify future commands. 

---

#### Locate the Data for the Workshop

```bash
WGBS_DATA=~/CourseData/EPI_data/module4/data
export WGBS_DATA

```

This will define a variable `$WGBS_DATA` that will simplify future commands. 

##### Check the Files

**Question 1** 

*Type the following command: `ls $WGBS_DATA`, what do you see?*
  
**Answer** 
  
You should see something similar to this:


```{:.output}
WGBS.A34002.137160.chr19.1.fastq.gz
WGBS.A34002.137160.chr19.2.fastq.gz
WGBS.A34002.137487.chr19.1.fastq.gz
WGBS.A34002.137487.chr19.2.fastq.gz
WGBS.A34002.137488.chr19.1.fastq.gz
WGBS.A34002.137488.chr19.2.fastq.gz

```


These are the files that will be used for the workshop. They contain a subset of WGBS reads from CEMT sample `CEMT0007`, which is a mammary gland epithelial cell line [(more information here)](https://ega-archive.org/datasets/EGAD00001000685). 

**Question 2**

*What do the ".1" and ".2" in the file names mean?*

**Answer**
  
They represent the `read1` and `read2` of the paired end reads. 


---

 <a name="map_bismark"></a>
 
### 2.2 Map Data using Bismark

We will now process and map the reads using the Bismark WGBS aligner [(more info here)](https://www.bioinformatics.babraham.ac.uk/projects/bismark/).

#### Map the first dataset using Bismark

To simplify the work, we will process the datasets one at a time. To align the first dataset, do the following: 


```{bash}
cd ~/workspace/module4
bismark --multicore 4 --bowtie2 $GENOME/genome/bismark_index \
    -1 $WGBS_DATA/WGBS.A34002.137160.chr19.1.fastq.gz -2 $WGBS_DATA/WGBS.A34002.137160.chr19.2.fastq.gz
    
```

**Question 3**

*What do all the options in the command mean?* (Hint check the help by using `bismark --help`)
  
**Answer** 
  
- The `--multicore 4` option is to do multithreaded processing to improve speed.

- The `--bowtie2` option is to use the mapping algorithm from bowtie2.

- The `$GENOME/genome/bismark_index` specifies the location of the index for the reference genome to use. This uses the `$GENOME` variable we defined previously. 

- The `-1 $WGBS_DATA/WGBS.A34002.137160.chr19.1.fastq.gz` specifies the location of read 1. Idem for `-2` which specifies read 2. This uses the `$WGBS_DATA` variable we defined previously. 


For more details, please refer to the Bismark [user guide](http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf). 


---

This step will take a few minutes to run for this reduced dataset. A dataset spanning a full genome will take several hours. 

For your own datasets, make sure you have enough computing walltime to run the alignment. 


**While you wait for the results, ask any questions you have up to this point to the instructors.**


#### Check files

At the end of the alignment, you should have the following files saved into your workshop folder: 

```{:output}
WGBS.A34002.137160.chr19.1_bismark_bt2_pe.bam  WGBS.A34002.137160.chr19.1_bismark_bt2_PE_report.txt
```

Let's look at the report: 

```{bash}
less WGBS.A34002.137160.chr19.1_bismark_bt2_PE_report.txt

```

**Question 4**

*What was the mapping efficiency? What percent of C's were methylated in CpG context?*
  
**Answer** 
  
According to the report: 


```{:output}
...
Mapping efficiency:     92.4% 
...
C methylated in CpG context:    57.4%
C methylated in CHG context:    0.6%
C methylated in CHH context:    0.5%
C methylated in unknown context (CN or CHN):    3.5%
...

```

Close the report by pressing `q`. 

---

#### Prepare files for loading in IGV

We need to sort the `bam` file and prepare an index so we will be able to load it in IGV. We will use the program `samtools` for this.

```{bash}
samtools sort WGBS.A34002.137160.chr19.1_bismark_bt2_pe.bam -o WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam 
samtools index WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam

```

#### Check Files

At the end, you should have the following files: 

```
WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam
WGBS.A34002.137160.chr19.1_bismark_bt2_pe.bam         
WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam.bai
WGBS.A34002.137160.chr19.1_bismark_bt2_PE_report.txt
```

---

<a name="repeat"></a>

### 2.3 Repeat Alignment for All Datasets

**Question 5**

*How would you repeat the alignment with the other datasets?*

**Answer** 

This is the command to run bismark on the two other samples: 

```{bash}
cd ~/workspace/module4

bismark --multicore 4 --bowtie2 $GENOME/genome/bismark_index \
    -1 $WGBS_DATA/WGBS.A34002.137487.chr19.1.fastq.gz -2 $WGBS_DATA/WGBS.A34002.137487.chr19.2.fastq.gz
    
bismark --multicore 4 --bowtie2 $GENOME/genome/bismark_index \
    -1 $WGBS_DATA/WGBS.A34002.137488.chr19.1.fastq.gz -2 $WGBS_DATA/WGBS.A34002.137488.chr19.2.fastq.gz
    
```

Remember, for the command to work, both `$GENOME` and `$WGBS_DATA` need to be defined. 

This is the command to prepare the samples for IGV (sort and index): 

```{bash}
samtools sort  WGBS.A34002.137487.chr19.1_bismark_bt2_pe.bam -o WGBS.A34002.137487.chr19.1_bismark_bt2_pe_sorted.bam
samtools index WGBS.A34002.137487.chr19.1_bismark_bt2_pe_sorted.bam

samtools sort  WGBS.A34002.137488.chr19.1_bismark_bt2_pe.bam -o WGBS.A34002.137488.chr19.1_bismark_bt2_pe_sorted.bam
samtools index WGBS.A34002.137488.chr19.1_bismark_bt2_pe_sorted.bam

```


---

<a name="load_igv"></a>

### 2.4 Load Data and Explore using IGV

While you wait for the previous steps to finish executing, it is a good idea to begin exploring the alignments. 

#### Copy Files to Your Local Computer to View in IGV (optional)

Retrieve the files called `WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam` and `WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam.bai` from the server using your internet browser and the public IP address of your AWS instance. 


#### Launch IGV on your computer

If you haven't installed it yet, please get it here [IGV download](http://software.broadinstitute.org/software/igv/download).

Make sure that the human genome is selected in the top left corner. It should read: **Human (hg38)**.

Load your sorted `bam` file in IGV using `File -> Load from file`. *For this to work, you need to have the index file (`.bai`) in the same location as the `bam` file.* You can also load this directly from the internet using the URL of this file you got from your public IP page, using `File -> Load from URL`. 

Now, on IGV, go to the following location:

```
chr19:43,375,889-45,912,052
```

And zoom in until you see something.

For instance, try the following window:

```
chr19:44,527,387-44,536,873
```

You should see something like this:

<img src="https://github.com/bioinformatics-ca/EPI_2021/blob/master/img/region2.png?raw=true" alt="Region" width="750" />

*If it looks different, can you change the way the colors are displayed?*

*Which section of which chromosome is covered by this dataset?* 

*Can you see any interesting patterns in the coverage?* 

---

<a name="methylation_profiles"></a>

### 2.5 Generate Methylation Profiles

So far we have only mapped the reads using Bismark. We can generate methylation profiles using the following command:

```{bash}
cd ~/workspace/module4

bismark_methylation_extractor --bedGraph WGBS.A34002.137160.chr19.1_bismark_bt2_pe.bam

```

**Question 6** 

*How would you do the same for the other replicates?*

**Answer**

These are the commands that you should use:

```{bash}
cd ~/workspace/module4

bismark_methylation_extractor --bedGraph WGBS.A34002.137487.chr19.1_bismark_bt2_pe.bam
bismark_methylation_extractor --bedGraph WGBS.A34002.137488.chr19.1_bismark_bt2_pe.bam

```


---

Make sure that *all* the files produced so far are displayed in your internet browser. 

**While you wait for all the steps to finish, you can ask the instructors any questions you might have up until this point.**

---

Load all the files in IGV using `File -> Load from file` or `File -> Load from URL`.

At this point, if you load the region `chr19:44,527,387-44,536,873` you should see something like

<img src="https://github.com/bioinformatics-ca/EPI_2021/blob/master/img/region2_full.png?raw=true" alt="Region" width="750" />

This promoter looks to be hypomethylated. 

*Can you find a promoter that is hypermethylated?*

How about `chr19:45,637,715-45,657,380`?

*How would you look for a CpG island using this view of the data?*

Keep exploring the files and see if you can find differences in the methylation profiles of the samples. 

---

<a name="methylkit"></a>

## 3. Differential Methylation Analysis in MethylKit

The following section will use the Bioconductor package `methylKit` to do a differential methylation analysis. 
You can do it in your own computer (if you have installed `R` and `methylKit`) or in the AWS instance. 

To install `methylKit` locally on your computer, make sure you have a recent version of `R` and 
follow the instructions in this [page](https://bioconductor.org/packages/release/bioc/html/methylKit.html). 

---

<a name="load_r"></a>

### 3.1 Load R and MethylKit

If you are working in AWS, you will need to load `R`. The image we provide already has the libraries we need. 

To launch `R` simply type the following to your terminal: 

```{bash}
cd ~/workspace/module4
R

```

If you did this properly, the following message will be displayed and your prompt will change from **`ubuntu@ip-00-00-00-0:~/workspace/module4$`** to **`>`**: 

```{:output}
R version 4.0.5 (2021-03-31) -- "Shake and Throw"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
...

```

Once you have successfully launched `R`, you can load `methylKit` with the following command: 

```{r}
library("methylKit") 

```

---

<a name="r_import"></a>

### 3.2 Import the Alignment Data into methylKit

#### Process Bismark Alignments

To read the alignment data into `methylKit`, run the following command: 

```{r}
methRaw.160 = processBismarkAln( location = "WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam",
                         sample.id="A34002.137160", assembly="hg38", 
                         read.context="CpG", save.folder="methylkit")

```

This command will import the data into a format that is readable by `methylKit`. At the same time, it will save two files under the `methylkit` directory with the information so that it is easy to load again at any time: 

```
methylkit/A34002.137160_CpG_conversionStats.txt
methylkit/A34002.137160_CpG.txt

```

If everything goes well and you see the files, do the same for the other two samples: 

```{r}
methRaw.487 = processBismarkAln( location = "WGBS.A34002.137487.chr19.1_bismark_bt2_pe_sorted.bam",
                         sample.id="A34002.137487", assembly="hg38", 
                         read.context="CpG", save.folder="methylkit")
                         
methRaw.488 = processBismarkAln( location = "WGBS.A34002.137488.chr19.1_bismark_bt2_pe_sorted.bam",
                         sample.id="A34002.137488", assembly="hg38", 
                         read.context="CpG", save.folder="methylkit")

```

---

#### Create a MethylKit Object

Now that all the samples have been read with `methylKit`, you can create a file list to make it easier to load the full dataset as a `methylkit` object. For the purposes of this tutorial, we will consider that samples belong to two experimental groups: `A34002.137160` as the control group (treatment = 0) and `A34002.137487` & `A34002.137488` as the treatment group (treatment = 1). We use the `methRead()` function to create our object, as shown below: 

```{r}
file.list = list( file.path("methylkit", "A34002.137160_CpG.txt"),
                  file.path("methylkit", "A34002.137487_CpG.txt"),
                  file.path("methylkit", "A34002.137488_CpG.txt") )

myobj = methRead(file.list,
           sample.id=list("A34002.137160","A34002.137487","A34002.137488"),
           assembly="hg38",
           treatment=c(0,1,1),
           context="CpG",
           mincov = 10
           )

```

**Question 7** 

*What do all the options in the `methRead()` command mean?*
  
**Answer** 
   
  
- `file.list` object points to the location of the input data in a MethylKit format. 

- `sample.id` points to a list with the appropriate sample name for each file. 

- `assembly` specifies which build of the human reference genome is used. 

- `treatment` specifies which sample belongs to each experimental group. 

- `context` specifies the methylation context. 

- `mincov` specifies the minimum coverage required to be included in the object. 


For more details, please refer to the MethylKit [user guide](https://bioconductor.org/packages/release/bioc/manuals/methylKit/man/methylKit.pdf).


---


If the files were loaded properly, you can check the object you just created by running the following command: 


```{r}
myobj

```

Which should output the following message followed by previews of the contents of the object: 

```
methylRawList object with 3 methylRaw objects
...
```

You can also get basic statistics on your object by using the following command: 

```{r}
getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)

```

--- 

### 3.3 Find Differentially Methylated Regions

<a name="r_dmr"></a>

#### Merge Samples

Before doing any additional analysis, `methylKit` needs to determine which methylated bases have sufficient coverage in all samples so they can be compared. To do that, the samples should be merged with the `unite()` function. This function has a parameter `destrand=` that is turned off by default. We will set the `destrand` option to `TRUE` which will merge the coverage of both strands. When doing your own analyses, be aware that for some kinds of methylation analyses (such as CpH methylation) results are strand-specific, so this option should be used carefully. 

```{r}
meth = unite(myobj, destrand=TRUE)

```

#### Perform Differential Methylation Analysis

The standard function for Differential Methylation Analysis on `methylKit` is `calculateDiffMeth()`. It takes any **merged** methylkit object as input. Depending on the number of replicates, it uses either Fisher’s exact or logistic regression to calculate P-values. It also, automatically produces Q-values, which are a kind of adjusted P-value. To use it with the results we obtained before, run the following command: 

```{r}
myDiff = calculateDiffMeth(meth)

```

To check the output, just type `myDiff` and read the summary. If you want an example of the output, check the solution below. 


**Example**  

This is what the output looks like: 

```{:output}
methylDiff object with 2941 rows
--------------
    chr    start      end strand       pvalue     qvalue  meth.diff
1 chr19 42002896 42002896      + 3.271268e-01 0.69569299   8.951407
2 chr19 42002978 42002978      + 1.912989e-01 0.60732656 -21.666667
3 chr19 42007251 42007251      + 6.999764e-05 0.03228847 -55.681818
4 chr19 42007255 42007255      + 3.958578e-01 0.75196047 -11.835106
5 chr19 42007283 42007283      + 8.451850e-01 0.91347038  -2.457757
6 chr19 42007314 42007314      + 9.102723e-01 0.92865750  -1.604278
--------------
sample.ids: A34002.137160 A34002.137487 A34002.137488 
destranded TRUE 
assembly: hg38 
context: CpG 
treament: 0 1 1 
resolution: base 
```

---

To filter results by their statistical significance, `methylKit` provides the `getMethylDiff()` function which allows you to extract only the deferentially methylated CpG's that meet a specific Q-value threshold. Additionally, it is also possible to specify whether to keep `hypo` or `hyper` methylated CpG's only. Finally, the `bedgraph()` function allows you to save the the `methylDiff` object into a BedGraph file so you can open it with your genome browser of choice. Let's create two BedGraph files with hypo and hyper methylated CpG's with a Q-value below 0.05 based on the data above: 

```{r}
myDiff.hyper = getMethylDiff(myDiff,qvalue=0.05,difference=10,type="hyper")
bedgraph(myDiff.hyper, file.name = "hyper.CpG.bedGraph", col.name = "qvalue")

myDiff.hypo = getMethylDiff(myDiff,qvalue=0.05,difference=10,type="hypo")
bedgraph(myDiff.hypo, file.name = "hypo.CpG.bedGraph", col.name = "qvalue")

```

Two new files should appear now in your workshop folder: 

```{:output}
~/workspace/module4/hyper.CpG.bedGraph
~/workspace/module4/hypo.CpG.bedGraph
```

#### Bin Results to Obtain Differentially Methylated Regions 

By default, `methylKit` will compute results with an individual CpG resolution. To get Differentially Methylated Regions (**DMR**), you have to bin your results first, using a window size of your choice. The function to do this is `tileMethylCounts()`, which takes a regular methylkit object as input.  In this case, we will create 1000bp bins using the following command: 

```{r}
tiles = tileMethylCounts(myobj,win.size=1000,step.size=1000,cov.bases = 10)

```

As with CpG level results, samples need to be merged before the analysis can continue: 

```{r}
meth.tiles = unite(tiles, destrand=TRUE) 

```

Now, we will use the `calculateDiffMeth()` and `getMethylDiff()` functions to get the DMRs. 

**Questions 8 & 9**

*Do you know how to do it, based on the information above?*

*Based on the number of differentially methylated CpGs you found above, do you anticipate many statistically significant DMRs in your analysis?*
  
**Answer** 
  
  
Use the following commands to perform a `DMR` analysis: 

```{r}
myDiff.tiles = calculateDiffMeth(meth.tiles)

myDiff.tiles.hyper = getMethylDiff(myDiff.tiles,qvalue=0.1,difference=10,type="hyper")
bedgraph(myDiff.tiles.hyper, file.name = "hyper.DMR.bedGraph", col.name = "qvalue")

myDiff.tiles.hypo = getMethylDiff(myDiff.tiles,qvalue=0.1,difference=10,type="hypo")
bedgraph(myDiff.tiles.hypo, file.name = "hypo.DMR.bedGraph", col.name = "qvalue")


```


---

Using the navigation pane, download the bedGraph files you just produced and try to open them with IGV. 

*Do the statistical results match what you had seen before when exploring the data?*

*What interesting genomic features are found close to the DMRs? What could this mean?*


### Congrats, you're done!

You can quit R using the `quit()` or `q()` command. Remember to stop your AWS instance after this lab to avoid unnecessary costs. 


Once you are finished make sure you download all the files you need and continue exploring on IGV. 



