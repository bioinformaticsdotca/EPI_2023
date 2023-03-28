---
layout: tutorial_page
permalink: /epi2020_module4_lab
title: Epigenomics Lab 4
header1: Workshop Pages for Students
header2: Whole Genome Bisulfite Sequencing
image: /site_images/CBW_Epigenome-data_icon.jpg
home: https://bioinformaticsdotca.github.io/epigenomics_2020
---

-----------------------

# Introduction to WGBS and Analysis 

*by Guillaume Bourque, PhD and Jose Hector Galvez*

#### *Important notes:*

* The username **user01** in all of the commands below is provided as an example. You should always replace it with the user that was assigned to you for this workshop. 

---

## Contents: 

1. [Introduction](#introduction)

    1.1 [Connecting to the Secure Cloud](#connect)

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

* An internet browser or `ssh` (if you want to follow using only a terminal)
* IGV

---
<a name="connect"></a>
### 1.1 Connecting to the Secure Cloud

The easiest way to connect to the secure cloud is using `jupyterlab` through a web browser. To do that, open your browser and type `jupyter.cbw-oct-2020.calculquebec.cloud` in the URL field. Type your username and password in the appropriate fields. Then you should see a window like this: 

<img src="https://github.com/bioinformatics-ca/Epigenomics_2020/blob/master/JupyterLab_job.png?raw=true" alt="JupyterLab Job" width="750" />

Select the options shown in the image above (make sure you select `JupyterLab` as the **User Interface**). Once you do that you should see a window like this: 

<img src="https://github.com/bioinformatics-ca/Epigenomics_2020/blob/master/JupyterLab_launcher.png?raw=true" alt="JupyterLab Launcher" width="750" />

Click on the `Terminal` button from the options presented. You should see a window like this:

<img src="https://github.com/bioinformatics-ca/Epigenomics_2020/blob/master/JupyterLab_terminal.png?raw=true" alt="JupyterLab Terminal" width="750" />

---

Optionally, you can connect directly from the terminal of your computer using `ssh`. If you feel comfortable doing this, follow the instructions below. 

<details>
  <summary>
**Optional connection option** (click here)
  </summary>

Using `ssh` command to connect to the Calcul Quebec Secure Cloud. 

```{bash}
ssh user01@login1.cbw-oct-2020.calculquebec.cloud

```

You will be in your home folder. At this step, before continuing, request an interactive compute node for the workshop using the following command: 

```{bash}
salloc -N 1 -n 2  --mem 7000

```

</details>


---
<a name="map_tutorial"></a> 
## 2. Mapping Tutorial
<a name="getting_started"></a>
### 2.1 Getting Started

#####  Connect to the Calcul Quebec Secure Cloud

Now that you have a terminal open in the Secure Cloud, please make sure that you have access to the required modules for this workshop. To do so, copy and paste the following commands into the terminal window: 

```{bash}
source ~/.bashrc
module avail

```

You should see a list of modules display in your command line. If you don't, please notify the instructor so they can add required lines to your `.bashrc` file, then rerun the `source ~/.bashrc` command. 


##### Prepare directory for the workshop

```{bash}
rm -rf ~/WGBSworkshop
mkdir -p ~/WGBSworkshop
cd ~/WGBSworkshop

```

##### Save genome location
```{bash}
GENOME=/cvmfs/soft.mugqic/CentOS6/genomes/C3G_workshop/Homo_sapiens.GRCh38.chr19
export GENOME

```

This will define a variable `$GENOME` that will simplify future commands. 

---

#### Copy data for workshop

```{bash}
cp -r ~/projects/def-sponsor00/SHARE/WGBS/* .

```

##### Check the files

Type the following command: `ls data`. What do you see? 

<details>
  <summary>
**Solution** (click here)
  </summary>
  
You should see something similar to this

```{bash}
[user01@login1 WGBSworkshop]$ ls data
WGBS.A34002.137160.chr19.1.fastq.gz
WGBS.A34002.137160.chr19.2.fastq.gz
WGBS.A34002.137487.chr19.1.fastq.gz
WGBS.A34002.137487.chr19.2.fastq.gz
WGBS.A34002.137488.chr19.1.fastq.gz
WGBS.A34002.137488.chr19.2.fastq.gz

```


These are the files that will be used for the workshop. They contain a subset of reads from CEMT sample `CEMT0007`, which is a mammary gland epithelial cell line. You should also be able to see them in the left side navigation pane. Additionally, a copy of this tutorial was also copied to the `WGBSworkshop` folder. You can open it by double clicking it on the navigation pane if it is more convenient for you. 



**What do the ".1" and ".2" in the file names mean?**
<details>
  <summary>
**Solution** (click here)
  </summary>
  
They represent the `read1` and `read2` of the paired end reads. 

</details> 
</details> 

---
 <a name="map_bismark"></a>
### 2.2 Map Data using Bismark

We will now process and map the reads using Bismark.

#### Map the first dataset using Bismark
We will process the datasets one at a time. To align the first dataset, 


```{bash}
module load mugqic/bismark/0.16.1 mugqic/bowtie2/2.2.4 mugqic/samtools/1.3 

bismark --multicore 2 --bowtie2 $GENOME/genome/bismark_index \
    -1 data/WGBS.A34002.137160.chr19.1.fastq.gz -2 data/WGBS.A34002.137160.chr19.2.fastq.gz
    
```

What do all the options in the command mean? 

<details>
  <summary>
**Solution** (click here)
  </summary>
  
The `--multicore 2` is to do multithreaded processing to improve speed.

The `--bowtie2` is to use the mapping algorithm from bowtie2.

The `$GENOME/genome/bismark_index` specifies the location of reference genome to use.

The `-1 data/WGBS.A34002.137160.chr19.1.fastq.gz` specifies the location of read 1. Idem for `-2` which specifies read 2.


For more details, please refer to the Bismark [user guide](http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf). 
</details>


---

This step will take a few minutes to run for this reduced dataset. A dataset spanning a full genome will take several hours. 
For your own datasets, make sure you have enough computing walltime to run the alignment. 


**While you wait for the results, ask any questions you have up to this point to the instructors.**


#### Check files
At the end of the alignment, you should have the following files saved into your workshop folder: 

```
WGBS.A34002.137160.chr19.1_bismark_bt2_pe.bam  WGBS.A34002.137160.chr19.1_bismark_bt2_PE_report.txt
```

Let's look at the report: 

```{bash}
less WGBS.A34002.137160.chr19.1_bismark_bt2_PE_report.txt

```

Close the report by pressing `q`. 

---

#### Prepare files for loading in IGV

We need to sort the `bam` file and prepare an index so we will be able to load it in IGV. We will use the program `samtools` for this.

```{bash}
samtools sort WGBS.A34002.137160.chr19.1_bismark_bt2_pe.bam -o WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam 
samtools index WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam

```

#### Check files
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

**How would you repeat the alignment with the other datasets?**

<details>
  <summary>
**Solution** (click here)
  </summary>

This is the command to run bismark on the two other samples: 

```{bash}
module load mugqic/bismark/0.16.1 mugqic/bowtie2/2.2.4 mugqic/samtools/1.3 

bismark --multicore 2 --bowtie2 $GENOME/genome/bismark_index \
    -1 data/WGBS.A34002.137487.chr19.1.fastq.gz -2 data/WGBS.A34002.137487.chr19.2.fastq.gz
    
bismark --multicore 2 --bowtie2 $GENOME/genome/bismark_index \
    -1 data/WGBS.A34002.137488.chr19.1.fastq.gz -2 data/WGBS.A34002.137488.chr19.2.fastq.gz
    
```

This is the command to prepare the samples for IGV (sort and index): 

```{bash}
samtools sort  WGBS.A34002.137487.chr19.1_bismark_bt2_pe.bam -o WGBS.A34002.137487.chr19.1_bismark_bt2_pe_sorted.bam
samtools index WGBS.A34002.137487.chr19.1_bismark_bt2_pe_sorted.bam

samtools sort  WGBS.A34002.137488.chr19.1_bismark_bt2_pe.bam -o WGBS.A34002.137488.chr19.1_bismark_bt2_pe_sorted.bam
samtools index WGBS.A34002.137488.chr19.1_bismark_bt2_pe_sorted.bam

```


*Did you need to repeat the module load commands?*

*And what context would you need to repeat them?*

</details>

---
<a name="load_igv"></a>
### 2.4 Load Data and Explore using IGV

While you wait for the previous steps to finish executing, it is a good idea to begin exploring the alignments. 

#### Copy files to your local computer to view in IGV

Retrieve the files called `WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam` and `WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam.bai` from the server using the navigation pane and selecting `Download` from the right-click option menu. 

Optionally, you can retrieve the files directly using the terminal and `scp` using the commands below. 

<details>
  <summary>
**Optional download option** (click here)
  </summary>
  
  
```{bash}
scp user01@login1.cbw-oct-2020.calculquebec.cloud:/home/user01/WGBSworkshop/WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam .
scp user01@login1.cbw-oct-2020.calculquebec.cloud:/home/user01/WGBSworkshop/WGBS.A34002.137160.chr19.1_bismark_bt2_pe_sorted.bam.bai .

```

_**Important note:**_

Remember to substitute `user01` with your own username in the commands above. *Look carefully!* It appears twice: _`user01@login1...`_ and _`/home/user01`_; **both** need to be changed. 

</details>

#### Launch IGV on your computer

If you haven't installed it yet, please get it here [IGV download](http://software.broadinstitute.org/software/igv/download).

Make sure that the human genome is selected in the top left corner. It should read: **Human (hg38)**.

Load your sorted bam and index file in IGV using `File -> Load from file`.

Now go to:

```
chr19:43,375,889-45,912,052
```

And zoom in until you see something.

For instance, try the following window:

```
chr19:44,527,387-44,536,873
```

You should see something like this:

<img src="https://github.com/bioinformatics-ca/Epigenomics_2020/blob/master/region2.png?raw=true" alt="Region" width="750" />

*If it looks different, can you change the way the colors are displayed?*

*Which section of which chromosome is covered by this dataset?* 

*Can you see any interesting patterns in the coverage?* 

---
<a name="methylation_profiles"></a>
### 2.5 Generate Methylation Profiles

So far we have only mapped the reads using Bismark. We can generate methylation profiles using the following command:

```{bash}
bismark_methylation_extractor --bedGraph WGBS.A34002.137160.chr19.1_bismark_bt2_pe.bam

```

**How would you do the same for the other replicates?**

<details>
  <summary>
  **Solution:**
  </summary>

These are the commands that you should use:

```{bash}
bismark_methylation_extractor --bedGraph WGBS.A34002.137487.chr19.1_bismark_bt2_pe.bam
bismark_methylation_extractor --bedGraph WGBS.A34002.137488.chr19.1_bismark_bt2_pe.bam

```

</details>

---

Download *all* the files produced so far to your local computer using the navigation pane (or `scp`) as explained above. 

**While you wait for all the steps and downloads to finish, you can ask the instructors any questions you might have up until this point.**

---

Load all the downloaded files in IGV using `File -> Load from file`.

At this point, if you load the region `chr19:44,527,387-44,536,873` you should see something like

<img src="https://github.com/bioinformatics-ca/Epigenomics_2020/blob/master/region2_full.png?raw=true" alt="Region" width="750" />

This promoter looks to be hypomethylated. 

*Can you find a promoter that is hypermethylated?*

How about `chr19:45,637,715-45,657,380`?

*How would you look for a CpG island using this view of the data?*

---
<a name="methylkit"></a>
## 3. Differential Methylation Analysis in MethylKit

The following section will use the Bioconductor package `methylKit` to do a differential methylation analysis. 
You can do it in your own computer (if you have installed `R` and `methylKit`) or in the Secure Cloud. 

To install `methylKit` locally on your computer, make sure you have a recent version of `R` and 
follow the instructions in this [page](https://bioconductor.org/packages/release/bioc/html/methylKit.html). 

---
<a name="load_r"></a>
### 3.1 Load R and MethylKit

If you are working from the Secure Cloud, you will need to load a version of `R` that has the libraries we need. 
To do that, run the following commands: 

```{bash}
module purge
module load mugqic/R_Bioconductor/3.6.1_3.10

```

Then simply launch `R` by typing the following to your terminal: 

```{bash}
R

```

If you did this properly, the following message will be displayed and your prompt will change from `$` to `>`: 

```
R version 3.6.1 (2019-07-05) -- "Action of the Toes"
Copyright (C) 2019 The R Foundation for Statistical Computing
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

Now that all the samples have been read with `methylKit`, you can create a file list to make it easier to load the full dataset as a `methylkit` object. For the purposes of this tutorial, we will consider that samples belong to two experimental groups: `A34002.137160` as the control group (treatment = 0) and `A34002.137487` & `A34002.137488` as the treatment group (treatment = 1). We use the `methRead()` funciton to create our object, as shown below: 

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

What do all the options in the `methRead()` command mean? 

<details>
  <summary>
**Solution** (click here)
  </summary>
  
- `file.list` object points to the location of the input data in a MethylKit format. 

- `sample.id` points to a list with the appropriate sample name for each file. 

- `assembly` specifies which build of the human reference genome is used. 

- `treatment` specifies which sample belongs to each experimental group. 

- `context` specifies the methylation context. 

- `mincov` specifies the minimum coverage required to be included in the object. 


For more details, please refer to the MethylKit [user guide](https://bioconductor.org/packages/release/bioc/manuals/methylKit/man/methylKit.pdf) .



</details>

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

To check the ouptut, just type `myDiff` and read the summary. If you want an example of the output, check the solution below. 

<details>
  <summary>
  **Solution** (click here)
  </summary>

This is what the output looks like: 

```
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

</details>

---

To filter results by their statistical significance, `methylKit` provides the `getMethylDiff()` function which allows you to extract only the diferentially methylated CpG that meet a specific Q-value threshold. Additionally, it is also possible to specify whether to keep `hypo` or `hyper` methylated CpGs only. Finally, the `bedgraph()` function allows you to save the the `methylDiff` object into a BedGraph file so you can open it with your genome browser of choice. Let's create two BedGraph files with hypo and hyper methylated CpGs with a Q-value below 0.05 based on the data above: 

```{r}
myDiff.hyper = getMethylDiff(myDiff,qvalue=0.05,difference=10,type="hyper")
bedgraph(myDiff.hyper, file.name = "hyper.CpG.bedGraph", col.name = "qvalue")

myDiff.hypo = getMethylDiff(myDiff,qvalue=0.05,difference=10,type="hypo")
bedgraph(myDiff.hypo, file.name = "hypo.CpG.bedGraph", col.name = "qvalue")

```

Two new files should appear now in your workshop folder: 

```
WGBSworkshop/hyper.CpG.bedGraph
WGBSworkshop/hypo.CpG.bedGraph
```

#### Bin Results to Obtain Differentially Methylated Regions 

By default, `methylKit` will compute results with a CpG resolution. To get Differentially Methylated Regions, you have to bin your results first using a window size of your choice. The function to do this is `tileMethylCounts()`, which takes a regular methylkit object as input.  In this case, we will create 1000bp bins using the following command: 

```{r}
tiles = tileMethylCounts(myobj,win.size=1000,step.size=1000,cov.bases = 10)

```

As with CpG level results, samples need to be merged before the analysis can continue: 

```{r}
meth.tiles = unite(tiles, destrand=TRUE) 

```

Now, we will use the `calculateDiffMeth()` and `getMethylDiff()` functions to get the DMRs. 


**Do you know how to do it, based on the information above?**

**Based on the number of diferentially methylated CpGs you found above, do you anticipate many statistically significant DMRs in your analysis?**

<details>
  <summary>
  **Solution** (click here)
  </summary> 
  
Use the following commands to perform a `DMR` analysis: 

```{r}
myDiff.tiles = calculateDiffMeth(meth.tiles)

myDiff.tiles.hyper = getMethylDiff(myDiff.tiles,qvalue=0.1,difference=10,type="hyper")
bedgraph(myDiff.tiles.hyper, file.name = "hyper.DMR.bedGraph", col.name = "qvalue")

myDiff.tiles.hypo = getMethylDiff(myDiff.tiles,qvalue=0.1,difference=10,type="hypo")
bedgraph(myDiff.tiles.hypo, file.name = "hypo.DMR.bedGraph", col.name = "qvalue")


```

</details> 

---

Using the navigation pane, download the bedGraph files you just produced and try to open them with IGV. 

**Do the statistical results match what you had seen before when exploring the data?** 

**What interesting genomic features are found close to the DMRs? What could this mean?** 


### Congrats, you're done!

You can quit R using the `quit()` command. To close your connection to the Secure Cloud, close the terminal tab, and browser window. 


Once you are finished make sure you download all the files you need and continue exploring on IGV. 



