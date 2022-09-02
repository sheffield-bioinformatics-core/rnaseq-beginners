---
title: "RNA-seq for Beginners"
author: "Sheffield Bioinformatics Core"
output:
  html_document:
    toc: yes
    toc_float: yes
    css: stylesheets/styles.css
editor_options:
  chunk_output_type: inline
---

### Sheffield Bioinformatics Core

<img src="media/logo-sm.png" align="right"/>

web : [sbc.shef.ac.uk](https://sbc.shef.ac.uk)\
twitter: [\@SheffBioinfCore](https://twitter.com/SheffBioinfCore)\
email: [bioinformatics-core\@sheffield.ac.uk](bioinformatics-core@sheffield.ac.uk)

------------------------------------------------------------------------

# Workshop Overview

An overall workflow for the processing and analysis of RNA-seq data is given in the image below from [Ting-You Wang's RNA-seq analysis page](https://databeauty.com/blog/tutorial/2016/09/13/RNA-seq-analysis.html).

![](https://databeauty.com/figures/2016-09-13-RNA-seq-analysis/rna_seq_workflow.png)

In this workshop we are going to concentrate on differential-expression and pathways analysis - which are the **least computationally-intensive** parts of the workflow and require the **least Bioinformatics experience**.

For those interested in alignment and QC steps, we have some materials available that use the Galaxy online resource

-   [Pre-processing materials](https://sbc.shef.ac.uk/ngs_intro_workshop/03-rna-seq.nb.html)

**We will also have some courses on using the command-line and the nextflow workflow manager to process raw RNA-seq data.**

However, regardless of whatever method you use to process the data, decisions that you make before commencing sequencing can have a huge impact on the results.

## Experimental Design

Before embarking on any high-throughput experiment, it is important to pay due attention to the *experimental design*. This famous quote from the statistician R.A Fisher in the 1938 is still applicable to modern technologies

> **To call in the statistician after the experiment is done may be no more than asking him to perform a postmortem examination: he may be able to say what the experiment died of.**

Experimental Design encompasses questions such as

-   which controls to use?
    -   positive / negative controls
    -   healthy controls
-   what experimental conditions?
-   what technical and biological factors are present?

When performing a high-throughput experiment, our measurements will be subject to biological variation (which we may be interested in) and technical variation (which we probably won't be). Being able to control these factors and minimise biases is key to experimental design.

Confounding factors in our design may arise by accident, or might be caused by not considering all possible sources of variation:-

<img src="media/confounding_factor.PNG"/> (image from Cancer Research UK Cambridge Institute course on Experimental Design)

We can also introduce so-called "batch-effects" by our choice of when samples are prepared for sequencing. Large experiments may necessitate multiple runs or batches, and we should try and minimize the possible impact of batches but including a good representation of each condition of interest in each batch.

<img src="media/blocking.PNG"/> (image from Cancer Research UK Cambridge Institute course on Experimental Design)

When planning next-generation sequencing experiments, you will also need to consider

-   the type of sequencing (e.g. whole-genome, exome, RNA-seq)
    -   will largely be dictated by your biological question
-   single-end or paired-end
-   how many reads (10 Million? 20 Million?, 100 Million?)

Some recommendations on these questions and more are provided by the [Cancer Research Uk Cambridge Institute Genomics Core](https://www.cruk.cam.ac.uk/core-facilities/genomics-core/sequencing). Often the sequencing vendor performing your experiment will have some default options available.

The vendor may not advise on the *sample-size*; how many samples you will be sequencing to address your biological hypothesis of interest. This is a complex question and is often influenced by practical and financial constraints. The Sheffield Bioinformatics Core is able to advise on this, and any of the other issues above. `bioinformatics-core@sheffield.ac.uk`

# Differential expression

## Input Data

A differential expression analysis requires two input files to be created.

-   a count matrix
-   a sample information table

The count matrix can be obtained by performing quantification (outside the scope of this workshop...). This will usually be generated for you by the sequencing vendor or Bioinformatics Core. The structure of the count matrix is shown below.

| **gene** | **sampleA** | **sampleB** |
|:--------:|:-----------:|:-----------:|
|    A     |    1500     |     900     |
|    B     |     20      |     10      |
| **...**  |   **...**   |   **...**   |

The gene named A was sequenced 1500 times on sampleA and 900 times on sampleB etc. These are referred to as *raw counts*, and cannot just put this numbers into a standard statistical test (e.g. t-test) to assess significance. There are several reasons for this.

-   The count values do not follow a normal-distribution so cannot be analysed using traditional methods
-   There are many, many genes being measured in the dataset leading to a multiple testing problem.
-   The count values are influenced by technical (as opposed to biological) variation that need to be accounted for:-
    -   size of gene; *longer* genes will have more reads assigned to them
    -   library size; for a sample that is sequenced to a higher depth it will seem as though all genes are more highly-expressed.

The term *differential expression* was first used to refer to the process of finding statistically significant genes from a *microarray* gene expression study.

![](media/de_explained.png)

Such methods were developed on the premise that microarray expression values are approximately *normally-distributed* when appropriately transformed (e.g. by using a log$_2$ transformation) so that a modified version of the standard *t-test* can be used. The same test is applied to each gene under investigation yielding a *test statistic*, *fold-change* and *p-value*. Similar methods have been adapted to RNA-seq data to account for the fact that the data are *count-based* and do not follow a normal distribution.

## Interactive exploration of the results with *DEGUST*

![](media/rna_advanced_degust_1.png)

<font size="8"><http://degust.erc.monash.edu/></font>

`Degust` is a web tool that can analyse counts files and test for differential gene expression. It offers and interactive view of the differential expression results and also sample quality assessment.

R-based methods such as `edgeR` (implemented in Degust) and `DESeq2` have their own method of normalising counts. You will probably encounter other methods of normalising RNA-seq reads such as *RPKM*, *CPM*, *TPM* etc. [This blog](https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/) provides a nice explanation of the current thinking. As part of the `Degust` output, you have the option of downloading normalised counts in various formats. Some other online visualisation tools require normalised counts as input, so it is good to have these to-hand.

We will use a previously-published count matrix. This was downloaded from the Gene Expression Omnibus (GEO) under the accession number [GSE60450](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450). Note that we have shortened the column headings and added gene symbols to help with visualisation and annotation

::: exercise
Download the counts from [this link](GSE60450_Lactation-GenewiseCounts_rename_symbol.csv)
:::

### Uploading the count matrix to Degust

-   From the main degust page, click *Upload your counts file*
-   Click on Browse
-   Select the location of the file `GSE60450_Lactation-GenewiseCounts_rename_symbol.csv`, and click *Open*.
-   Click *Upload*
-   A Configuration page will appear.

![](media/degust_config.png)

-   For Name type "*GSE60450*" (or whatever you want to call the analysis)
-   For Info columns select *SYMBOL*
-   Click Add condition
    -   Referring to the experiment design (below), select the Basal samples and call the condition Basal
    -   Repeat for the Luminal samples
-   Save the settings and then View the results

| Run        | Name    | CellType | Status    |
|------------|---------|----------|-----------|
| SRR1552444 | MCL1-LA | basal    | virgin    |
| SRR1552445 | MCL1-LB | luminal  | virgin    |
| SRR1552446 | MCL1-LC | Luminal  | pregnancy |
| SRR1552447 | MCL1-LD | Luminal  | pregnancy |
| SRR1552448 | MCL1-LE | luminal  | lactation |
| SRR1552449 | MCL1-LF | luminal  | lactation |
| SRR1552450 | MCL1-DG | basal    | virgin    |
| SRR1552451 | MCL1-DH | luminal  | virgin    |
| SRR1552452 | MCL1-DI | basal    | pregnancy |
| SRR1552453 | MCL1-DJ | basal    | pregnancy |
| SRR1552454 | MCL1-DK | basal    | lactation |
| SRR1552455 | MCL1-DL | basal    | lactation |

### Overview of Degust sections

-   Top black panel with Configure settings at right.
-   Left: Conditions: Control and Treatment.
-   Left: Method selection for DGE. 
-   Top centre: Plots, with options at right.
-   When either of the expression plots are selected, a heatmap appears below.
-   A table of genes (or features); expression in treatment relative to control (Treatment column); and significance (FDR column).

(**Not that the screenshots are for illustration purposes and taken from a different dataset to that being analysed in the tutorial**)

![](http://sepsis-omics.github.io/tutorials/modules/dge/images/image12.png)

### MDS plot

This is a multidimensional scaling plot which represents the variation between samples. It is a similar concept to a Principal Components Analysis (PCA) plot. The x-axis is the dimension with the highest magnitude. In a standard control/treatment setup, samples should be split along this axis. A desirable plot is shown below:-

![](media/degust_mds.png)

## MA-plot

![](media/degust_ma.png)

Each dot shows the change in expression in one gene.

-   The average expression (over both condition and treatment samples) is represented on the x-axis.
    -   Plot points should be symmetrical around the x-axis.
    -   We can see that many genes are expressed at a low level, and some are highly expressed.
-   The fold change is represented on the y axis.
    -   If expression is significantly different between batch and chem, the dots are red. If not, they are blue. (In Degust, significant means FDR \<0.05).
    -   At low levels of gene expression (low values of the x axis), fold changes are less likely to be significant.

Click on the dot to see the gene name.

## Parallel coordinates and heatmap

![](media/degust_parallel_heatmap.png)

Each line shows the change in expression in one gene, between control and treatment.

-   Go to Options at the right.
    -   For FDR cut-off set at 0.001.
    -   This is a significance level (an adjusted p value). We will set it quite low in this example, to ensure we only examine key differences.
-   Look at the Parallel Coordinates plot. There are two axes:
    -   Left: Control: Gene expression in the control samples. All values are set at zero.
    -   Right: Treatment Gene expression in the treatment samples, relative to expression in the control.
-   The blocks of blue and red underneath the plot are called a heatmap.
    -   Each block is a gene. Click on a block to see its line in the plot above.
    -   Look at the row for the chem. Relative to batch, genes expressed more are red; genes expressed less are blue.

## Table of genes

![](media/degust_gene_table.png)

Table of genes

-   gene_id: names of genes. Note that gene names are sometimes specific to a species, or they may be only named as a locus ID (a chromosomal location specified in the genome annotation).
-   FDR: False Discovery Rate. This is an adjusted p value to show the significance of the difference in gene expression between two conditions. Click on column headings to sort. By default, this table is sorted by FDR.
-   basal and luminal: log2(Fold Change) of gene expression. The default display is of fold change in the treatment relative to the control. Therefore, values in the batch column are zero. This can be changed in the Options panel at the top right.
    -   In some cases, a large fold change will be meaningful but in others, even a small fold change can be important biologically.

The table can be sorted according to any of the columns (e.g. fold-change or p-value)

## Download and R code

Above the genes table is the option to download the results of the current analysis to a csv file. You can also download the *R* code required to reproduce the analysis by clicking the *Show R code* box underneath the Options box.

Plots such as the MDS, MA and heatmap can also be exported by right-clicking on the plot.

## Exercise

::: exercise
**Question:** Do the sample groupings in the MDS plot make sense? Do any samples appear to be mislabeled? What effect might this have on the analysis?
:::

::: exercise
**Question:** Having identified the problem with the analysis, modify the configuration and repeat. How many genes are differentially expressed this time?
:::

## Analysing a different contrast

Comparing Basal vs Luminal wasn't really the main question of interest in the dataset, but it serves to illustrate the importance of checking QC plots.

-   Create conditions *Basal.Pregnant*, *Basal.Lactation*, etc using the corrected experimental design
-   Make sure that *Basal.Pregnant* and *Basal.Lactation* are both ticked as initial select

![](media/degust-correct-config.png)

Take some time to understand the various parts of the report

::: exercise
**Exercise:** Make sure the FDR cut-off and abs LogFC cutoffs are set to default and *download* the file and rename to `background.csv`. We will use this later.
:::

::: exercise
**Exercise**: How many genes are differentially-expressed with an FDR \< 0.05 and abs logFC \> 1. Download this file and rename it to `B.preg_vs_lactation.csv`.
:::

::: exercise
**Exercise**: Repeat the analysis for Luminal.Pregnant vs Luminal.Lactation and download the table of differentially-expressed results (same FDR and log fold-change).
:::

### File Downloads

::: information
If you didn't manage to complete these analyses, you can download the files from here by right-clicking on each link and selecting "Save Link as" (or equivalent). They are also available in the course google drive.

-   [B.preg_vs_lactation.csv](B.preg_vs_lactation.csv)
-   [L.preg_vs_lactation.csv](L.preg_vs_lactation.csv)
-   [background.csv](background.csv)
:::

## Overlapping Gene Lists

We might sometimes want to compare the lists of genes that we identify using different methods, or genes identified from more than one contrast. In our example dataset we can compare the genes in the contrast of pregnant vs luminal in basal and luminal cells

The website *venny* provides a really nice interface for doing this.

![](media/venny-config.png)

-   Open both your *Basal Pregnant vs Basal Lactation* and *Luminal Pregnant vs Luminal Lactation* results files in Excel
-   Go to the venny website
    -   <http://bioinfogp.cnb.csic.es/tools/venny/>
-   Copy the names of genes with adjusted p-value less than 0.05 in the Basal analysis into the **List 1** box on the venny website. **List 1** can be renamed to *Basal*
    -   *You can select all entries in a column with the shortcut CTRL + SPACE*
-   Copy the names of genes with adjusted p-value less than 0.05 in the Luminal analysis into the **List 2** box on the venny website. **List 2** can be renamed to **Luminal**
-   venny should now report the number of genes found in each list, the size of the intersection, and genes unique to each method

## Modified analysis using Hidden Factors

Let's consider the situation where we want to identify genes that are different between pregnancy and lactating samples *regardless of the cell type*. We have already done this using the venn diagram approach above, but in this final analysis we will include all the pregnancy and lactating samples, but correct for the differences in cell type. This is more efficient than analysing each cell type separately and comparing the results. Each statistical test we do will involve calling many false positives. We can achieve the same outcome by performing just one statistical test.

This can be done by telling Degust about the *hidden factors* in our dataset. The hidden factor in this dataset is whether the sample is from the `basal` or `luminal` samples. In other words, this is a technical factor that influences our results but not a factor that we wish to compare. We only need to specify which samples are from `basal` and DEGUST will infer that the other samples belong to a different cell type.

See below for the correct configuration to include the hidden factors.

![](media/hidden_factor.PNG)

You should see that on the MDS plot the samples cluster according to cell type. However, this is fine because we are going to incorporate this hidden factor in the analysis

![](media/hidden_factor_batches.PNG)

::: exercise
**Exercise**: How many genes are differentially-expressed with an adjusted p-value cut-off of 0.05 and log2 cutoff of 1. How does this compare the number of intersecting genes in your venn diagram above.
:::

The hidden factor method can be used to analyse datasets where samples cluster on technical rather than biological factors. e.g.

-   sample batch
-   gender

![](media/batch_effect.png)

However, this is only true if your **experimental design was correct** and the technical variation is not confounded with biological groups. e.g. treated and untreated samples should be in all batches

# Enrichment and Pathways Analysis

In this section we will use the following files

-   [`background.csv`](background.csv) containing one row for each gene in the comparison Basal.pregnant vs Basal.lactation (27,179 rows).
-   [`B.preg_vs_lactation.csv`](B.preg_vs_lactation.csv) containing one row for each found to be DE in the contrast Basal.pregnant vs Basal.Lactation.

It will be helpful to have both these files open in Excel.

In this section we move towards discovering if our results are ***biologically significant***. Are the genes that we have picked statistical flukes, or are there some commonalities.

There are two different approaches one might use, and we will cover the theory behind both. The distinction is whether you are happy to use a hard (and arbitrary) threshold to identify DE genes.

## Over-representation analysis

"Threshold-based" methods require defintion of a statistical threshold to define list of genes to test (e.g. FDR \< 0.01). Then a *hypergeometric* test or *Fisher's Exact* test is generally used. These are typically used in situations where plenty of DE genes have been identified, and people often use quite relaxed criteria for identifying DE genes (e.g. raw rather than adjusted p-values or FDR value)

The question we are asking here is;

> ***"Are the number of DE genes associated with Theme X significantly greater than what we might expect by chance alone?"***

We can answer this question by knowing

-   the total number of DE genes
-   the number of genes in the gene set (pathway or process)
-   the number of genes in the gene set that are found to be DE
-   the total number of tested genes (background)

The formula for Fishers exact test is;

$$ p = \frac{\binom{a + b}{a}\binom{c +d}{c}}{\binom{n}{a +c}} = \frac{(a+b)!(c+d)!(a+c)!(b+d)!}{a!b!c!d!n!} $$

with:-

| **Differentially Expressed** | **Not Differentially Expressed** | **Total** |                        |
|:----------------:|:----------------:|:----------------:|:----------------:|
|         In Gene Set          |                a                 |     b     |         a + b          |
|       Not in Gene Set        |                c                 |     d     |         c + d          |
|          **Total**           |            **a + c**             | **b +d**  | **a + b + c + d (=n)** |

In this first test, our genes will be grouped together according to their Gene Ontology (GO) terms:- <http://www.geneontology.org/>

## Using GOrilla

There are several popular online tools for performing enrichment analysis

We will be using the online tool [GOrilla](http://cbl-gorilla.cs.technion.ac.il/) to perform the pathways analysis. It has two modes; the first of which accepts a list of *background* and *target* genes.

1.  Go to <http://cbl-gorilla.cs.technion.ac.il/>
2.  Read the "Running Example"

![](media/gorilla-example.png)

3.  Choose Organism: `Mus Musculus`
4.  Choose running mode: `Two unranked lists of genes`
5.  Paste the gene symbols corresponding to DE genes in *Basal pregant vs Basal Lactation* into the Target set.

-   **The shortcut CTRL + SPACE will let you select an entire column**

6.  Paste the gene symbols from the Background set into the other box.
7.  Choose an Ontology: `Process`
8. Under **Advanced parameters** you can tick the option to *Output results in Microsoft Excel format*
9.  `Search Enriched GO terms`

You should be presented with a graph of enriched GO terms showing the relationship between the terms. Each GO term is coloured according to its statistical significance.

![](media/GOrilla-network.PNG)

Below the figure is the results table. This links to more information about each GO term, and lists each gene in the category that was found in your list. The enrichment column gives 4 numbers that are used to determine enrichment (similar to the Fisher exact test we saw earlier)

-   N, total number of genes (should be the same in all rows)
-   B, total number of genes annotated with the GO term
-   n, total number of genes that were found in the list you uploaded (same for all rows)
-   b, number of genes in the list you uploaded that intersect with this GO term

![](media/GOrilla-table.PNG)

::: exercise
**Exercise:** Use GOrilla to find enriched pathways in the Basal pregnant vs lactation analysis
:::



## Threshold-free analysis

This type of analysis is popular for datasets where differential expression analysis does not reveal many genes that are differentially-expressed on their own. Instead, it seeks to identify genes that as a group have a tendancy to be near the extremes of the log-fold changes. The results are typically presented in the following way.

![](media/overexpressed-gsea.png)

The "barcode"-like panel represents where genes from a particular pathway (**HALLMARK_E2F_TARGETS** in this case) are located in a gene list *ranked* from most up-regulated to most down-regulated. The peak in the green curve is used to indicate where the majority of genes are located. If this is shifted to the left or the right it indicates that genes belonging to this gene set have a tendancy to be up- or down-regulated.

As such, it does not rely on having to impose arbitrary cut-offs on the data. Instead, we need to provide a measure of the importance of each gene such as it's fold-change. These are then used the rank the genes.

The Broad institute has made this analysis method popular and provides [a version of GSEA](http://software.broadinstitute.org/gsea/index.jsp) that can be run via a java application. However, the application can be a bit fiddly to run, so we will use the GeneTrail website instead

<https://genetrail.bioinf.uni-sb.de/start.html>

-   Open the file `background.csv` in Excel and delete all columns except the `SYMBOL` and `basal.lactation` column. <img src="media/GeneTrail_prep.png"/>
-   Go to the GeneTrail website, and select Transcriptomics from the front page
-   Select the **Paste the content of a text file in a tabular format option** and the contents of your modified excel file into the box. **Do not paste the column headings**
-   Click Upload

Hopefully it should recognise your input without any errors, and on the next screen the **Set-level statistic** should be automatically set to **GSEA**

::: warning
If your data does not get uploaded, double-check that the column heading **basal.lactation** has not been pasted into the text box
:::

To make the analysis run faster, you can de-select the GO pathways (biological processes, molecular function and cellular compartment)

<img src="media/genetrail_setup.PNG"/>

After a short wait, you will be able to view and download the results. The tested pathways are grouped into different sources (Kegg, Reactome or Wikipathways)

<img src="media/genetrail_KEGG.PNG"/>

Each of the significant pathways can be explored in detail by clicking the **More..** link; such as showing which genes in that pathways are up- or downregulated.

<img src="media/genetrail_KEGG_result.PNG"/>

The Rank of the gene shown is the position of the gene in the ranked list; with 1 being most up-regulated gene. The score is the score used to rank the genes (fold-change in our example).

::: exercise
**Exercise:** Use GeneTrail to analyse
:::