Ribo TIS Hunter
==================================
README for Ribo TIS Hunter (0.1.0)
==================================
Time-stamp: <2016-09-29 Peng Zhang>

Introduction
============

Translation is a critical step in gene regulation that synthesizes proteins from a given RNA template. The development of the ribosome profiling (riboseq) technique has enables the measurement of translation at a genome-wide level. The basic idea of ribosome profiling is to perform deep-sequencing of the ribosome-protected mRNA fragment (~30 nts), termed ribosome footprints to determine the occupancy of translating ribosomes on a given mRNA. There are several variants of the ribosome profiling technique that are based on the use of different translation inhibitors. The regular ribo-seq utilizes Cycloheximide (CHX), a translation elongation inhibitor to freeze all translating ribosomes. In contrast to CHX, the translation inhibitor lactimidomycin (LTM) and harringtonine (Harr) have a much stronger effect on initiating ribosomes. The use of these two inhibitors allows for the global mapping of translating initiating sites (TISs) when they are coupled with with ribosome profiling (TI-Seq). In addition, when LTM is used sequentially with puromycin (PMY), the TISs can be mapped quantitatively and can be compared between different conditions.
we present a novel algorithm, named Ribo TIS Hunter (ribotish), for identifying translation activities using ribosome profiling data. Ribo TIS Hunter captures significant TISs using negative binomial model, and frame biased open reading frames (ORFs) using rank sum test. Ribo TIS Hunter can also do differential analysis between two TI-Seq data.

Install
=======

Please check the file 'INSTALL' in the distribution.

Usage of Ribo TIS Hunter
========================

::

  ribotish [-h] [--version] {quality,predict,tisdiff}

:Example for quality control: ``ribotish quality -b ltm.bam -g gene.gtf -t -o ltm_qual.txt``

:Example for prediction: ``ribotish predict -t ltm.bam -b chx.bam -g gene.gtf -f genome.fa -o pred.txt``

:Example for differential TIS: ``ribotish tisdiff -1 pred1.txt -2 pred2.txt -a ltm1.bam -b ltm2.bam -g gene.gtf -o diff.txt --plotout diff.pdf``

There are 3 functions available as sub-commands.

:quality:	Quality control for riboseq bam data.
:predict:	Main function to predict ORF/TIS.
:tisdiff:	Call diffential TIS between two TIS data

The main input data should be in bam file format. Reads should be trimmed and aligned to genome. Intron splicing is supported. Some attributes are needed such as NM, NH and MD. For STAR, ```--outSAMattributes All``` should be set. bam file should be sorted and indexed by samtools.

quality
~~~~~~

Quality control of riboseq bam data. This function checks reads distribution around annotated protein coding regions on user provided transcripts, show frame bias and estimate P-site offset for different group of reads. Reads are grouped by read length as well as 5' end match or mismatch. 5' end mismatch ('m0') reads often have different distribution from matched reads. To turn off 5' end mismatch grouping, use ```--nom0```. 

There are 3 output files: a txt file recording all distribution data, a pdf figure file and a python file for P-site offset parameters. 

Quick examples:
For regular riboseq::
ribotish quality -b chx.bam -g gene.gtf
For TI-Seq data::
ribotish quality -b ltm.bam -g gene.gtf -t

Options
--------------

-b RIBOBAMPATH
``````````````

Riboseq bam data file. Reads should be trimmed and aligned to genome.
-g GENEPATH
```````````

Gene annotation file. Acceptable formats include gtf, gff, bed and genepred with gene names. Input file format can be auto detected or specified by ```--geneformat``` option


-o OUTPUT
`````````

Output all distribution data. Default: bampath[:-4]+'_qual.txt'. Quality and offset estimation is based on this distribution. User can save this file for further quick estimation trying different thresholds by ```-i``` option.

-t/--tis
````````

This data is TIS enriched, for LTM and Harr. Quality will pay more attention to TIS sites.

-i INPUT
````````

Input previous output file, do not read gene file and bam file again.

--geneformat GENEFORMAT
```````````````````````

Gene annotation file format (gtf, bed, gpd, gff, default: auto)

-f FIGPDFPATH
`````````````

Output pdf figure file. Default: bampath[:-4]+'_qual.pdf'

-r PARAPATH
```````````

Output offset parameter file. Default: bampath+'.para.py'. This file saves P-site offsets for different reads lengths in python code dict format, and can be used in further analysis.

-l LENS
```````

Range of tag length Default: 25,35. The last number (35) is not included, i.e. the longest length considered is 34.

-d DIS
``````

Position range near start codon or stop codon Default: -40,20

--bins BINS
```````````

Bins for cds profile Default: 20

--nom0
```````````

Do not consider reads with mismatch at position 0 (5' end mismatch) as a new group.

--th TH
````

Threshold for quality. Default: 0.5. Group that frame bias ratio < TH will be considered as low quality and this group of reads will not be used in further analysis. The offset for low quality groups will not be set in parameter file.

-p NUMPROC
``````````

Number of processes. Default: 1

-v, --verbose
```````````

Increase output verbosity.


Output files
~~~~~~~~~~~

OUTPUT
``````

OUTPUT is a txt file recording all distribution data in python format for each group of reads. These distributions are shown in pdf figure file. Quality and offset estimation is based on this distribution. User can save this file for further quick estimation trying different thresholds by ```-i``` option.

Pdf figure
``````````

Pdf figure file is plot of all the distributions and illustration of quality and P-site offset. The left part is for 5' end matched reads and the right part is for 5' end mismatch reads if ```--nom0``` is not set. 

Upper panel: the length distribution of RPFs uniquely mapped to annotated protein-coding regions.

Lower panel: different quality metrics for RPFs uniquely mapped to annotated protein-coding regions.
Each row shows the RPFs with different lengths.
 - Column 1: distribution of RPF 5’ end in 3 frames in all annotated codons. The percentage of the reads from the dominant reading frame is shown. 
 - Column 2: the distribution of RPF 5’end count near annotated TIS. The estimate of the P site offset and TIS accuracy are also shown. The RPFs of a specific length that do not pass threshold are considered as low quality and removed.              
 - Column 3: the distribution of RPF 5’end count near annotated stop codon. 
 - Column 4: The RPF profile throughout the protein-coding regions in 3 frames. TIS enrich score (TIS count / CDS average) is also shown for TIS data.

Offset parameter file
`````````````````````

This file saves P-site offsets for different reads lengths in python code dict format, and can be used in further analysis. The default offset file name is bampath+'.para.py' accompanied with the input bam file, and this default file name will be auto-recognized in further analysis. The offset parameter file is easy to interpret and can be edited by user if auto estimated offsets are not satisfying. If the bam file is in a different directory and user do not want to create a parameter file in that directory, we recommend creating a link for the bam file in current working directory, e.g. ```ln -s original/dir/ribo.bam```



Other useful links
==================

:Cistrome: http://cistrome.org/ap/
:bedTools: http://code.google.com/p/bedtools/
:UCSC toolkits: http://hgdownload.cse.ucsc.edu/admin/exe/


