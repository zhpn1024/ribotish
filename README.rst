README for Ribo-TISH (0.1.5)
==================================
<2017-3-15 Peng Zhang>

Introduction
============

Translation is a critical step in gene regulation that synthesizes proteins from a given RNA template. The development of the ribosome profiling (riboseq) technique has enables the measurement of translation at a genome-wide level. The basic idea of ribosome profiling is to perform deep-sequencing of the ribosome-protected mRNA fragment (~30 nts), termed ribosome footprints to determine the occupancy of translating ribosomes on a given mRNA. There are several variants of the ribosome profiling technique that are based on the use of different translation inhibitors. The regular ribo-seq utilizes Cycloheximide (CHX), a translation elongation inhibitor to freeze all translating ribosomes. In contrast to CHX, the translation inhibitor lactimidomycin (LTM) and harringtonine (Harr) have a much stronger effect on initiating ribosomes. The use of these two inhibitors allows for the global mapping of translating initiating sites (TISs) when they are coupled with with ribosome profiling (TI-Seq). In addition, when LTM is used sequentially with puromycin (PMY), the TISs can be mapped quantitatively and can be compared between different conditions.
we present a novel algorithm, named Ribo TIS Hunter (Ribo-TISH), for identifying translation activities using ribosome profiling data. Ribo-TISH uses statistical tests to assess the significance of translation activities. It captures significant TISs using negative binomial test, and frame biased open reading frames (ORFs) using rank sum test. Ribo-TISH can also perform differential analysis between two TI-Seq data.

Install
=======

Please check the file 'INSTALL.rst' in the distribution.

Usage of Ribo-TISH
========================

::

  ribotish [-h] [--version] {quality,predict,tisdiff}

:Example for quality control: ``ribotish quality -b ltm.bam -g gene.gtf -t -o ltm_qual.txt``

:Example for prediction: ``ribotish predict -t ltm.bam -b chx.bam -g gene.gtf -f genome.fa -o pred.txt``

:Example for differential TIS: ``ribotish tisdiff -1 pred1.txt -2 pred2.txt -a qti1.bam -b qti2.bam -g gene.gtf -o diff.txt --plotout diff.pdf``

There are 3 functions available as sub-commands.

:quality:	Quality control for riboseq bam data.
:predict:	Main function to predict ORF/TIS.
:tisdiff:	Call diffential TIS between two TIS data

The main input data should be in bam file format. Reads should be trimmed and aligned to genome. Intron splicing is supported. Some attributes are needed such as NM, NH and MD. For STAR, ```--outSAMattributes All``` should be set. bam file should be sorted and indexed by samtools_.

.. _samtools: https://github.com/samtools/samtools


quality
~~~~~~~

Quality control of riboseq bam data. This function checks reads distribution around annotated protein coding regions on user provided transcripts, show frame bias and estimate P-site offset for different group of reads. Reads are grouped by read length as well as 5' end match or mismatch. 5' end mismatch ('m0') reads often have different distribution from matched reads. To turn off 5' end mismatch grouping, use ```--nom0```. 

There are 3 output files: a txt file recording all distribution data, a pdf figure file and a python file for P-site offset parameters. 

Quick examples:

For regular riboseq
::

  ribotish quality -b chx.bam -g gene.gtf

For TI-Seq data
::

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
```````

Threshold for quality. Default: 0.5. Group that frame bias ratio < TH will be considered as low quality and this group of reads will not be used in further analysis. The offset for low quality groups will not be set in parameter file.

-p NUMPROC
``````````

Number of processes. Default: 1

-v/--verbose
`````````````

Increase output verbosity.


Output files
------------

OUTPUT
```````

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

Ribo-TISH does not guarantee that it can always find best P-site offset values. Users should check the quality figures and edit the parameter file if necessary. 

predict
~~~~~~~

This is the main function of Ribo-TISH. This function predicts ORF/TIS with riboseq bam files. This function uses negative binomial model to fit TI-Seq background and test significance of TIS sites. For regular riboseq data, rank sum test between in frame reads and out frame reads inside the ORF is tested.

Quick examples:

Combine TI-Seq and regular riboseq data
::

  ribotish predict -t ltm.bam -b chx.bam -g gene.gtf -f genome.fa -o pred.txt

For TI-Seq data only
::

  ribotish predict -t ltm.bam -g gene.gtf -f genome.fa -o pred.txt

User provided candidates with two regular riboseq data
::

  ribotish predict -b chx1.bam,chx2.bam -g gene.gtf -f genome.fa -i cand.txt -o pred.txt

Options
--------------

-t TISBAMPATHS
``````````````

Input TI-seq bam data files, comma seperated.

-b RIBOBAMPATHS
```````````````

Regular riboseq bam data files, comma seperated. 

At least one bam file should be provided by either ```-t``` or ```-b```.

-g GENEPATH
```````````

Gene annotation file for TIS background estimation and ORF prediction. Acceptable formats include gtf, gff, bed and genepred with gene names. Input file format can be auto detected or specified by ```--geneformat``` option. 
If user need to use different gene annotation files for background estimation and prediction, use ```-a``` option to provide another gene annotation for prediction. If user provided candidates ```-i``` option is set, the transcript annotation for the candidates should be found in gene annotation file.

-f GENOMEFAPATH
```````````````

Genome fasta file. The fasta file should has a .fai index file accompanied with genome fasta file (indexed) or indexable (fasta sequences have fixed length in each line). This program will index the genome file before prediction if .fai index file can not be found.

-o OUTPUT
`````````

Output all possible ORF results that fit the thresholds. 


-i INPUT
````````

Only test input candidate ORFs, format: 

=======  =====  =====
transID  start  stop 
=======  =====  =====

Start, stop position is 0 based, half open. Stop - start should be multiples of 3. Transcript should be found in gene annotation file.

--geneformat GENEFORMAT
```````````````````````

Gene annotation file format (gtf, bed, gpd, gff, default: auto)

--tispara TISPARA
`````````````````

Input P-site offset parameter files for ```-t``` bam files. The default parameter files are bampath+'.para.py' for each bam file, which is generated in ```ribotish quality``` function. To use this option, each bam file should be provided with a file, and file names are separated with comma. If no parameter file is found, default offset 12 will apply for all reads in the bam data.

--ribopara RIBOPARA
```````````````````

Input P-site offset parameter files for ```-b``` bam files. Same as ```--tispara``` option.

--nparts NPARTS
```````````````

Group transcript according to TIS reads density quantile. Default: 10.

TIS background estimation uses ORF in-frame read counts to estimate negative binomial parameters. Since different transcripts have different expression levels, the background is different for highly expressed and lowly expressed transcripts. Ribo-TISH groups expressed transcripts into N parts based on TIS reads density of the transcript. Each transcript group have same total number of TIS reads.

-e ESTPATH
``````````

Output TIS background estimation result. If only one bam file is provided by ```-t``` option, the default file name is tisbampath+'.bgest.txt'. If multiple TIS data provided, the default file name is tisBackground.txt
The result file contains negative binomial parameters, group levels and thresholds for each group.

-s INESTPATH
````````````

Input background estimation result file instead of instant estimation. By default, if only one bam file is provided by ```-t``` option, the program will first look for file name tisbampath+'.bgest.txt'. If this file exists, background parameters in this file will be used. Otherwise, TIS background estimation will run and generate a result file according to ```-e``` option.


-a AGENEPATH
````````````

Another gene annotation file for ORF prediction instead of ```-g``` gene file

--alt
`````

Use alternative start codons. If set, all codons with 1 base different from ATG will be considered as start codon in ORF finding. Affect both TIS background estimation and prediction. Do not affect ```-i``` mode prediction. To customize alt start codons, use ```--altcodons```.


--altcodons ALTCODONS
`````````````````````

Use provided alternative start codons, comma seperated, e.g. ```--altcodons CTG,GTG,ACG```. Turn on ```--alt``` option. Do not need to provide 'ATG'. Do not support 'N' bases.

--tis2ribo
``````````

Add TIS bam counts to regular riboseq counts. Use TIS data also for ORF frame test. This option will be turned on automatically if ```-b``` is not provided.

--harr
``````

The data is treated with harringtonine (instead of LTM). For Harr data, the reads at TIS sites are not as focus as LTM reads. Reads in flanking region (default 15 codons) of TIS will not be used for TIS background estimation. To customize flanking size, use ```--harrwidth```.


--harrwidth HARRWIDTH
`````````````````````

Flanking region for harr data, in codons. Default: 15. Turn on ```--harr``` option.

--enrichtest
````````````

Use enrich test instead of frame test. Enrich test is rank sum test between in-frame reads inside ORF and same frame reads outside ORF.

--nocompatible
``````````````

Do not require reads compatible with transcript splice junctions. 

--minaalen MINAALEN
```````````````````

Minimum amino acid length of candidate ORF, Default: 6.

--genefilter GENEFILTER
```````````````````````

Only process given genes. Comma separated. 

--tpth TPTH
```````````

TIS p value threshold. Default: 0.05.

--fpth FPTH
```````````

Frame p value threshold. Default: 0.05.

--minpth MINPTH
```````````````

At least one of TIS or frame p value should be lower than this threshold. Default: 0.05.

--fspth FSPTH
`````````````

Fisher's p value threshold. Default: 0.05.

--fsqth FSQTH
`````````````

Fisher's FDR q value threshold. Default: 1.

-p NUMPROC
``````````

Number of processes. Default: 1

-v/--verbose
`````````````

Increase output verbosity.


Output files
------------

OUTPUT
```````
The output is a txt file all possible ORF results that fit the thresholds. Some of the columns are:

:GenomePos:	Genome position and strand of TIS site, 0 based, half open
:Start:		TIS of the ORF on transcript
:stop:		3' end of stop codon on transcript
:TisType:	Relative position of this TIS to annotated ORF of the transcript. 'Novel' if no ORF annotation.
:TISGroup:	Group of the transcript for TIS background estimation
:TISCount:	Number of reads with P-site at TIS site
:TISPvalue:	One tailed negative binomial test p-value for TISCount (TIS test)
:RiboPvalue:	One tailed rank sum test p-value for regular riboseq frame bias inside ORF (frame test)
:RiboPStatus:	For all ORFs sharing same stop codon, 'T' means top (best) p-value, 'L' means local best p-value, 'N' means other. All 'N' in ```-i``` mode.
:FisherPvalue:	Combination of TIS and Ribo p-values using Fisher's method
:TISQvalue:	BH correction q-value of TIS test
:RiboQvalue:	BH correction q-value of frame test
:FisherQvalue:	BH correction q-value of Fisher's p-value
:AALen:		Amino acid length of the ORF

tisdiff
~~~~~~~

This is the function for differential TIS dentification. This function uses two different TIS test results generated by ```ribotish predict``` using different QTI-Seq data. First a normalization factor is estimated by Trimmed Mean of M values (TMM) method on common significant TIS counts in the two results. Then binomial test p-value and fold change are calculated.

Quick examples:

::

  ribotish tisdiff -1 pred1.txt -2 pred2.txt -a qti1.bam -b qti2.bam -g gene.gtf -o diff.txt --plotout diff.pdf

Options
--------------

-1 TIS1PATH, -2 TIS2PATH
````````````````````````

Predict result of group 1 & 2 TIS data

-a TIS1BAMPATHS, -b TIS1BAMPATHS
````````````````````````````````

Group 1 & 2 TIS riboseq bam files, comma seperated

-g GENEPATH
```````````

Gene annotation file. Acceptable formats include gtf, gff, bed and genepred with gene names. Input file format can be auto detected or specified by ```--geneformat``` option. 

-o OUTPUT
`````````

Output result file


--geneformat GENEFORMAT
```````````````````````

Gene annotation file format (gtf, bed, gpd, gff, default: auto)

--tis1para TIS1PARA, --tis2para TIS2PARA
````````````````````````````````````````

Input P-site offset parameter files for group 1 & 2 bam files. The default parameter files are bampath+'.para.py' for each bam file, which is generated in ```ribotish quality``` function. To use this option, each bam file should be provided with a file, and file names are separated with comma. If no parameter file is found, default offset 12 will apply for all reads in the bam data.


--nocompatible
``````````````

Do not require reads compatible with transcript splice junctions. 

--plotout PLOTOUT
`````````````````

Scatter plot output pdf file.

--figsize FIGSIZE
`````````````````

Scatter plot figure size. Default: 8,8.

-f FOLDCHANGE
`````````````

Minimum fold change threshold. Default: 1.5.

--ipth IPTH
```````````

Input TIS p value threshold. Default: 0.05.

--iqth IQTH
```````````

Input TIS q value threshold. Default: 0.05.

--opth OPTH
```````````

Output TIS diff p value threshold. Default: 0.05.

--oqth OQTH
```````````

Output TIS diff q value threshold. Default: 0.05.

-p NUMPROC
``````````

Number of processes. Default: 1

-v/--verbose
`````````````

Increase output verbosity.


Output files
------------

OUTPUT
```````
The output is a txt file all differential TIS results that fit the thresholds. Some of the columns are:

:FoldChange:	Fold change (b/a) value after normalization
:DiffPvalue:	Binomial differential test p-value, one tailed.
:DiffQvalue:	BH correction q-value of DiffPvalue
