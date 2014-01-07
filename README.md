
zCall: A Rare Variant Caller for Array-Based Genotyping
=======================================================

1. Overview
-----------

zCall is a variant caller specifically designed for calling rare single 
nucleotide polymorphisms (SNPs) from array-based technology. This caller is 
implemented as a post-processing step after a default calling algorithm has 
been applied such as Illumina's GenCall algorithm. zCall uses the intensity 
profile of the common allele homozygote cluster to define the location of the 
other two genotype clusters. 

The zCall code includes three prototype versions, and an extended version; 
these are documented respectively in README_prototypes and README_extended.md. 
The prototypes are run in a series of steps with some manual intervention; 
extended zCall can be run with a single command and has other added 
capabilities.

2. Publication and downloads
-----------------------------

The paper describing zCall is:
Goldstein JI, Crenshaw A, Carey J, Grant GB, Maguire J, Fromer M, 
O'Dushlaine C, Moran JL, Chambert K, Stevens C; Swedish Schizophrenia 
Consortium; ARRA Autism Sequencing Consortium, Sklar P, Hultman CM, Purcell S, 
McCarroll SA, Sullivan PF, Daly MJ, Neale BM. zCall: a rare variant caller 
for array-based genotyping: Genetics and population analysis. Bioinformatics. 
2012 Oct 1;28(19):2543-2545. Epub 2012 Jul 27. PubMed PMID: 22843986.

zCall is hosted on Github (https://github.com/jigold/zCall). Prototype versions 
are available as .zip files, while the extended version is in the src 
directory. The entire zCall repository can be cloned using Git or downloaded 
from Github as a .zip file (approximately 0.8 MB). Extended zCall can be 
installed from a download of the full repository, using the included Makefile.

3. Disclaimer
---------------

The prototype and extended editions of zCall include code provided by Illumina. 
The Illumina provided Code was provided as-is and with no warranty as to 
performance and no warranty against it infringing any other party's 
intellectual property rights.

4. Contacts
------------

The original zCall method and prototype implementations were developed by 
Jackie Goldstein et al.  For questions about prototype zCall or reporting 
problems with the code, please send an email to Jackie Goldstein 
(jigold@broadinstitute.org). For all other inquiries, please send an email to 
both Ben Neale (bneale@broadinstitute.org) and Jackie Goldstein 
(jigold@broadinstitute.org).

Any queries concerning the extended version of zCall in the src directory 
should be directed to Iain Bancarz (ib5@sanger.ac.uk).

This document was written by Iain Bancarz.