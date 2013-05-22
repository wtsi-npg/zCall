#! /usr/bin/python

# Copyright (c) 2013 Genome Research Ltd. All rights reserved.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Publication of the zCall algorithm:
# Goldstein JI, Crenshaw A, Carey J, Grant GB, Maguire J, Fromer M, 
# O'Dushlaine C, Moran JL, Chambert K, Stevens C; Swedish Schizophrenia 
# Consortium; ARRA Autism Sequencing Consortium, Sklar P, Hultman CM, 
# Purcell S, McCarroll SA, Sullivan PF, Daly MJ, Neale BM. 
# zCall: a rare variant caller for array-based genotyping: Genetics and 
# population analysis. Bioinformatics. 2012 Oct 1;28(19):2543-2545. 
# Epub 2012 Jul 27. PubMed PMID: 22843986.

# Author: Iain Bancarz, ib5@sanger.ac.uk


"""Write Plink .ped file from calls in GTC files; use to generate test data

Produces a single line, representing one sample, and appends to .ped file

Adapted from legacy zCall.py
Usage: writePEDfile.py [bpm path] [gtc path] [ped path]"""

import sys
sys.path.append(sys.path[0]+'/../zcall')
from BPM import *
from GTC import *

def main():

    bpm = BPM(sys.argv[1])
    gtc = GTC(sys.argv[2], bpm.normID)
    outPath = sys.argv[3]
    useManifest = False
    
    ### Parse sample name from input gtc file name
    sampleName = sys.argv[2].split("/")
    sampleName = sampleName[len(sampleName) - 1]
    sampleName = sampleName.split(".")[0]

    out = [sampleName, sampleName, "0", "0", "-9", "-9"] #output holder in python list; have no sample information so use "0" and "-9" for mid, pid, gender, case/control status
    for i in range(gtc.getTotalSNPs()):
        if useManifest: 
            alleleA = bpm.A[i]
            alleleB = bpm.B[i]
        else:
            alleleA = "A"
            alleleB = "B"
        origCall = gtc.genotypes[i]     
        if origCall == 1: ## AA is original call
            out.append(alleleA)
            out.append(alleleA)
        elif origCall == 2: ## AB is original call
            out.append(alleleA)
            out.append(alleleB)
        elif origCall == 3: ## BB is original call
            out.append(alleleB)
            out.append(alleleB)
        else: ## NC is original call
            out.append("0")
            out.append("0") 
    ## Output to std out the new calls in PED format
    outFile = open(outPath, 'a') # append to existing file
    outFile.write(" ".join(out)+"\n")
    outFile.close()

if __name__ == "__main__":
    main()
