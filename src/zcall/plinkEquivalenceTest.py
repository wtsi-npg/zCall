#! /software/bin/env python

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


"""Standalone script to compare Plink binary datasets

Test datasets for equivalence within flip of major/minor alleles
Use to validate test data for zcall

Inputs: .bpm.csv path, and two Plink filename prefixes to compare
"""

import re, struct, sys
from plink import PlinkHandler
from BPM import *

class PlinkEquivalenceTester(PlinkHandler):

    def alleleFlips(self, bimPath1, bimPath2):
        """Count allele flips (eg. AG -> GA) in two .bim files"""
        lines1 = open(bimPath1).readlines()
        lines2 = open(bimPath2).readlines()
        flips = [None]*self.snpTotal
        if self.snpTotal!=len(lines2):
            raise ValueError("SNP totals in .bim files do not match!")
        for i in range(self.snpTotal):
            (a1a, a1b) = self.parseAlleles(lines1[i])
            (a2a, a2b) = self.parseAlleles(lines2[i])
            if a1a==a2a and a1b==a2b:
                flip = 0
            elif a1a==a2b and a1b==a2a:
                flip = 1
            else:
                raise ValueError("Allele values incompatible with flip!")
            flips[i] = flip
        print len(flips), "total SNPs"
        print sum(flips), "flipped alleles in .bim files"
        print round(float(sum(flips))/len(flips), 3), "flip rate"
        return flips

    def compareFam(self, famPath1, famPath2):
        """Check number of lines in .fam files; should be one per sample"""
        len1 = len(open(famPath1).readlines())
        len2 = len(open(famPath2).readlines())
        if len1 != len2:
            raise ValueError(".fam files of unequal size!")
        print len1, "samples found."
        return len1

    def compareGenotypes(self, g1, g2):
        """Compare genotype calls in two .bed files

        A 'reversal' is 1 (major homozygote) to 3 (minor homozygote) or 3 to 1
        Any other call difference is counted as an error
        Mismatch rate = (errors+reversals)/total
        """
        callTotal = len(g1)
        if callTotal!=len(g2): 
            print callTotal, len(g2)
            raise ValueError("Total numbers of calls do not match!")
        elif callTotal % self.snpTotal != 0:
            raise ValueError("Total calls are not a multiple of total SNPs!")
        equivalent = True
        reversals = 0
        errors = 0
        for i in range(callTotal):
            j = i % self.snpTotal
            if g1[i]==g2[i]:
                continue
            elif (g1[i]==3 and g2[i]==1) or (g1[i]==1 and g2[i]==3):
                reversals += 1
            else:
                errors += 1
            equivalent = False
        print reversals, "reversed genotypes"
        print errors, "total errors"
        print callTotal, "total calls"
        print round((reversals+errors)/float(callTotal), 6), "mismatch rate"
        return equivalent

    def plinkBinaryEquivalent(self, stem1, stem2):
        """Are the given plink binary datasets equivalent?

        Arguments: Two Plink file stems (without .bed, .bim, .map suffix)
        Plink may arbitrarily flip major/minor alleles
        Read .bim files to identify flips
        Compare .bed files, allowing for flip events

        (Used to validate test data for zCall)
        """
        sampleTotal = self.compareFam(stem1+'.fam', stem2+'.fam')
        gtypes1 = self.readBedFile(stem1+'.bed', sampleTotal)
        gtypes2 = self.readBedFile(stem2+'.bed', sampleTotal)
        flips = self.alleleFlips(stem1+'.bim', stem2+'.bim')
        return self.compareGenotypes(gtypes1, gtypes2)

    def parseAlleles(self, bimLine):
        """Parse allele designations from a .bim file line"""
        words = re.split("\s+", bimLine.strip())
        return (words[4], words[5])

def main():
    if len(sys.argv)!=4:
        print "Usage: "+sys.argv[0]+" [.bpm.csv path] [stem 1] [stem 2]\n"+\
            "'stem' is path to a Plink binary dataset without .bed, .bim, .fam suffix"
        sys.exit(0)

    bpmPath = sys.argv[1]
    stem1 = sys.argv[2]
    stem2 = sys.argv[3]    
    bpm = BPM(bpmPath)
    PlinkEquivalenceTester(bpm).plinkBinaryEquivalent(stem1, stem2)

if __name__ == "__main__":
    main()
