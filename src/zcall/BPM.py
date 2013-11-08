#! /usr/bin/python
import sys


# Original code supplied by Illumina, Inc. subsequently modified by Broad 
# Institute and Genome Research Ltd. The Illumina provided code was provided 
# as-is and with no warranty as to performance and no warranty against it 
# infringing any other party's intellectual property rights. All contributions 
# are copyright their respective authors. See https://github.com/wtsi-npg/zCall
# for revision history.
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


class BPM:
    ''' Python class to parse a .bpm.csv file '''
    def __init__(self, bpmFile, normalize=True):
        self.names = []
        self.chr = []
        self.pos = []
        self.normID = []
        self.A = []
        self.B = []
        self.ilmnStrand = [] # illumina strand designation
        self.numSNPs = 0
        self.readFile(bpmFile)
        if normalize: self.normalizeStrand()

    def readFile(self, bpmFile):
        """ read input from .bpm.csv path and update instance variables """
        for line in open(bpmFile, 'r'):
            line = line.replace("\n", "")
            line = line.replace("\r", "")

            fields = line.split(",")

            if line.find("Chromosome") != -1:
                continue

            else:
                self.names.append(fields[1])
                self.chr.append(fields[2])
                self.pos.append(fields[3])
                alleles = fields[5].replace("[", "")
                alleles = alleles.replace("]", "")
                alleles = alleles.split("/")
                self.A.append(alleles[0]) # allele A
                self.B.append(alleles[1]) # allele B
                self.ilmnStrand.append(fields[6][0]) # one of T,B,M,P
                self.normID.append(int(fields[8])) # normalization ID for that snp
        self.numSNPs = len(self.names)

    def normalizeStrand(self):
        """ Implement Illumina normalization to top strand 

        Equivalent normalization done by simtools for other callers at WTSI
        See zCall/src/doc/illumina_strand_normalization.pdf

        Strand is usually TOP (T) or BOT (B), but may be MINUS (M) or PLUS (P)

        Commentary from simtools/Manifest.cpp:

  // input_snp will be something like "[A/C]"; we want to store it as "AC"
  // Illumina method aims to designate A as Allele A on TOP, 
  // and the T as Alelle A on BOT.
  // 
  // Read as [A/B]    Strand      Store as Allele A, Allele B for TOP
  
  // [C/A]             BOT (B)            GT    
  // [G/A]             BOT                CT
  // [G/C]             BOT                CG
  // [T/G]             BOT                AC 
  // [T/C]             BOT                AG
  // [T/A]             BOT                AT
  
  // So BOT to TOP is C->G, A->T, G->C, T->A

  // [A/C]             TOP (T)            AC
  // [A/G]             TOP                AG
  // [A/T]             TOP                AT
  // [C/G]             TOP                CG
  // [C/T]             TOP                CT
  // [G/T]	           TOP                GT
  
  // weird cases
  // [D/I]             M                  DI
  // [D/I]             P                  DI
  // [I/D]             M                  ID
  // [I/D]             P                  ID
  // [N/A]             P                  NA
        """
        
        # apply normalization to alleles in self.A and self.B
        # Use '?' for unknown or missing data
        for i in range(self.numSNPs):
            strand = self.ilmnStrand[i]
            if strand=='T' or strand=='M' or strand=='P':
                continue
            elif strand=='B':
                # normalize A allele
                if self.A[i]=='C': self.A[i] = 'G'
                elif self.A[i]=='G': self.A[i] = 'C'
                elif self.A[i]=='T': self.A[i] = 'A'
                elif self.A[i]=='A': self.A[i] = 'T'
                else: self.A[i] = '?'
                # normalize B allele
                if self.B[i]=='C': self.B[i] = 'G'
                elif self.B[i]=='G': self.B[i] = 'C'
                elif self.B[i]=='T': self.B[i] = 'A'
                elif self.B[i]=='A': self.B[i] = 'T'
                else: self.B[i] = '?'
            else:
                self.A[i] = '?' 
                self.B[i] = '?'


    def getChromosomes(self):
        return self.chr

    def getPositions(self):
        return self.pos

    def getTotalSNPs(self):
        return self.numSNPs
