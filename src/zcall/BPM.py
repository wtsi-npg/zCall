#! /usr/bin/python

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
    def __init__(self, bpmFile):
        self.names = []
        self.chr = []
        self.pos = []
        self.normID = []
        self.A = []
        self.B = []
        self.ilmnStrand = [] # illumina strand designation
        self.numSNPs = 0
        self.readFile(bpmFile)

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


    def getChromosomes(self):
        return self.chr

    def getPositions(self):
        return self.pos

    def getTotalSNPs(self):
        return self.numSNPs
