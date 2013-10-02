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


import struct
import math
from cStringIO import StringIO

class GTC:
    ''' Class for parsing GTC file'''

    def __init__(self, inputPath, bpmNormIDs):
        ''' Init function for class. Input is a .gtc file '''
        self.inputPath = inputPath
        self.f = StringIO(open(self.inputPath, 'rb').read())
        #self.f = open(self.inputPath, 'rb') # open file handler for binary file
        self.BPMnormIDs = bpmNormIDs # list with norm ID for each snp
        self.TOC = self.parseTOC() # parse table of contents to get location of other information
        self.numSNPs = self.readNumSNPs()
        self.sampleName = self.parseString(10) # parse sample name
        self.samplePlate = self.parseString(11) # parse sample plate
        self.sampleWell = self.parseString(12) # parse sample well
        self.clusterFile = self.parseString(100) # parse what cluster file was used
        self.snpManifest = self.parseString(101) # parse what snp manifest was used
        self.imagingDate = self.parseString(200) # parse imaging date
        self.autoCallDate = self.parseString(201) # parse autocall date
        self.autoCallVersion = self.parseString(300) # parse autocall version
        self.rawXintensities = self.extractIntensities(1000) # parse raw x intensities into python list object
        self.rawYintensities = self.extractIntensities(1001) # parse raw y intensities into python list object
        self.normalizationTransformations = self.extractNormalizationTransformations(400) # parse normalization transformation arrays into python dictionary where key is the order they appeared in the gtc file and the value is a dictionary with keys offset_x, offset_y,scale_x, scale_y, shear, theta and values are floats        
        self.genotypes = self.extractGenotypes(1002) # parse genotypes (0,1,2,3) into python list object
        self.baseCalls = self.extractBaseCalls(1003) # parse basecalls (AT,TT,AT,--) into python list object

        self.normXintensities, self.normYintensities = self.normalizeIntensities()


    def getInputPath(self):
        """Return original input path"""
        return self.inputPath
          
    def parseTOC(self):
        '''Parse Table of Contents of GTC file
        No input
        Output is a dictionary where the ID for that entry is the key and the value is the offset for that variable in the GTC file
        '''
        self.f.seek(4,0)
        line = self.f.read(4)

        count = struct.unpack("i",line)[0] 
        TOC = {}

        i = 0
        while i < count:
            line = self.f.read(2)
            id = struct.unpack("h",line)
            line = self.f.read(4)            
            offset = struct.unpack("I",line)
            TOC[id[0]] = offset[0]
            i+=1
        return TOC


    def parseString(self,id):
        '''
        Extract a string variable from GTC file such as SampleName.
        Input is ID for that variable in TOC
        Output is a string
        '''        
        offset = self.TOC[id]
        self.f.seek(offset,0)

        line = self.f.read(1)
        nbytes = struct.unpack("b",line)[0]
        line = self.f.read(nbytes)
        type = nbytes * "s"
        x = "".join(list(struct.unpack(type, line)))
        
        return x

    def extractIntensities(self, id):
        '''
        Extract intensity values (x or y depending on input ID).
        Input is ID for variable of interest in TOC
        Output is a list with integer intensity values in the order they were parsed
        '''        
        intensities = [0]*self.numSNPs
        offset = self.TOC[id]

        self.f.seek(offset,0)
        line = self.f.read(4)
        count = struct.unpack("i",line)[0]
        if count != self.numSNPs:
            msg = "Inconsistent SNP count in GTC file intensity group: "+\
                "Expected %s, found %s" % (self.numSNPs, count)
            raise ValueError(msg)
        #self.numSNPs = count
        i = 0
        while i < count:
            line = self.f.read(2)
            y = struct.unpack("H",line)
            intensities[i] = y[0]
            i += 1
        return intensities

    def extractNormalizationTransformations(self, id):
        '''
        Extract normalization transformation arrays
        Input is ID for Normalization Transformations in TOC.
        Output is dictionary where keys are the order xForm array appears in gtc file (ex: 1,2,3...).
        The values of the dictionary are another dictionary
        where the keys are shear, offset_x, offset_y, theta, scale_x, scale_y and the values are floats
        '''
        normTransforms = {}
        offset = self.TOC[id]

        self.normIDlist = list(set(self.BPMnormIDs)) # ordered list of unique normIDs
        self.normIDlist.sort()
        
        self.f.seek(offset,0)
        line = self.f.read(4)
        count = struct.unpack("i",line)[0]

        i = 0
        while i < count:
            line = self.f.read(4)
            line = self.f.read(48)
            x = struct.unpack("<12f", line)
            normTransforms[self.normIDlist[i]] = {"offset_x":x[0],"offset_y":x[1],"scale_x":x[2],"scale_y":x[3],"shear":x[4],"theta":x[5]}
            i += 1
        return normTransforms
    
    def extractBaseCalls(self, id):
        '''
        Extract base calls.
        Input is id for BaseCalls in TOC
        Output is a list with one basecall for each SNP (ex: AT, GT,AA...)
        '''
        baseCalls = [None]*self.numSNPs
        offset = self.TOC[id]
        self.f.seek(offset,0)
        line = self.f.read(4)
        count = struct.unpack("i",line)[0]
        if count != self.numSNPs:
            msg = "Inconsistent SNP count in GTC file intensity group: "+\
                "Expected %s, found %s" % (self.numSNPs, count)
        i = 0
        while i < count:
            line = self.f.read(2)
            calls = struct.unpack("ss",line)
            baseCalls[i] = calls[0] + calls[1]
            i += 1
        return baseCalls

    def extractGenotypes(self, id):
        '''
        Extract genotypes.
        Input is ID for Genotypes in TOC
        Output is a list with one genotype per SNP (0,1,2,3)
        '''
        genotypes = [None]*self.numSNPs
        offset = self.TOC[id]
        self.f.seek(offset,0)
        line = self.f.read(4)
        count = struct.unpack("i",line)[0]
        if count != self.numSNPs:
            msg = "Inconsistent SNP count in GTC file intensity group: "+\
                "Expected %s, found %s" % (self.numSNPs, count)
        i = 0
        while i < count:
            line = self.f.read(1)
            gt = struct.unpack("b",line)
            genotypes[i] = gt[0]
            i += 1
        return genotypes

    def getTotalSNPs(self):
        return self.numSNPs

    def normalizeIntensities(self):
        '''
        Use Normalization transformations to convert raw intensities to normalized intensities
        No Input
        Outputs are normalized x and y intensities in python lists
        '''
        normXIntensities = [0]*self.numSNPs
        normYIntensities = [0]*self.numSNPs
        
        i = 0
        while i < self.numSNPs:
            xraw = self.rawXintensities[i]
            yraw = self.rawYintensities[i]
            normID = self.BPMnormIDs[i]

            offset_x = self.normalizationTransformations[normID]["offset_x"]
            offset_y = self.normalizationTransformations[normID]["offset_y"]
            scale_x = self.normalizationTransformations[normID]["scale_x"]
            scale_y = self.normalizationTransformations[normID]["scale_y"]
            theta = self.normalizationTransformations[normID]["theta"]
            shear = self.normalizationTransformations[normID]["shear"]

            tempx = xraw - offset_x
            tempy = yraw - offset_y

            tempx2 = math.cos(theta) * tempx + math.sin(theta) * tempy
            tempy2 = -1 * math.sin(theta) * tempx + math.cos(theta) * tempy

            tempx3 = tempx2 - (shear * tempy2)
            tempy3 = tempy2

            try: xn = tempx3 / float(scale_x)
            except ZeroDivisionError: xn = 0.0
            try: yn = tempy3 / float(scale_y)
            except ZeroDivisionError: yn = 0.0

            if xn < 0:
                xn = 0.0
            if yn < 0:
                yn = 0.0

            normXIntensities[i] = xn
            normYIntensities[i] = yn
            i += 1

        return (normXIntensities, normYIntensities)            

    def readNumSNPs(self):
        """Read total number of SNPs from first block of intensity data

        Adapted from extractIntensities"""
        offset = self.TOC[1000]
        self.f.seek(offset,0)
        line = self.f.read(4)
        count = struct.unpack("i",line)[0]
        return count
