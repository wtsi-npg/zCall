#! /usr/bin/env python

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

import os, json, sys
from math import isinf, isnan
from ConfigParser import ConfigParser
from GTC import *
from BPM import *
from EGT import *

class SharedBase:
    """Constants and simple functions shared across multiple classes"""

    Z_KEY = 'BEST_Z'
    T_KEY = 'BEST_THRESHOLDS'
    M_KEY = 'SAMPLE_METRICS'

    def readSampleJson(self, inPath):
        """Read sample GTC paths from .json file used by genotyping pipeline"""
        contents = open(inPath).read()
        try: 
            samples = json.loads(contents)
        except ValueError: 
            msg = "JSON parser error! Input to failed parse attempt:\n"
            sys.stderr.write(msg)
            sys.stderr.write(contents)
            raise
        gtc = []
        for sample in samples: gtc.append(sample['result'])
        return gtc
        

class CallingBase(SharedBase):
    """ 'Base' class containing useful methods for zcall subclasses """

    def __init__(self, bpmPath, egtPath, threshPath=None):
        self.bpm = BPM(bpmPath)
        self.egt = EGT(egtPath)
        self.snpTotal = self.egt.getTotalSNPs()
        if self.bpm.getTotalSNPs() != self.snpTotal:
            raise ValueError("ERROR: SNP totals in .egt and .bpm inputs differ")
        if threshPath!=None:
            self.thresholds = ThresholdContainer(threshPath)
        else:
            self.thresholds = None

    def call(self, gtc, i):
        """ re-call ith SNP in GTC file, using zcall thresholds
        
        call codes: 0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        """
        if self.thresholds==None:
            raise ValueError("Must specify thresholds before calling!")
        normX = gtc.normXintensities[i]
        normY = gtc.normYintensities[i]
        Tx = self.thresholds.getX(i)
        Ty = self.thresholds.getY(i)
        call = None
        inputSum = normX+normY+Tx+Ty # use sum to check for NaN or infinity
        if (isnan(inputSum) or isinf(inputSum)):
            ## NaN or infinity present in inputs -- no call
            call = 0
        elif normX < Tx and normY < Ty: ## Lower left quadrant
            call = 0
        elif normX >= Tx and normY <= Ty: ## Lower right quadrant
            call = 1
        elif normX < Tx and normY >= Ty: ## Upper left quadrant
            call = 3
        else: ## Upper right quadrant
            call = 2
        return call

    def findMAF(self, nAA, nBB, nAB):
        """ Find minor allele frequency """
        maf = None
        try:
            if nAA > nBB:
                maf = (nAB + 2 * nBB) / float(2*(nAA + nAB + nBB))
            else:
                maf = (nAB + 2 * nAA) / float(2*(nAA + nAB + nBB))
        except ZeroDivisionError:
            maf = 0
        return maf

    def normalizeCall(self, call, nAA, nBB):
        """Normalization:  Flip genotype call so 1 is always the common allele homozygote and 3 is the minor allele homozygote 

        Enforces convention that major allele is on X intensity axis
        Allele counts are taken from EGT object
        """
        if nBB > nAA: 
            if call == 1:
                call = 3
            elif call == 3:
                call = 1
        return call

    def setThresholds(self, thresholds):
        """Set thresholds to given ThresholdContainer object"""
        self.thresholds = thresholds

class ConfigReader:
    """Convenience wrapper for the ConfigParser class to read .ini files"""

    def __init__(self, configPath):
        configFile = open(configPath, 'r')
        self.cp = ConfigParser()
        self.cp.readfp(configFile)
        configFile.close()

    def getParser(self):
        return self.cp

class ArgParserExtra:
    """Class for additional parsing/validation of command-line arguments"""
    def __init__(self, args):
        self.args = args
    
    def enableProfile(self):
        """Check if script should be run with Python profiler enabled

        Note that parser converts - to _ in --no-profile option key"""
        cp = ConfigReader(self.args['config']).getParser()
        if cp.has_option('zcall', 'profile'): default = True
        else: default = False
        if self.args['profile'] and self.args['no_profile']:
            raise ValueError("Cannot specify both --profile and --no-profile!")
        elif self.args['profile']:
            enable = True
        elif self.args['no_profile']:
            enable = False
        else:
            enable = default
        return enable

    def validateInputs(self, keys):
        """Validate input paths with given keys

        Convert paths to absolute paths, and return revised args"""
        for key in keys:
            if not os.path.exists(self.args[key]):
                raise OSError("Input path '"+self.args[key]+"' does not exist!")
            elif not os.access(self.args[key], os.R_OK):
                raise OSError("Cannot read input path '"+self.args[key]+"'!")
            else:
                self.args[key] = os.path.abspath(self.args[key])
        return self.args

    def validateOutputDirPath(self, dirPath):
        """Validate output directory path, and convert to absolute path"""
        if not os.path.exists(dirPath):
            raise OSError("Output '"+dirPath+"' does not exist!")
        elif not os.path.isdir(dirPath):
            raise OSError("Output '"+dirPath+"' is not a directory!")
        elif not os.access(dirPath, os.W_OK):
            raise OSError("Cannot write to output '"+dirPath+"'!")
        else:
            dirPath = os.path.abspath(dirPath)
        return dirPath

    def validateOutputDir(self, key="out"):
        """Validate output directory argument

        Convert path to absolute path, and return revised args"""
        self.args[key] = self.validateOutputDirPath(self.args[key])
        return self.args

    def validateOutputFile(self, key="out"):
        """Validate an output file path

        Directory must exist and be writable, file may not be"""
        if os.path.isdir(self.args[key]):
            raise OSError("Output file '"+self.args[key]+"' is a directory!")
        (dirName, fileName) = os.path.split(self.args[key])
        dirName = self.validateOutputDirPath(dirName)
        self.args[key] = os.path.join(dirName, fileName)
        return self.args
        

class ThresholdContainer:
    """Class to contain thresholds, read from thresholds.txt file"""

    def __init__(self, inPath):
        (self.x, self.y) = self.readThresholds(inPath)
        
    def getX(self, i):
        """Get x threshold for the ith SNP"""
        return self.x[i]

    def getY(self, i):
        """Get y threshold for the ith SNP"""
        return self.y[i]

    def readThresholds(self, inPath):
        """ Read a thresholds.txt file; return lists of x and y thresholds """
        thresholdsX = []
        thresholdsY = []
        for line in open(inPath, 'r'):
            line = line.replace("\n", "")
            if line.find("Tx") != -1:
                continue
            else:
                fields = line.split("\t")
                if fields[1] != "NA":
                    tx = float(fields[1])
                else:
                    tx = fields[1]
                if fields[2] != "NA":
                    ty = float(fields[2])
                else:
                    ty = fields[2]
                thresholdsX.append(tx)
                thresholdsY.append(ty)
        return (thresholdsX, thresholdsY)
