#! /usr/bin/env python

"""Process Plink format genotyping data

See http://pngu.mgh.harvard.edu/~purcell/plink/
"""

import struct
from BPM import *

class PlinkHandler:
    """Class to handle Plink format data"""

    def __init__(self, bpm):
        """Initialize with a BPM object"""
        self.bpm = bpm
        self.snpTotal = self.bpm.getTotalSNPs()
        self.sortMap = self.snpSortMap()

    def callsToBinary(self, calls, reorder=True):
        """Translate genotype calls for one sample to Plink binary

        4 genotype calls are packed into one byte of output
        Returns a list of struct.pack strings corresponding to output bytes"""
        if len(calls) != self.snpTotal:
            raise ValueError("Number of calls is not equal to SNP total!")
        if reorder:
            sortedCalls = [None]*self.snpTotal
            for i in range(self.snpTotal):
                sortedCalls[self.sortMap[i]] = calls[i] 
            calls = sortedCalls
        if self.snpTotal % 4 != 0:
            # if not an integer number of bytes, pad with no calls
            calls.extend([0]*(self.snpTotal % 4)) 
        output = []
        i = 0
        while i < self.snpTotal:
            byte = struct.pack('B', self.callsToByte(calls[i:i+4]))
            output.append(byte)
            i += 4
        return output

    def callsToByte(self, calls):
        """Convert list of 4 calls to an integer in Plink binary format 

        Create byte string of the form '01001101', convert to integer
        See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
        """
        if len(calls) != 4:
            raise ValueError("Must have exactly 4 calls for byte conversion!")
        byte = []
        for call in calls:
            if call==1: bcall = '00' # major homozygote, 'AA'
            elif call==2: bcall = '01' # heterozygote, 'AB'
            elif call==3: bcall = '11' # minor homozygote, 'BB'
            else: bcall = '10' # missing genotype, call=0
            byte.append(bcall)
        byteString = ''.join(byte)
        byteString = byteString[::-1] # reverse order of string characters
        return int(byteString, 2)

    def numericChromosomes(self, chroms):
        """Convert to numeric chromosome IDs used by Plink"""
        for i in range(len(chroms)):
            if chroms[i]=='X': chroms[i] = 23
            elif chroms[i]=='Y': chroms[i] = 24
            elif chroms[i]=='XY': chroms[i] = 25
            elif chroms[i]=='MT': chroms[i] = 26
            else: chroms[i] = int(chroms[i])
        return chroms

    def snpSortMap(self):
        """Sort snps into (chromosome, position) order

        Ensures compatibility with sorted .bim files generated by Plink
        Return a map from original position to sorted position"""
        chroms = self.numericChromosomes(self.bpm.getChromosomes())
        pos = self.bpm.getPositions()
        coords = [None]*self.snpTotal
        for i in range(self.snpTotal):
            coords[i] = (chroms[i], int(pos[i]), i)
        coords.sort()
        sortMap = {}
        for i in range(self.snpTotal):
            [chrom, pos, orig] = coords[i]
            sortMap[orig] = i
        return sortMap

    def writeBed(self, binaryCalls, outPath, verbose=False):
        """Write output for one or more samples in Plink .bed format

        Input: List of call bytes, and output path
        Output file:  First 2 bytes are Plink magic number
        3rd byte is flag for an individual-major file
        Subsequent bytes represent genotype calls
        """
        header = [0b01101100, 0b00011011, 0b00000000]
        output = []
        for byte in header: output.append(struct.pack('B', byte))
        output.extend(binaryCalls)
        out = open(outPath, 'w')
        for byte in output: out.write(byte)
        out.close()
        if verbose: print len(output), "bytes written."

    def writeBim(self):
        """Write a Plink .bim file to accompany .bed output

        Similar to Plink .map format, except:
        - 2 additional columns for allele names (use A and B as dummy values)
        - Entries are sorted into (chromosome, position) order
        """
        pass
