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

"""Standalone script to convert plain text to input for zCall

Input: Text file
       One sample path per line (no spaces in pathnames!)
       Optional second column is sample name; default to file name
       Optional uri prefix

Output: JSON file with sample names/uri's in correct format
"""

import os, re, sys
import shared
try: 
    import argparse, json
except ImportError: 
    sys.stderr.write(shared.importErrorMessage) 
    sys.exit(1)

class sampleTextParser:

    def __init__(self):
        pass

    def readText(self, inPath):
        """Read GTC paths and (optional) names from plain text file"""
        lines = open(inPath).readlines()
        gtcPaths = {}
        for line in lines:
            if re.match('#', line): continue
            words = re.split('\s+', line.strip())
            gtc = words[0]
            if len(words) > 1: 
                name = words[1]
            else: 
                items = re.split('\.', os.path.basename(gtc))
                items.pop() # remove .gtc suffix
                name = '.'.join(items)
            gtcPaths[name] = gtc
        return gtcPaths

    def writeJson(self, inPath, outPath, uriPrefix):
        """Read plain text and output sample .json"""
        gtcPaths = self.readText(inPath)
        output = []
        names = gtcPaths.keys()
        names.sort()
        uriset = set()
        for name in names:
            sample = {}
            if uriPrefix==None: uri = name
            else: uri = uriPrefix+name
            if uri in uriset:
                raise ValueError("Sample URI \""+uri+"\" is not unique")
            else:
                uriset.add(uri)
            sample['uri'] = uri
            sample['result'] = gtcPaths[name]
            sample['gender'] = -9
            output.append(sample)
        out = open(outPath, 'w')
        out.write(json.dumps(output)+"\n")
        out.close()


def main():
    """Method to run as script from command line.  Run with --help for usage."""
    description = "Process a plain text file listing GTC paths and (optional) sample names. Outputs a JSON file for input to zCall scripts. JSON file is in a similar format to the samples file in the WTSI Genotyping Pipeline."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input', required=True, metavar="PATH", 
                        help="Path to .txt file with sample GTC paths, one per line. Optional second column contains sample names. Paths must not contain whitespace. Input lines starting with # are ignored.")
    parser.add_argument('--output', required=False, metavar="PATH", 
                        help="Path for .json output", default="samples.json")
    parser.add_argument('--uri', required=False, metavar="STRING",
                       help="URI prefix to be prepended to all sample names, eg. uri:my_institution:'", default=None)
    args = vars(parser.parse_args())
    sampleTextParser().writeJson(args['input'], args['output'], args['uri']) 

if __name__ == "__main__":
    main()
