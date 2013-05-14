#! /software/bin/env python

"""Standalone script to convert plain text to input for zCall

Input: Text file
       One sample path per line (no spaces in pathnames!)
       Optional second column is sample name; default to file name
       Optional uri prefix

Output: JSON file with sample names/uri's in correct format
"""

import os, re, sys
try: 
    import argparse, json
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
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
                       help="URI prefix to be prepended to all sample names, eg. uri:my_institution'", default=None)
    args = vars(parser.parse_args())
    sampleTextParser().writeJson(args['input'], args['output'], args['uri']) 

if __name__ == "__main__":
    main()
