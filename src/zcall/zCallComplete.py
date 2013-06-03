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


"""Standalone, self-contained script to run the complete zcall process.

zcall steps:
1. Generate thresholds
2. Evaluate thresholds
3. Merge evaluations to find best threshold
4. Apply zcall with given threshold to no-calls in input data

This standalone script does not allow parallelization of evaluation or calling, 
or reuse of thresholds. It is intended as a convenience script for small
datasets where parallelization is not required, and an illustration of the 
zcall method.
"""
import cProfile, os, sys, time
try: 
    import argparse, json
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)
from calibration import ThresholdFinder, SampleEvaluator, MetricEvaluator
from runZCall import SampleCaller

class ZCallComplete:

    EVALUATION = 'threshold_evaluation.json'
    MERGED = 'merged_evaluation.json'
    LOG = 'zcall_log.json'

    def __init__(self, args):
        self.args = args

    def call(self, bpm, egt, thresholdPath, sJson, outDir, prefix, binary,
             verbose):
        if verbose: print "Running zcall with thresholds", thresholdPath
        logPath = os.path.join(outDir, self.LOG)
        caller = SampleCaller(bpm, egt, thresholdPath)
        caller.run(sJson, outDir, prefix, logPath, binary, verbose)

    def evaluate(self, bpm, egt, thresholds, gtc, start, end, outDir, verbose):
        """Evaluate threshold.txt files by concordance and gain metrics"""
        outPath = os.path.join(outDir, self.EVALUATION)
        eva = SampleEvaluator(bpm, egt)
        eva.run(thresholds, gtc, start, end, outPath, verbose)
        return outPath

    def merge(self, metricPath, thresholdJson, outDir, textPath, config, 
              verbose):
        outPath = os.path.join(outDir, self.MERGED)
        eva = MetricEvaluator(config)
        results = eva.writeBest([metricPath,], thresholdJson, outPath, 
                                textPath, verbose)
        thresholdPath = results[eva.getBestThresholdKey()]
        return thresholdPath

    def prepare(self, zstart, ztotal, egtPath, outDir, config, verbose=False):
        """ Prepare threshold.txt files for given range of z scores"""
        tf = ThresholdFinder(egtPath, config)
        return tf.runMultiple(zstart, ztotal, outDir, verbose)
        

    def run(self):
        tJson = self.prepare(self.args['zstart'],
                             self.args['ztotal'],
                             self.args['egt'],
                             self.args['out'],
                             self.args['config'],
                             self.args['verbose'])
        mJson = self.evaluate(self.args['bpm'],
                              self.args['egt'],
                              tJson,
                              self.args['samples'],
                              self.args['gtc_start'],
                              self.args['gtc_end'],
                              self.args['out'],
                              self.args['verbose'])
        thresholdPath = self.merge(mJson, 
                                   tJson, 
                                   self.args['out'], 
                                   self.args['text'], 
                                   self.args['config'],
                                   self.args['verbose'])
        self.call(self.args['bpm'],
                  self.args['egt'],
                  thresholdPath,
                  self.args['samples'],
                  self.args['out'],
                  self.args['plink'],
                  self.args['binary'],
                  self.args['verbose'])
        

def main():
    args = parseArgs()
    start = time.time()
    if args['profile']==True:
        cProfile.run('ZCallComplete('+str(args)+').run()')
    else:
        ZCallComplete(args).run()
    if args['verbose']==True:
        duration = time.time() - start
        print "zCall finished. Duration:", round(duration, 2), "seconds."

def parseArgs():
    description = "Standalone script to run the complete zcall process: Generate and evaluate thresholds, and apply zcall to no-calls in the input data."
    parser = argparse.ArgumentParser(description=description)
    configDefault = os.path.join(sys.path[0], '../etc/config.ini')
    configDefault = os.path.abspath(configDefault)
    parser.add_argument('--bpm', required=True, metavar="PATH", 
                        help="BPM .csv manifest file")
    parser.add_argument('--egt', required=True, metavar="PATH", 
                        help="Path to .egt input file.")
    parser.add_argument('--config', metavar="PATH", default=configDefault,
                        help="Path to .ini config file. Default = etc/config.ini")
    parser.add_argument('--out', metavar="DIR", default=".",
                        help="Directory for output; defaults to current working directory.")
    parser.add_argument('--plink', default='zcall', metavar="STRING", 
                        help="Prefix for Plink output files")
    parser.add_argument('--binary', action='store_true', default=False,
                        help="Write Plink binary output. If this option is not given, output is in Plink text format.")
    parser.add_argument('--zstart', metavar="INT", default=7, type=int,
                    help='Starting z score. Default = %(default)s')
    parser.add_argument('--ztotal', metavar="INT", default=1, type=int,
                        help='Total number of integer z scores to generate. Default = %(default)s')
    parser.add_argument('--gtc_start', metavar="INT", default = 0,
                        help="Starting index in GTC .json file for threshold evaluation")
    parser.add_argument('--gtc_end', metavar="INT", default = -1,
                        help="Ending index in GTC .json file for threshold evaluation")
    parser.add_argument('--samples', required=True, metavar="PATH", 
                        help="Path to .json file containing sample URIs (unique identifiers), gender codes (optional), and .gtc data paths")
    parser.add_argument('--text', required=False, metavar="PATH", 
                       help="Path for text summary of calibration metrics. Optional.", 
                       default=None)
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print status information to standard output")
    parser.add_argument('--profile', action='store_true', default=False,
                        help="Use cProfile to profile runtime operation")
    args = vars(parser.parse_args())
    inputKeys = ('bpm', 'egt', 'config')
    for key in inputKeys:
        if not os.access(args[key], os.R_OK):
            msg = "Cannot read path: \""+args[key]+"\"\n"
            sys.stderr.write(msg)
            sys.exit(1)
        else:
            args[key] = os.path.abspath(args[key])
    if not os.path.isdir(args['out']) or  not os.access(args['out'], os.W_OK):
        msg = "Output path \""+args['out']+"\" is not a writable directory!\n"
        sys.stderr.write(msg)
        sys.exit(1)
    

    return args

if __name__ == "__main__":
    main()
