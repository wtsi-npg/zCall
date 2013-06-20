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


""" Find thresholds for given .egt file and Z score(s)
Combines findMeanSD.py, findBetas.r, findThresholds.py from original zCall

Iain Bancarz, ib5@sanger.ac.uk
January 2013
"""

import cProfile, os, sys, time
try: 
    import argparse, json
    from tempfile import NamedTemporaryFile
    from calibration import ThresholdFinder
    from utilities import ArgParserExtra
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)

"""
Calibration procedure:
1. Run findMeanSD.py on given EGT file.  Outputs mean_sd.txt
2. Run findBetas.r on output from (1). Outputs betas.txt
3. Run findThresholds.py on EGT file and betas.txt, with given Z score(s).
Outputs from (1) and (2) are written to a temporary directory, deleted on exit.

Recommended default Z score = 7.  Suggested range of alternatives = 3 to 15.
"""

def main():
    # 'main' method to run script from command line
    start = time.time()
    args = get_args()
    verbose = args['verbose']
    if verbose: print "Starting prepareThresholds.py"
    # validate arguments
    if args['ztotal']<1 or args['zstart']<1:
        raise ValueError("Invalid zstart or ztotal option.")
    parserExtra = ArgParserExtra(args)
    args = parserExtra.validateInputs(['egt', 'config'])
    args = parserExtra.validateOutputDir()
    profile = parserExtra.enableProfile()
    if profile:
        pstats = NamedTemporaryFile(prefix="prepareThresholds_", 
                                    suffix=".pstats", 
                                    dir=out, delete=False).name
        cmd = "ThresholdFinder('"+args['egt']+"', '"+config+\
            "').runMultiple(%s, %s, '%s', %s, %s)" % \
            (args['zstart'], args['ztotal'], args['out'], args['verbose'], 
             args['force'])
        cProfile.run(cmd, pstats)
    else:
        tf = ThresholdFinder(args['egt'], args['config'])
        tf.runMultiple(args['zstart'], args['ztotal'], args['out'], 
                       args['verbose'], args['force'])
        duration = time.time() - start
        if verbose: print "Finished. Duration:", round(duration, 1), "s"
    
def get_args():
    # parse command-line arguments and return dictionary of params
    description = "Generates threshold files for use with the zCall genotype caller.  Inputs are an .egt file and one or more Z score values.  The .egt file is a proprietary Illumina binary file format, containing typical means and standard deviations for intensity clusters.  An .egt file is supplied by Illumina for its own genotyping chips, or it may be generated using the GenomeStudio software for custom probe sets."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--egt', required=True, metavar="PATH", 
                        help="Path to .egt input file.")
    parser.add_argument('--out', metavar="DIR", default=".",
                        help="Directory for output; defaults to current working directory.  Filename(s) will be of the form <prefix>_z<zscore>_thresholds.txt, for an input file of the form <prefix>.egt")
    configDefault = os.path.join(sys.path[0], '../etc/config.ini')
    configDefault = os.path.abspath(configDefault)
    parser.add_argument('--config', metavar="PATH", default=configDefault,
                        help="Path to .ini config file. Default = etc/config.ini")
    parser.add_argument('--zstart', metavar="INT", default=7, type=int,
                    help='Starting z score. Default = %(default)s')
    parser.add_argument('--ztotal', metavar="INT", default=1, type=int,
                        help='Total number of integer z scores to generate. Default = %(default)s')
    parser.add_argument('--index_name', metavar="STRING", 
                        default="threshold_index.json",
                        help='Name for .json index file with paths to thresholds.txt output, written to output directory')
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print status information to standard output")
    parser.add_argument('--profile', action='store_true', default=False,
                        help="Use cProfile to profile runtime operation, if not activated by default in config.ini")
    parser.add_argument('--no-profile', action='store_true', default=False,
                        help="Do not use cProfile to profile runtime operation. Overrides default in config.ini.")
    parser.add_argument('--force', action='store_true', default=False,
                        help="Force overwrite of existing threshold files (if any)")
    args = vars(parser.parse_args())
    return args

if __name__ == "__main__":
    main()
