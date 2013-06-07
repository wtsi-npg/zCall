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


"""Evaluate concordance and gain for given GTC files and Z scores.

Author:  Iain Bancarz, ib5@sanger.ac.uk, January 2013
"""

import cProfile, os, sys
from calibration import SampleEvaluator
try: 
    import argparse     # optparse is deprecated, using argparse instead
    from tempfile import NamedTemporaryFile
    from utilities import ConfigReader
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)

def main():
    """Method to run as script from command line.  Run with --help for usage."""
    description = "Evaluate concordance/gain for given z scores and GTC files."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--thresholds', required=True, metavar="PATH", 
                        help="Path to .json file containing threshold .txt paths indexed by z score")
    parser.add_argument('--bpm', required=True, metavar="PATH", 
                        help="BPM .csv manifest file")
    parser.add_argument('--egt', required=True, metavar="PATH", 
                        help="EGT input file")
    parser.add_argument('--gtc', required=True, metavar="PATH", 
                        help="Path to .json file containing .gtc input paths")
    parser.add_argument('--out', required=True, metavar="PATH", 
                        help="Path for .json output")
    parser.add_argument('--start', metavar="INT", default = 0, type=int,
                        help="Starting index in GTC .json file")
    parser.add_argument('--end', metavar="INT", default = -1, type=int,
                        help="Ending index in GTC .json file")
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print status information to standard output")
    parser.add_argument('--profile', action='store_true', default=False,
                        help="Use cProfile to profile runtime operation, if not activated by default in config.ini")
    configDefault = os.path.join(sys.path[0], '../etc/config.ini')
    configDefault = os.path.abspath(configDefault)
    parser.add_argument('--config', metavar="PATH", default=configDefault,
                        help="Path to .ini config file. Default = etc/config.ini")
    args = vars(parser.parse_args())
    inputKeys = ['thresholds', 'bpm', 'egt', 'gtc']
    for key in inputKeys:
        if not os.access(args[key], os.R_OK):
            raise OSError("Cannot read input path \""+args[key]+"\"")
        else:
            args[key] = os.path.abspath(args[key])
    (dirName, fileName) = os.path.split(os.path.abspath(args['out']))
    if fileName=='' or not os.access(dirName, os.R_OK):
        raise OSError("Invalid output path \""+args['out']+"\"")
    cp = ConfigReader(os.path.abspath(args['config'])).getParser()
    if args['profile'] or cp.has_option('zcall', 'profile'):
        pstats = NamedTemporaryFile(prefix="evaluateThresholds_", 
                                    suffix=".pstats", 
                                    dir=os.path.dirname(args['out']), 
                                    delete=False).name
        cmd0 = "SampleEvaluator('%s', '%s')" % (args['bpm'], args['egt'])
        args1 = (args['thresholds'], args['gtc'], args['start'], 
                 args['end'], args['out'], args['verbose'])
        cmd1 = ".run('%s', '%s', %s, %s, '%s', %s)" % args1
        cProfile.run(cmd0+cmd1, pstats)
    else:
        eva = SampleEvaluator(args['bpm'], args['egt'])
        eva.run(args['thresholds'], args['gtc'], args['start'], 
                args['end'], args['out'], args['verbose'])

if __name__ == "__main__":
    main()
