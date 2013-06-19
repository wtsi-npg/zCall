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


"""Combine concordance and gain metric results for multiple samples and zscores.

Use to find 'best' z score for calling.

Author:  Iain Bancarz, ib5@sanger.ac.uk, January 2013
"""

import cProfile, os, sys
try: 
    import argparse, json
    from tempfile import NamedTemporaryFile
    from utilities import ConfigReader
    from calibration import MetricEvaluator
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)


def main():
    """Method to run as script from command line.  Run with --help for usage."""
    description = "Evaluate concordance/gain results for multiple thresholds and samples; find the 'best' z score for subsequent use of zCall."
    configDefault = os.path.join(sys.path[0], '../etc/config.ini')
    configDefault = os.path.abspath(configDefault)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--metrics', required=True, metavar="PATH", 
                        help="Path to .json file containing paths to .json metrics files")
    parser.add_argument('--thresholds', required=True, metavar="PATH", 
                        help="Path to .json file containing threshold .txt paths indexed by z score")
    parser.add_argument('--out', required=True, metavar="PATH", 
                        help="Path for .json output")
    parser.add_argument('--text', required=False, metavar="PATH", default=None,
                       help="Path for text summary of metrics. Optional.")
    parser.add_argument('--config', required=False, metavar="PATH", 
                        help="Path to .ini config file. Optional, defaults to "+configDefault, default=configDefault)
    parser.add_argument('--profile', action='store_true', default=False,
                        help="Write cProfile runtime profile stats to output directory, if not activated by default in config.ini. Profile output can be read with 'python -m pstats PATH'")
    args = vars(parser.parse_args())
    metricPaths = json.loads(open(args['metrics']).read())
    if not os.access(args['config'], os.R_OK):
        raise ValueError("Cannot read .ini path "+args['config'])
    cp = ConfigReader(args['config']).getParser()
    if args['profile'] or cp.has_option('zcall', 'profile'):
        pstats = NamedTemporaryFile(prefix="mergeEvaluation_", 
                                    suffix=".pstats", 
                                    dir=os.path.dirname(args['out']), 
                                    delete=False).name
        a = (args['config'], str(metricPaths), args['thresholds'], 
             args['out'], args['text'])
        if args['text']==None:
            cmd = "MetricEvaluator('%s').writeBest(%s, '%s', '%s', %s)" % a
        else:
            cmd = "MetricEvaluator('%s').writeBest(%s, '%s', '%s', '%s')" % a
        cProfile.run(cmd, pstats)
    else:
        MetricEvaluator(args['config']).writeBest(metricPaths, 
                                                  args['thresholds'], 
                                                  args['out'], args['text'])

if __name__ == "__main__":
    main()
