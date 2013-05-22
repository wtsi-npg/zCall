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

"""Convenience script to generate HTML documentation

Use instead of command-line pydoc. Allows user to specify Python version 2.7, even if installation of pydoc is compiled with an earlier Python version.

Run with -h or --help for command-line help.
"""

import os, pydoc, re, sys
from importlib import import_module
from modulefinder import ModuleFinder
try: 
    import argparse    
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)

def main():

    description="Convenience script to generate HTML documentation using pydoc"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--out', required=True,  metavar="PATH", 
                        help="Directory in which to write HTML output files.")
    parser.add_argument('--recursive', action='store_true', default=False,
                        help="Recursively import documentation for dependencies. If not recursive, zcall documents will contain broken links to standard modules. Recursive mode generates approximately 180 HTML files comprising 6 MB of data.")
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Write pydoc information to stdout.")
    args = vars(parser.parse_args())
    recursive = args['recursive']
    verbose = args['verbose']
    if not verbose: # suppress stdout chatter from pydoc.writedoc
        sys.stdout = open('/dev/null', 'w')

    localDir = os.path.dirname(os.path.realpath(__file__))
    sys.path.append(os.path.abspath(localDir+"/..")) # import from zcall dir
    zcallDir = os.path.abspath(localDir+"/../zcall")
    outDir = os.path.abspath(args['out'])
    if not (os.access(outDir, os.W_OK) and os.path.isdir(outDir)):
        msg = "ERROR: Output path "+outDir+" is not a writable directory.\n"
        sys.stderr.write(msg)
        sys.exit(1)
    os.chdir(outDir)
    import zcall
    pydoc.writedoc(zcall)
    modules = set()
    zcall = set()
    scripts = []
    mf = ModuleFinder()
    for script in os.listdir(zcallDir):
        if re.search("\.py$", script) and script!="__init__.py":
            words = re.split("\.", script)
            words.pop()
            scriptName = (".".join(words)) # name without .py suffix
            modules.add("zcall."+scriptName)
            zcall.add(scriptName)
            scripts.append(script)
    if recursive:
        for script in scripts:
            mf.run_script(os.path.join(zcallDir, script))
            for name, mod in mf.modules.iteritems(): 
                if name not in zcall: modules.add(name)
    for module in modules:
        pydoc.writedoc(import_module(module))

# NB findThresholds and findMeanSD files can only be run as scripts, not imported.  Omitting the .py extension from these files prevents pydoc from creating broken links in the main zcall page.


if __name__ == "__main__":
    main()
