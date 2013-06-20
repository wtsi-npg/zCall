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

"""Unit tests for zcall module

Test cases:
- prepareThresholds.py creates correct thresholds.txt files
- evaluateThresholds.py creates correct evaluation in .json format
- mergeEvaluation.py correctly merges results and identifies best Z score
- call.py produces Plink binary data

Required input data:
- Example BPM, EGT files
- Example GTC files -- may have privacy issues
- (TODO Test on public GTC data, with corresponding BPM/EGT)
- (TODO Create and test a one-step "wrapper" script to calibrate, evaluate and call without any parallelization)

Author:  Iain Bancarz, ib5@sanger.ac.uk
"""

import json, os, sys, unittest
from ConfigParser import ConfigParser
from hashlib import md5
from tempfile import mkdtemp

class TestScripts(unittest.TestCase):

    """Test command-line python scripts used in WTSI genotyping pipeline"""

    def call(self, binary):
        """Re-call GTC files using zCall and choice of output"""
        outStem = os.path.join(self.outDir, self.prefix)
        tPath = os.path.join(self.bigData, 'thresholds_HumanExome-12v1_z07.txt')
        logPath = os.path.join(self.outDir, 'zcall_log.json')
        args = ['zcall/runZCall.py',
                '--thresholds', tPath,
                '--bpm', self.bpmPath,
                '--egt', self.egtPath,
                '--samples', self.sampleJson,
                '--out', self.outDir,
                '--plink', self.prefix,
                '--log', logPath,
            ]
        if binary: args.append('--binary')
        self.assertEqual(os.system(' '.join(args)), 0) # run script
        if binary:
            suffixes = ['.bed', '.bim', '.fam']
            expected = ['55fa3cfd960d43366cb506ab004ac300',
                        '19dd8929cd63e3ee906e43b9bb59cd02',
                        'b836bd45459de6a9bc4b8d92e8c9e298']
        else:
            suffixes = ['.ped', '.map', '.fam']
            expected = ['731e0199b1433228abf80ec8126089d1',
                        'df64da4ef76b2df1040df375a8933f45',
                        'b836bd45459de6a9bc4b8d92e8c9e298']
        for i in range(len(suffixes)):
            self.assertTrue(os.path.exists(outStem+suffixes[i]))
            checksum = self.getMD5hex(outStem+suffixes[i])
            self.assertEqual(checksum, expected[i])  
        self.assertTrue(os.path.exists(logPath))
        self.validatePlink(self.prefix, binary)


    def getMD5hex(self, inPath):
        """Get MD5 checksum for contents of given file, in hex format"""
        m = md5()
        m.update(open(inPath).read())
        checksum = m.hexdigest()
        return checksum

    def readConfig(self, configPath=None):
        """Read local params from config file

        - bigdata: Directory for test files too big to upload to github"""
        if configPath==None:
            configPath = os.path.abspath('etc/config.ini')
        if not os.access(configPath, os.R_OK):
            raise IOError("Cannot read config path '"+configPath+"'")
        config = ConfigParser()
        config.readfp(open(configPath))
        bigData = config.get('test', 'bigdata')
        return bigData

    def validateThresholds(self, jsonPath):
        """Check that thresholds.txt files are correct

        Use for test_prepareThresholds, and to validate input for other tests
        """
        index = json.loads(open(jsonPath).read())
        for z in index.keys():
            self.assertEqual(self.getMD5hex(index[z]), self.expectedT[z])

    def validatePlink(self, prefix, binary):
        """Check that Plink can parse dataset without crashing"""
        startDir = os.getcwd()
        os.chdir(self.outDir)
        if binary: opt = '--bfile'
        else: opt = '--file'
        cmd = 'plink '+opt+' '+prefix+' > /dev/null'
        try:
            self.assertEqual(0, os.system(cmd))
        except AssertionError:
            os.chdir(startDir)
            raise
        os.chdir(startDir)

    def setUp(self):
        """Check for valid input/output files and directories"""
        self.prefix = 'test'
        self.dataDir = 'data'
        self.bigData = self.readConfig()
        for d in (self.dataDir, self.bigData):
            if not os.path.exists(d) or not os.path.isdir(d):
                msg = "Invalid test directory: \""+d+"\"\n"
                sys.stderr.write(msg)
                sys.exit(1)
        self.expectedT = {"6":"1a53e8cbba6750d43d5ff607cf616beb",
                          "7":"a8d8b62be728b62fc986230da13f4ef7",
                          "8":"1f14419d0053841cfa8ab3fb994de1c1"}
        self.outDir = mkdtemp(dir=self.dataDir)
        print "Created output directory", self.outDir
        self.bpmPath = os.path.join(self.bigData, 'HumanExome-12v1_A.bpm.csv')
        self.egtPath = os.path.join(self.bigData, 'HumanExome-12v1.egt')
        # TODO automatically generate threshold.txt files, if not present
        self.gtcPath = os.path.join(self.dataDir, 'gtc.json')
        self.metricIndex = os.path.join(self.bigData,
                                        'evaluation_metric_index.json')
        self.thresholdJsonName = 'thresholds.json'
        self.thresholdJson = os.path.join(self.bigData, self.thresholdJsonName)
        if os.path.exists(self.thresholdJson):
            self.validateThresholds(self.thresholdJson)
        else:
            sys.stderr.write("WARNING: Missing thresholds, see test/README.")
        self.sampleJson = os.path.join(self.dataDir, 'test_sample.json')

    def test_prepareThresholds(self):
        """Prepare thresholds.txt files

        Run as part of normal test suite
        Can also be used to generate input thresholds for other tests
        Checksums should ensure that generated thresholds are correct
        """
        zstart = 6
        ztotal = 3
        outPaths = []
        for i in range(zstart, zstart+ztotal):
            name = 'thresholds_HumanExome-12v1_z0'+str(i)+'.txt'
            outPaths.append(os.path.join(self.outDir, name))
        args = ['zcall/prepareThresholds.py',
                '--egt', self.egtPath,
                '--out', self.outDir,
                '--config etc/config.ini',
                '--zstart', str(zstart),
                '--ztotal', str(ztotal),
                '--index_name', self.thresholdJsonName]
        self.assertEqual(os.system(' '.join(args)), 0) # run script
        jsonOut = os.path.join(self.outDir, self.thresholdJsonName)
        self.assertTrue(os.access(jsonOut, os.R_OK))
        self.validateThresholds(jsonOut)

    def test_evaluateThresholds(self):
        """Evaluate thresholds for collections of sample GTC files"""
        argsBase = ('zcall/evaluateThresholds.py',
                    '--egt', self.egtPath,
                    '--bpm', self.bpmPath,
                    '--thresholds', self.thresholdJson)
        ranges = ((0,4), (4,8))
        for i in range(len(ranges)):
            (start, end) = ranges[i]
            args = list(argsBase)
            name = 'metrics0'+str(i)+'.json'
            outPath = os.path.join(self.outDir, name)
            args.extend(['--gtc', self.gtcPath, '--start', str(start),
                         '--end', str(end), '--out', outPath ])
            self.assertEqual(os.system(' '.join(args)), 0) # run script
            metricsNew = json.loads(open(outPath).read())
            oldPath = os.path.join(self.dataDir, name)
            metricsOld = json.loads(open(oldPath).read())
            self.assertEqual(metricsOld, metricsNew)

    def test_mergeEvaluation(self):
        """Merge evaluation results and find best Z score"""
        outPath = os.path.join(self.outDir, 'merged_evaluation.json')
        args = ['zcall/mergeEvaluation.py',
                '--metrics', self.metricIndex,
                '--thresholds', self.thresholdJson,
                '--out', outPath,
                '--text', os.path.join(self.outDir, 'metric_summary.txt'), ]
        self.assertEqual(os.system(' '.join(args)), 0) # run script
        resultsNew = json.loads(open(outPath).read())
        oldPath = os.path.join(self.dataDir, 'zEvaluation.json')
        resultsOld = json.loads(open(oldPath).read())
        for key in ('BEST_Z', 'SAMPLE_METRICS'):
            self.assertEqual(resultsOld[key], resultsNew[key])
        newT = resultsNew['BEST_THRESHOLDS'] # threshold.txt path
        self.assertEqual(self.expectedT["8"], self.getMD5hex(newT))

    def test_call_binary(self):
        """Re-call GTC files using zCall with binary output"""
        self.call(True)

    def test_call_text(self):
        """Re-call GTC files using zCall with text output"""
        self.call(False)

    def test_complete(self):
        """Test self-contained zcall script"""
        zstart = 6
        ztotal = 5
        args = ['zcall/zCallComplete.py',
                '--bpm', self.bpmPath,
                '--egt', self.egtPath,
                '--out', self.outDir,
                '--zstart', str(zstart),
                '--ztotal', str(ztotal),
                '--samples', self.sampleJson,
                '--text', os.path.join(self.outDir, 'metric_summary.txt'),
                '--plink', self.prefix,
                '--binary'
                ]
        self.assertEqual(os.system(' '.join(args)), 0)
        outStem = os.path.join(self.outDir, self.prefix)
        suffixes = ['.bed', '.bim', '.fam']
        expected = ['8e222b46b0760cba5de1d2bded337c76',
                    '19dd8929cd63e3ee906e43b9bb59cd02',
                    'b836bd45459de6a9bc4b8d92e8c9e298']
        for i in range(len(suffixes)):
            self.assertTrue(os.path.exists(outStem+suffixes[i]))
            checksum = self.getMD5hex(outStem+suffixes[i])
            self.assertEqual(checksum, expected[i])  
        self.validatePlink(self.prefix, True)

if __name__ == "__main__":
    unittest.main(verbosity=2)
