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

"""Classes to find thresholds and evaluate z scores, as parameters for zCall.

Threshold finding combines: findMeanSD.py, findBetas.r, findThresholds.py
Evaluation looks for optimal concordance and gain metrics on given dataset.
Import classes to front-end scripts for calibration and calling.

Contents:
- ThresholdFinder
- MetricEvaluator
- MetricFinder
- SampleEvaluator

Author: Iain Bancarz, ib5@sanger.ac.uk, January 2013
"""

import json, os, re, sys, tempfile
from ConfigParser import ConfigParser
from GTC import *
from BPM import *
from EGT import *
from utilities import CallingBase, ThresholdContainer, SharedBase

class ThresholdFinder:
    """Class to write threshold.txt files for given EGT input and z score(s).

    Threshold finding procedure:
    1. Run findMeanSD on given EGT file.  Outputs mean_sd.txt
    2. Run findBetas.r on output from (1). Outputs betas.txt
    3. Run findThresholds on EGT file and betas.txt, with given Z score(s).
    Outputs from (1) and (2) are written to a temporary directory, deleted on successful exit.

    Recommended default Z score = 7.  Suggested range of alternatives = 3 to 15.
"""

#    TODO  Modify findMeanSD and findThresholds so they can be imported, instead of being run in a subshell
 
    INDEX_NAME = "thresholds.json"

    def __init__(self, egtPath, configPath=None):
        self.egtPath = egtPath
        self.egt = EGT(egtPath)
        if configPath==None:
            configPath = os.path.join(sys.path[0], '../etc/config.ini')
            configPath = os.path.abspath(configPath)
        self.snpTotal = self.egt.getTotalSNPs()
        config = ConfigParser()
        config.readfp(open(configPath))
        self.rScript = config.get('zcall', 'rscript')
        self.digits = int(config.get('zcall', 'digits'))

    def findMeanSD(self, outPath):
        """Find cluster mean/sd from EGT object

        Based on findMeanSD script in original zcall implementation
        By default, apply 'sanity check' filter to SNPs"""
        exclude = set()
        output = [None]*self.snpTotal
        maxExclusion = 0.95
        for i in range(self.snpTotal):
            # Get SNP Name
            snp = self.egt.names[i]
            # Get number of points in each genotype cluster
            nAA = self.egt.nAA[i]
            nAB = self.egt.nAB[i]
            nBB = self.egt.nBB[i]
            nTotal = nAA + nAB + nBB
            (meanX, meanY, devX, devY) = self.meanDevXY(i)
            if nAA >= nBB:
                out = [snp, meanX, meanY, devX, devY, nBB, nAA] 
            elif nBB > nAA:
                out = [snp, meanX, meanY, devX, devY, nAA, nBB] 
            out = [str(o) for o in out]
            output[i] = "\t".join(out)
            # apply sanity checks and flag for exclusion if appropriate
            if not self.snpSanityCheckOK(nAA, nAB, nBB, nTotal):
                exclude.add(i)
        excludeRate = len(exclude)/float(self.snpTotal)
        head = ["SNP", "meanX", "meanY", "sdX", "sdY", 
                "nMinorHom", "nCommonHom"] 
        outFile = open(outPath, 'w')
        outFile.write("\t".join(head)+"\n")
        if excludeRate <= maxExclusion:
            for i in range(self.snpTotal):
                if i not in exclude: outFile.write(output[i]+"\n")
        else:
            msg = "Warning: findMeanSD exclusion rate of "+str(excludeRate)+\
                " exceeds maximum of "+str(maxExclusion)+\
                "; omitting SNP filter.\n"
            sys.stderr.write(msg)
            for i in range(self.snpTotal):
                outFile.write(output[i]+"\n")
        outFile.close()

    def findThresholds(self, betasPath, zScore, outPath):
        """Find thresholds from findBetas.r output and write thresholds.txt

        Replaces findThresholds.py script in original zCall implementation"""
        (beta0, beta1) = self.readBetas(betasPath)
        output = [None]*(self.snpTotal+1)
        output[0] = ["SNP", "Tx", "Ty"]
        for i in range(self.snpTotal):
            # Get SNP Name
            snp = self.egt.names[i]
            (Tx, Ty) = self.findThresholdPair(i, zScore, beta0, beta1)
            out = [snp, Tx, Ty]
            out = [str(o) for o in out]
            output[i+1] = out
        outFile = open(outPath, 'w')
        for line in output: outFile.write("\t".join(line)+"\n")
        outFile.close()

    def findThresholdPair(self, i, z, beta0, beta1, minIntensity=0.2):
        """Find x,y thresholds for the ith snp"""
        nAA = self.egt.nAA[i]
        nBB = self.egt.nBB[i]

        # Extract mean and standard deviations from EGT class
        meanXAA = self.egt.meanXAA[i]
        meanXBB = self.egt.meanXBB[i]
        devXAA = self.egt.devXAA[i]
        devXBB = self.egt.devXBB[i]
        
        meanYAA = self.egt.meanYAA[i]
        meanYBB = self.egt.meanYBB[i]
        devYAA = self.egt.devYAA[i]
        devYBB = self.egt.devYBB[i]

        # Calculate Thresholds
        if nAA <= 2 and nBB <= 2: 
            # Not enough points in common allele homozygote cluster
            Tx = "NA"
            Ty = "NA"
        else:
            if nAA >= nBB: # AA cluster tags the common allele
                if meanXAA < minIntensity:
                    Tx = "NA" # skips new genotype calls in zCall
                    Ty = "NA"
                else:
                    Ty = meanYAA + z * devYAA
                    # Solve for mean, sd of the minor allele hom. cluster
                    # based on betas and mean of common allele hom. cluster
                    meanXBB = beta1[1]*meanYAA + beta0[1]
                    devXBB = beta1[3]*devYAA + beta0[3]
                    Tx = meanXBB + z * devXBB
            else: # BB cluster tags the common allele
                if meanYBB < minIntensity:
                    Tx = "NA"
                    Ty = "NA"
                else:
                    Tx = meanXBB + z * devXBB
                    # Solve for mean, sd as above
                    meanYAA = beta1[0] * meanXBB + beta0[0]
                    devYAA = beta1[2] * devXBB + beta0[2]
                    Ty = meanYAA + z * devYAA 
        if Tx!="NA": Tx = round(Tx, self.digits)
        if Ty!="NA": Ty = round(Ty, self.digits)
        return (Tx, Ty)

    def meanDevXY(self, i):
        """Find mean/sd for x,y for given snp"""
        # Extract the mean and sd for each common allele homozygote clusters in the noise dimension
        meanXAA = self.egt.meanXAA[i]
        meanXBB = self.egt.meanXBB[i]
        devXAA = self.egt.devXAA[i]
        devXBB = self.egt.devXBB[i]
        meanYAA = self.egt.meanYAA[i]
        meanYBB = self.egt.meanYBB[i]
        devYAA = self.egt.devYAA[i]
        devYBB = self.egt.devYBB[i]
        if meanXAA >= meanYAA: ## AA is in the lower right quadrant
            meanY = meanYAA
            devY = devYAA
            meanX = meanXBB
            devX = devXBB                
        else: ## AA is in the upper left quadrant; however, this should never be the case by definition X -> A, Y -> B
            meanY = meanYBB
            devY = devYBB
            meanX = meanXAA
            devX = devXAA
        return (meanX, meanY, devX, devY)
            
    def snpSanityCheckOK(self, nAA, nAB, nBB, nTotal):
        """Apply sanity checks to SNPs for inclusion in mean/SD calculation"""
        cr =  float(nTotal) / float(self.egt.numPoints)
        if cr < 0.99:
            return False
        # Make sure there are at least 10 points in each homozygote cluster
        if nAA < 10 or nBB < 10:
            return False
        # Calculate and check MAF
        if nAA > nBB:
            maf = (nAB + 2 * nBB) / float(2 * nTotal)
        else:
            maf = (nAB + 2 * nAA) / float(2 * nTotal)
        if maf < 0.05:
            return False
        # Hardy-Weinberg Equilibrium (don't use site if p_hwe < 0.00001)
        chiCritical = 19.5 # p = 0.00001 for 1 DOF
        if nAA > nBB:
            p = 1.0 - maf
            q = maf        
            expAA = p**2 * nTotal
            expAB = 2 * p * q * nTotal
            expBB = q**2 * nTotal
        else:
            p = 1.0 - maf
            q = maf        
            expAA = q**2 * nTotal
            expAB = 2 * p * q * nTotal
            expBB = p**2 * nTotal
        chiSquare = ((nAA - expAA)**2 / float(expAA)) + \
            ((nAB - expAB)**2 / float(expAB)) + \
            ((nBB - expBB)**2 / float(expBB))
        if chiSquare > chiCritical:
            return False
        return True

    def thresholdFileName(self, zScore):
        egtName = re.split('/', self.egtPath).pop()
        items = re.split('\.', egtName)
        items.pop() # remove .egt suffix
        name = '.'.join(items)
        return 'thresholds_'+name+'_z'+str(zScore).zfill(2)+'.txt'

    def readBetas(self, inPath):
        """Read output of findBetas.r"""
        beta0 = [] # list container for beta intercept
        beta1 = [] # list container for beta of slope
        lines = open(inPath).readlines()
        for i in range(1, len(lines)):
            fields = re.split("\t", lines[i].strip())
            beta0.append(float(fields[1]))
            beta1.append(float(fields[2]))
        return (beta0, beta1)
    
    def run(self, zScore=7, outDir='/tmp', verbose=True, force=False):
        outPath = os.path.join(outDir, self.thresholdFileName(zScore))
        if os.path.exists(outPath) and force==False:
            if verbose: print outPath+" already exists; omitting calibration."
            return outPath
        scriptDir = os.path.abspath(sys.path[0])
        tempDir = tempfile.mkdtemp(prefix='zcall_')
        if verbose:
            msg = "Calibrating zCall: zscore = "+str(zScore)+"\n"+\
                "Writing temporary files to "+tempDir
            print msg
        meanSd = tempDir+'/mean_sd.txt'
        betas = tempDir+'/betas.txt'
        self.findMeanSD(meanSd)
        # findBetas.r command uses bash to redirect stderr
        cmd = 'bash -c "'+self.rScript+' '+scriptDir+'/findBetas.r '+\
            meanSd+' '+betas+' 1 " &> '+tempDir+'/rscript.log'
        status = os.system(cmd)
        self.findThresholds(betas, zScore, outPath)
        if status==0: 
            if verbose: print "Cleaning up temporary directory."
            os.system('rm -Rf '+tempDir)
        elif verbose: 
            print "Possible error, retaining temporary directory."
        if verbose: print "Finished calibration."
        return outPath

    def runMultiple(self, zstart, ztotal, outDir, 
                    verbose=True, force=False):
        z = zstart
        thresholdPaths = {}
        for i in range(ztotal):
            thresholdPath = self.run(z, outDir, verbose, force)
            thresholdPaths[str(z)] = thresholdPath # .json needs string as key
            z += 1
        indexPath = os.path.join(outDir, self.INDEX_NAME)
        index = open(indexPath, 'w')
        index.write(json.dumps(thresholdPaths))
        index.close()
        return indexPath

class MetricEvaluator(SharedBase):
    """Class to assess concordance/gain metrics and choose best z score"""
    
    def __init__(self, configPath=None):
        if configPath==None:
            configPath = os.path.join(sys.path[0], '../etc/config.ini')
            configPath = os.path.abspath(configPath)
        config = ConfigParser()
        config.readfp(open(configPath))
        self.minConcord = float(config.get('zcall', 'min_concord'))

    def getBestThresholdKey(self):
        return self.T_KEY

    def findBestZ(self, concords, gains):
        """Find best z score from normalized mean concordance/gain values

        Normalize concordance/gain by dividing each value by the maximum. Optionally, set a minimum normalized concordance. The 'best' is defined as the smallest z s.t. (concordance > gain) AND (concordance > minimum); or if none exists, return z with maximum concordance.  Note that z keys in hash are strings (for use in .json output), but need to sort them in integer order."""
        if concords.keys() != gains.keys():
            raise ValueError("Concordance and gain inputs do not match")
        cMax = max(concords.values())
        gMax = max(gains.values())
        zScores = []
        for z in concords.keys():
            concords[z] = concords[z]/cMax
            gains[z] = gains[z]/gMax
            zScores.append(int(z))
        zScores.sort()
        best = None
        zMax = None
        for i in range(len(zScores)):
            z = str(zScores[i])
            if concords[z] > gains[z] and concords[z] > self.minConcord:
                best = z
                break
            elif concords[z]==1:
                zMax = z
        if best == None: best = zMax
        return str(best)

    def findMeans(self, inPaths, verbose=False, outPath=None):
        """Read JSON result paths, find mean concordance/gain by z score

        Return concatenation of results and concordance/gain metrics
        """
        rows = []
        for inPath in inPaths:
            rows.extend(json.loads(open(str(inPath)).read()))
        zCounts = {}
        concords = {}
        gains = {}
        zIntegers = set()
        for row in rows:
            [gtc, z, concord, gain, counts] = row
            try: 
                zCounts[z] += 1
                concords[z] += concord
                gains[z] += gain
            except KeyError: 
                zCounts[z] = 1
                concords[z] = concord
                gains[z] = gain
            zIntegers.add(int(z))
        zIntegers = list(zIntegers)
        zIntegers.sort()
        if verbose: print "z\tconcord\tgain"
        for zInt in zIntegers:
            z = str(zInt)
            concords[z] = concords[z] / float(zCounts[z])
            gains[z] = gains[z] / float(zCounts[z])
            if verbose: print z+"\t"+str(round(concords[z], 4))+"\t"+\
                    str(round(gains[z], 4))
        if outPath!=None:
            out = open(outPath, 'w')
            out.write(json.dumps(rows))
            out.close()
        return (rows, concords, gains)

    def writeBest(self, inPaths, thresholdPath, outPath, textPath=None,
                  verbose=False):
        """Find best z score & thresholds.txt, write to file for later use

        Arguments:
        - List of SampleEvaluator output paths
        - JSON file with hash of thresholds.txt paths by z score
        - Output directory
        """
        (metrics, concords, gains) = self.findMeans(inPaths, verbose)
        best = self.findBestZ(concords, gains)
        if verbose: print self.Z_KEY+":"+best
        thresholdPaths = json.loads(open(thresholdPath).read())
        results = { self.Z_KEY:best, self.T_KEY:thresholdPaths[best],
                    self.M_KEY:metrics}
        out = open(outPath, 'w')
        out.write(json.dumps(results))
        out.close()
        if textPath!=None:
            self.writeMeanText(concords, gains, textPath)
        return results

    def writeMeanText(self, concords, gains, outPath, digits=6):
        """Write plain text file with mean concordance/gain by z score"""
        zRange = []
        for z in concords.keys(): zRange.append(int(z))
        zRange.sort()
        output = []
        for zScore in zRange:
            z = str(zScore)
            concord = round(concords[z], digits)
            gain = round(gains[z], digits)
            output.append("%s\t%s\t%s\n" % (z, concord, gain) )
        out = open(outPath, 'w')
        out.write(''.join(output))
        out.close()

class MetricFinder(CallingBase):
    """Class to evaluate GTC objects by concordance and gain metrics

    Initialize with egt path, bpm path
    Inherits common "calling" functions from CallingBase
"""

    def concordanceRate(self, counts):
        """Find concordance rate between original and new calls

        Ignores SNPs where original was a 'no call'"""
        [match, total] = [0,0]
        for i in range(1,4):
            for j in range(0,4):
                count = counts[(i,j)]
                total += count
                if i==j: match += count
        concord = float(match)/float(total)
        return concord

    def countCallTypes(self, thresholds, gtc):
        """Count call types (0, AA, AB, BB) for given GTC and thresholds.

        call codes: 0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        based on method in sampleConcordance.py
        """
        self.setThresholds(thresholds)
        counts = {}
        for i in range(4):
            for j in range(4): counts[(i,j)] = 0
        includedSNPs = self.findIncludedSNPs(thresholds)
        for i in includedSNPs:
            nAA = self.egt.nAA[i]
            nBB = self.egt.nBB[i]
            nAB = self.egt.nAB[i]
            origCall = self.normalizeCall(gtc.genotypes[i], nAA, nBB)
            newCall = self.normalizeCall(self.call(gtc, i), nAA, nBB)
            counts[(origCall, newCall)] += 1
        return counts

    def getMetrics(self, thresholds, gtc):
        """Find call types, concordance and gain for given threshold and GTC

        Arguments are ThresholdContainer and GTC objects"""
        counts = self.countCallTypes(thresholds, gtc)
        concord = self.concordanceRate(counts)
        gain = self.gainRate(counts)
        return (counts, concord, gain)

    def gainRate(self, counts):
        """Find rate of call gain

        Defined as no calls in original which are called by zcall"""
        [gain, total] = [0,0]
        for i in range(4):
            count = counts[(0,i)] # no call in original GTC
            total += count
            if i!=0: gain += count
        gainRate = float(gain)/float(total)
        return gainRate

    def includeSNP(self, i, nAA, nBB, nAB, thresholds):
        """Should ith SNP be included in concordance calculation?

        Require autosomal SNP with MAF>=5%
        Want at least 10 points in each homozygote cluster
        Also exclude SNPs without defined zcall thresholds"""
        include = True
        chrom = self.bpm.chr[i]
        maf = self.findMAF(nAA, nBB, nAB)
        if maf < 0.05 or chrom == "X" or chrom == "Y" or nAA < 10 or nBB < 10 \
                or thresholds.getX(i)=="NA" or thresholds.getY(i)=="NA":
            include = False
        return include        

    def findIncludedSNPs(self, thresholds):
        """ Find set of included SNP indices """
        included = []
        for i in range(self.snpTotal):
            nAA = self.egt.nAA[i]
            nBB = self.egt.nBB[i]
            nAB = self.egt.nAB[i]
            if self.includeSNP(i, nAA, nBB, nAB, thresholds): 
                included.append(i)
        return included

class SampleEvaluator(SharedBase):
    """Evaluate z scores and thresholds for one or more GTC files."""

    def __init__(self, bpmPath, egtPath):
        self.bpmPath = bpmPath
        self.egtPath = egtPath
        self.bpm = BPM(bpmPath)
        self.metricFinder = MetricFinder(bpmPath, egtPath)
        self.evaluated = 0 # running total of evaluated samples
        self.sampleTotal = 0 # total number of samples to evaluate 

    def convertCountKeys(self, counts):
        """Convert keys in counts dictionary to string; required for JSON output

        Counts are indexed by (original call, final call)
        """
        output = {}
        for key in counts.keys():
            (i,j) = key
            output[str(i)+':'+str(j)] = counts[key]
        return output

    def evaluate(self, thresholds, gtc, verbose=False):
        """Evaluate z thresholds for given thresholds & sample GTC

        Inputs:
        - dictionary of threshold paths indexed by z score
        - GTC object
        """
        gtcName = os.path.split(gtc.getInputPath())[1]
        self.evaluated += 1
        if verbose: 
            print "Evaluating z scores for sample "+str(self.evaluated)+\
                " of "+str(self.sampleTotal)+": "+gtcName
        zList = thresholds.keys()
        zList.sort()
        results = {}
        for z in zList:
            if verbose: print "Finding metrics for z score", z
            results[z] = self.metricFinder.getMetrics(thresholds[z], gtc)
        output = []
        for z in zList:
            (counts, concord, gain) = results[z]           
            converted = self.convertCountKeys(counts)
            output.append([gtcName, z, concord, gain, converted])
        return output

    def run(self, thresholdPath, sampleJson, start, end, outPath,verbose=False):
        """Evaluate thresholds for given list of sample GTC paths

        Inputs:
        - Path to .json file with hash of paths to threshold.txt files
        - Path to .json file with paths of sample GTC files
        - Start index in GTC .json file (use to split GTC for batch processing)
        - End index in GTC .json file
        - Output path

        Output:
        - JSON file with GTC filename, z score, metrics, and call type counts

        GTC input is in sample JSON format used by genotyping pipeline """
        metrics = []
        if verbose: print "Reading thresholds."
        thresholdPaths = json.loads(open(thresholdPath).read())
        thresholds = {}
        for z in thresholdPaths.keys(): 
            thresholds[z] = ThresholdContainer(thresholdPaths[z])
        if verbose: print "Evaluating samples."
        output = []
        gtcPaths = self.readSampleJson(sampleJson)
        self.evaluated = 0
        if end==-1: end = len(gtcPaths)
        if start>=end:
            raise ValueError("Must have GTC start index < end index")
        elif start < 0:
            raise ValueError("Must have GTC start index > 0")
        elif end > len(gtcPaths):
            raise ValueError("Must have GTC end index <= total GTC paths")
        self.sampleTotal = end - start
        gtcPaths = gtcPaths[start:end]
        for gtcPath in gtcPaths:
            gtc = GTC(gtcPath, self.bpm.normID)
            output.extend(self.evaluate(thresholds, gtc, verbose))
        # write results
        out = open(outPath, 'w')
        out.write(json.dumps(output))
        out.close()
