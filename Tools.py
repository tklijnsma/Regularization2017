#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, tempfile, shutil, re, subprocess, sys, traceback
from os.path import *
from glob import glob

from time import strftime
datestr = strftime( '%b%d' )

import ROOT


########################################
# Plotting
########################################

ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gStyle.SetOptStat(0)
c = ROOT.TCanvas( 'ctc', 'ctc', 1000, 800 )

LeftMargin   = 0.15
RightMargin  = 0.03
BottomMargin = 0.15
TopMargin    = 0.03
def SetCMargins(
    LeftMargin   = 0.15,
    RightMargin  = 0.03,
    BottomMargin = 0.15,
    TopMargin    = 0.03,
    ):
    c.SetLeftMargin( LeftMargin )
    c.SetRightMargin( RightMargin )
    c.SetBottomMargin( BottomMargin )
    c.SetTopMargin( TopMargin )

PLOTDIR = 'plots_{0}'.format(datestr)
def SetPlotDir( newdir ):
    global PLOTDIR
    PLOTDIR = newdir
def SaveC( outname, asPNG=False, asROOT=False ):
    global PLOTDIR
    if not isdir(PLOTDIR): os.makedirs(PLOTDIR)
    subdir = ''
    if len(outname.split('/')) == 2:
        subdir = outname.split('/')[0]
        if not isdir( join( PLOTDIR, subdir ) ): os.makedirs( join( PLOTDIR, subdir ) )

    outname = join( PLOTDIR, subdir, basename(outname).replace('.pdf','').replace('.png','') )
    c.SaveAs( outname + '.pdf' )
    if asPNG:
        c.SaveAs( outname + '.png' )
    if asROOT:
        c.SaveAs( outname + '.root' )

ROOTCOUNTER = 1000
def GetUniqueRootName():
    global ROOTCOUNTER
    name = 'root{0}'.format(ROOTCOUNTER)
    ROOTCOUNTER += 1
    return name
def GetPlotBase(
        xMin = 0, xMax = 1,
        yMin = 0, yMax = 1,
        xTitle = 'x', yTitle = 'y',
        SetTitleSizes = True,
    ):
    base = ROOT.TH1F()
    ROOT.SetOwnership( base, False )
    base.SetName( GetUniqueRootName() )
    base.GetXaxis().SetLimits( xMin, xMax )
    base.SetMinimum( yMin )
    base.SetMaximum( yMax )
    base.SetMarkerColor(0)
    base.GetXaxis().SetTitle( xTitle )
    base.GetYaxis().SetTitle( yTitle )
    if SetTitleSizes:
        base.GetXaxis().SetTitleSize( 0.06 )
        base.GetYaxis().SetTitleSize( 0.06 )
    return base


########################################
# Help functions
########################################

TESTMODE = False
def TestMode( flag=True ):
    global TESTMODE
    TESTMODE = flag
def IsTestMode():
    global TESTMODE
    return TESTMODE

DEFAULTJOBDIR = abspath( 'Scan_{0}'.format(datestr) )
def SetDefaultJobDir( newdirname='Scan_{0}'.format(datestr) ):
    global DEFAULTJOBDIR
    DEFAULTJOBDIR = newdirname
def GetDefaultJobDir():
    global DEFAULTJOBDIR
    return DEFAULTJOBDIR


def executeCommand( cmd, captureOutput=False ):
    if not isinstance( cmd, basestring ):
        cmdStr = '\n    '.join( cmd )
        cmdExec = ' '.join(cmd)
    else:
        cmdStr = cmd
        cmdExec = cmd

    if TESTMODE:
        print '\nTESTMODE: ' + cmdStr + '\n'
    else:
        if not captureOutput:
            print '\nEXECUTING: ' + cmdStr + '\n'
            os.system( cmdExec )
        else:
            output = subprocess.check_output(
                cmd,
                shell=True,
                )
            return output


class AnalysisError(Exception):
    pass
def ThrowError( errstr ):
    raise AnalysisError( errstr )


########################################
# Physics functions
########################################

def BasicT2WS(
        datacard,
        maps         = None,
        options      = None,
        smartMaps    = None,
        outputWSfile = None,
        verbose      = False,
        ):

    outputDir = abspath( 'workspaces_{0}'.format(datestr) )
    if not isdir( outputDir ): os.makedirs( outputDir )
    
    if not outputWSfile: outputWSfile = basename(datacard).replace( '.txt', '.root' )
    outputWSfile = join( outputDir, outputWSfile )

    cmd = []
    cmd.append( 'text2workspace.py' )
    cmd.append( datacard )
    cmd.append( '-o {0}'.format(outputWSfile) )
    cmd.append( '-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel' )
    cmd.append( '--PO verbose' )
    cmd.append( '--PO \'higgsMassRange=123,127\'' )

    if smartMaps:

        # Example procPat:
        # '--PO \'map=.*/InsideAcceptance_genPt_200p0_350p0:r_xH_PTH_200_350[1.0,-1.0,4.0]\'',
        # Want to replace this by:
        # '--PO \'map=.*/InsideAcceptance_genPt_([\dpm]+)p0_([\dpm]+)p0:r_xH_PTH_\1_\2[1.0,-1.0,4.0]\'',

        # in smartMap form: ( r'.*/InsideAcceptance_genPt_([\dm]+)p0_([\dm]+)p0', r'r_xH_PTH_\1_\2[1.0,-1.0,4.0]' )

        # Manual maps should override smart maps; gather all the patters that are already in a manualMap
        if maps:
            manualMapPats = []
            for manualMap in maps:
                match = re.search( r'map=(.*):', manualMap )
                if not match: continue
                manualMapPats.append( match.group(1) )

        newMaps = []
        for binprocPat, yieldParPat in smartMaps:

            for proc in signalprocesses:
                for bin in bins:

                    binprocStr = '{0}/{1}'.format( bin, proc )

                    manualMapAvailable = False
                    if maps:
                        for manualMapPat in manualMapPats:
                            if re.match( manualMapPat, binprocStr ):
                                manualMapAvailable = True

                    if manualMapAvailable:
                        continue
                    elif not re.match( binprocPat, binprocStr ):
                        continue

                    if verbose: print 'Pattern \"{0}\" matches with \"{1}\"'.format( binprocPat, binprocStr )

                    yieldPar = re.sub( binprocPat, yieldParPat, binprocStr )
                    newMap = '--PO \'map={0}:{1}\''.format( binprocStr, yieldPar )

                    newMaps.append( newMap )

        for newMap in newMaps:
            cmd.append(newMap)

    if maps:
        for manualMap in maps:
            cmd.append( manualMap )

    if not options:
        pass
    elif isinstance( options, basestring ):
        cmd.append( options )
    else:
        cmd.extend( options )

    executeCommand( cmd )




def MultiDimCombineTool(
        datacard,
        nPoints       = 100,
        nPointsPerJob = 3,
        queue         = '1nh',
        notOnBatch    = False,
        jobDirectory  = None,
        fastscan      = False,
        asimov        = False,
        jobPriority   = 0,
        extraOptions  = [],
        ):

    datacard = abspath( datacard )

    scanName = basename(datacard).replace('.root','')
    scanName = re.sub( r'\W', '', scanName )
    scanName = 'SCAN_{0}_{1}'.format( datestr, scanName )

    currentdir = os.getcwd()
    if not TESTMODE:
        if not jobDirectory:
            jobDirectory = DEFAULTJOBDIR
        jobDirectory = AppendNumberToDirNameUntilItDoesNotExistAnymore( jobDirectory )
        if not isdir( jobDirectory ):
            print 'Creating directory {0}'.format( jobDirectory )
            os.makedirs( jobDirectory )
        print 'Moving to directory {0}'.format( jobDirectory )
        os.chdir( jobDirectory )


    cmd = [
        'combineTool.py',
        datacard,
        '-n {0}'.format( scanName ),
        '-M MultiDimFit',
        '--cminDefaultMinimizerType Minuit2',
        '--cminDefaultMinimizerAlgo migrad',
        '--algo=grid',
        '--saveNLL',
        '--saveInactivePOI 1',
        '--points={0} '.format(nPoints),
        '--minimizerStrategy 2',
        ]

    if fastscan:
        cmd.append( '--fastScan' )
    if asimov:
        cmd.append( '-t -1' )

    if not extraOptions:
        pass
    elif isinstance( extraOptions, basestring ):
        cmd.append( extraOptions )
    else:
        cmd.extend( extraOptions )


    if not notOnBatch:
        cmd.append(
            '--split-points {0} '.format(nPointsPerJob)
            )

        if 't3' in os.environ['HOSTNAME']:

            if not queue in [ 'all.q', 'long.q', 'short.q' ]:
                print 'Queue \'{0}\' is not available on PSI'.format(queue)
                return

            if jobPriority != 0:
                cmd.append(
                    '--job-mode psi --task-name {0} --sub-opts=\'-q {1} -p {2}\' '.format( scanName, queue, jobPriority ),
                    )
            else:
                cmd.append(
                    '--job-mode psi --task-name {0} --sub-opts=\'-q {1}\' '.format( scanName, queue ),
                    )

        else:

            if not queue in [ '8nm', '1nh', '8nh', '1nd', '2nd' ]:
                print 'Queue \'{0}\' is not available on lxplus'.format(queue)
                if not IsTestMode():
                    return

            cmd.append(
                '--job-mode lxbatch --task-name {0} --sub-opts=\'-q {1}\' '.format( scanName, queue ),
                )

    executeCommand( cmd )

    os.chdir( currentdir )



def ListOfPDFIndicesToFreeze( postfitFile, freezeAllIndices=False, verbose=False, snapshotName='MultiDimFit', returnAlsoFloating=False ):

    wsFp = ROOT.TFile.Open( postfitFile )
    ws   = wsFp.Get('w')
    loadedSnapshot = ws.loadSnapshot(snapshotName)

    if not loadedSnapshot:
        ThrowError( 'Could not load {0} snapshot - are you passing the post fit workspace?'.format( snapshotName ), throwException=True )

    varsToFreeze = []
    varsToFloat = []
    
    # Find all pdf indexes 
    allCats = ws.allCats()
    catitr = allCats.createIterator()
    cat = catitr.Next()
    while cat:
        if cat.GetName().startswith("pdfindex"):
            varsToFreeze.append( cat.GetName() )
        cat = catitr.Next()
        
    # Find all background pdfs
    allPdfs = ws.allPdfs()

    if verbose:
        nPdfs = allPdfs.getSize()
        print 'Looping over {0} pdfs in the workspace'.format(nPdfs)

    pdfitr = allPdfs.createIterator()
    pdf = pdfitr.Next()
    while pdf:
        if pdf.GetName().startswith("shapeBkg_bkg"):
            # bgks from hzz are RooHistPdfs, not RooMultiPdfs
            if hasattr( pdf, 'getNumPdfs' ):
                # Loop over all shapes in the envelope
                for ishape in xrange(pdf.getNumPdfs()):

                    shape = pdf.getPdf( ishape )
                    observables = ROOT.RooArgList(shape.getObservables( ws.allVars() ))
                    observables = filter(lambda x: not x.startswith('CMS_hgg_mass'), map(lambda x: observables[x].GetName(), xrange(observables.getSize()) ) )

                    # Freeze all pdf parameters except those from the best fit function
                    if ishape == pdf.getCurrentIndex() and not freezeAllIndices:
                        varsToFloat.extend( observables )
                    else:
                        varsToFreeze.extend( observables )
        pdf = pdfitr.Next()

    # For some reason, the first variable is never frozen; simply append it again at end of list
    varsToFreeze.append( varsToFreeze[0] )

    if verbose:
        print '\n\nFreezing the following variables:'
        for i in varsToFreeze:
            print i

        print '\n\nFloating the following variables:'
        for i in varsToFloat:
            print i

    wsFp.Close()

    if returnAlsoFloating:
        return varsToFreeze, varsToFloat
    else:
        return varsToFreeze


def ListSet(
        datacardRootFile,
        setName='POI',
        pattern='*',
        ):

    closeFp = False
    if isinstance( datacardRootFile, basestring ):
        datacardFp = ROOT.TFile.Open( datacardRootFile )
        w = datacardFp.Get('w')
        closeFp = True
    elif isinstance( datacardRootFile, ROOT.RooWorkspace ):
        w = datacardRootFile

    argset = w.set(setName)
    if not argset:
        print 'No set \'{0}\' in {1}'.format( setName, datacardRootFile )
        datacardFp.Close()
        return []

    arglist = ROOT.RooArgList( argset )

    varNames = []
    for i in xrange( arglist.getSize() ):
        varName = arglist[i].GetName()
        if pattern != '*':
            if re.search( pattern, varName ):
                varNames.append( varName )
        else:
            varNames.append( varName )

    if closeFp: datacardFp.Close()
    return varNames




def GetSignalIndices( dc, verbose=False ):

    if isfile( dc ):
        dcFile = dc
        with open( dcFile, 'r' ) as dcFp:
            dc = dcFp.read()

    # Find lines that start with 'process'
    lines = dc.split('\n')

    foundFirstProcessLine  = False
    foundSecondProcessLine = False
    for line in lines:
        if line.startswith('process '):
            if not foundFirstProcessLine:
                firstProcessLine = line
                foundFirstProcessLine = True
            elif foundFirstProcessLine and not foundSecondProcessLine:
                secondProcessLine = line
                break
    else:
        raise RuntimeError( 'Could not find two lines that start with \'process \' in the supplied datacard' )

    processNames   = firstProcessLine.split()[1:]
    processIndices = [ int(i) for i in secondProcessLine.split()[1:] ]

    signalDict = {}
    for name, index in zip( processNames, processIndices ):
        if index >= 0:
            # Skip backgrounds
            continue
        elif 'OutsideAcceptance' in name:
            # Skip outside acceptance
            continue    
        else:
            signalDict[index] = name

    if verbose:
        print '\nFound the following signals:'
        for key, value in signalDict.iteritems():
            print '  {0:4}: {1}'.format( key, value )

    return signalDict



def ConvertTChainToArray(
        rootFileList,
        treeName        = 'limit',
        variablePattern = '*',
        returnStyle     = 'dictPerPoint',
        verbose         = False
        ):

    if len(rootFileList) == 0:
        ThrowError( 'rootFileList has length 0' )
        sys.exit()

    if verbose: 'Looking for at least one filled root file to obtain the variable list from...'
    foundTree = False
    for rootFile in rootFileList:

        rootFp = ROOT.TFile.Open( rootFile )
        allKeys = rootFp.GetListOfKeys()

        if not allKeys.Contains( 'limit' ):
            if verbose: print '    No tree \'{0}\' in \'{1}\''.format( treeName, rootFile )
            rootFp.Close()
        else:
            tree = rootFp.Get( treeName )
            allVarObjArray = tree.GetListOfBranches()
            treeLoaded = True
            if verbose: print '    Found tree in {0}'.format( rootFile )
            break
    else:
        ThrowError( 'Not a single root file had a tree called \'{0}\'; Cannot extract any data', throwException=True )


    # Get the variable list from this one file
    nAllVars = allVarObjArray.GetEntries()

    useVars = []
    for iVar in xrange(nAllVars):
        varName = allVarObjArray[iVar].GetTitle()
        if not variablePattern == '*':
            if re.search( variablePattern, varName ):
                useVars.append( varName.split('/',1)[0] )
        else:
            useVars.append( varName.split('/',1)[0] )

    rootFp.Close()


    # Now read the entries from the chain
    chain = ROOT.TChain( treeName )
    for rootFile in rootFileList:
        chain.Add( rootFile )


    # Return an object

    if returnStyle == 'dict':

        res = {}
        for varName in useVars:
            res[varName] = []

        for event in chain:
            for varName in useVars:
                res[varName].append( getattr( event, varName ) )


    elif returnStyle == 'dictPerPoint':

        res = []
        for event in chain:
            entry = {}
            for varName in useVars:
                entry[varName] = getattr( event, varName )
            res.append( entry )

    return res



def AppendNumberToDirNameUntilItDoesNotExistAnymore(
        dirName,
        nAttempts = 100,
        ):

    dirName = abspath(dirName)

    if not isdir(dirName):
        return dirName

    dirName += '_{0}'
    for iAttempt in xrange(nAttempts):
        if not isdir( dirName.format(iAttempt) ):
            dirName = dirName.format(iAttempt)
            break
    else:
        ThrowError( 'Could not create a unique directory for {0}'.format(dirName.format('X')) )
        sys.exit()

    print '[info] New directory: {0}'.format( dirName )
    return dirName