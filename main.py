#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, argparse, re, sys, shutil
from os.path import *
from glob import glob
from array import array

from time import strftime
datestr = strftime( '%b%d' )
datestr_detailed = strftime( '%y%m%d_%H%M%S' )

import ROOT
import Tools
from Tools import c
from Tools import SetCMargins
from Tools import SetPlotDir
from Tools import SaveC
from Tools import GetUniqueRootName
from Tools import GetPlotBase


########################################
# Main
########################################

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument( '-v', '--verbose',        action='store_true' )
    parser.add_argument( '-t', '--test',           action='store_true' )
    parser.add_argument( '-i', '--input',     type=str )
    parser.add_argument( '--postfitWS',       type=str )
    parser.add_argument( '-b', '--onBatch',        action='store_true' )

    parser.add_argument( '--processTextDatacards', action='store_true')

    parser.add_argument( '--t2ws',                 action='store_true')
    parser.add_argument( '--makePostfit',          action='store_true')
    parser.add_argument( '--tauscan',              action='store_true')
    parser.add_argument( '--PlotTauScan',          action='store_true')

    parser.add_argument( '--tauscanToys',          action='store_true')
    parser.add_argument( '--makeToys',             action='store_true')
    parser.add_argument( '--fitToys',              action='store_true')
    parser.add_argument( '--cleanToys',            action='store_true')
    parser.add_argument( '--plotToys',             action='store_true')

    args = parser.parse_args()    

    if args.test: Tools.TestMode()

    # Overwrite for now
    args.verbose = True


    #____________________________________________________________________
    if args.processTextDatacards:

        dcFile = args.input
        with open( dcFile, 'r' ) as dcFp:
            dc = dcFp.read()

        nBins = len(Tools.GetSignalIndices( dc, verbose=args.verbose ))

        constraintLines = []

        # First constraint is always the same:
        constraintLines.append( 'constr0 constr RooFormulaVar::fconstr0("sqrt(2.0*tau)*(-r0+r1)",{tau,r0,r1}),const_mean[0.],const_std[1.]' )

        # Second derivatives in middle of spectrum
        for iBin in xrange(1,nBins-1):
            constraintLines.append(
                'constr{1} constr RooFormulaVar::fconstr{1}("sqrt(2.0*tau)*(r{0}-2*r{1}+r{2})",{{tau,r{0},r{1},r{2}}}),const_mean[0.],const_std[1.]'.format(
                    iBin-1, iBin, iBin+1
                    ))

        # Final bin slightly different again
        constraintLines.append(
            'constr{0} constr RooFormulaVar::fconstr{0}("sqrt(2.0*tau)*(r{1}-r{0})",{{tau,r{1},r{0}}}),const_mean[0.],const_std[1.]'.format(
                nBins-1, nBins-2
                ))

        if args.verbose:
            print '\nAdding the following lines to the text datacard:'
            print '\n'.join(constraintLines)
        dc += '\n\n' + '\n'.join(constraintLines)

        # wsDir = 'workspaces_{0}'.format(datestr)
        # if not isdir( wsDir ): os.makedirs(wsDir)
        # out = join( wsDir, basename(dcFile).replace( '.txt', '_regularised.txt' ) )
        out = dcFile.replace( '.txt', '_regularised.txt' )

        if args.verbose: print '\nWriting edited datacard to {0}'.format(out)
        with open( out, 'w' ) as outFp:
            outFp.write(dc)


    #____________________________________________________________________
    if args.t2ws:

        dcFile = args.input

        signalDict = Tools.GetSignalIndices( dcFile, verbose=args.verbose )

        sortedIndices = signalDict.keys()
        sortedIndices.sort( reverse=True )
        sortedSignals = [ signalDict[i] for i in sortedIndices ]

        if args.verbose:
            print '\nAssuming the following order of the signals:\n  ' + '\n  '.join( sortedSignals )

        maps = []

        # Add dummy map to create the variable OneOverSqrtTau
        maps.append(
            '--PO \'map=dummyTau:tau[0.0]\''
            )

        for iSignal, signal in enumerate(sortedSignals):
            maps.append(
                '--PO \'map=.*/{0}:r{1}[1.0,-1.0,4.0]\''.format( signal, iSignal )
                )

        Tools.BasicT2WS(
            dcFile,
            maps    = maps,
            verbose = args.verbose,
            )


    #____________________________________________________________________
    if args.makePostfit:

        ws              = abspath(args.input)
        wsTag           = basename(ws).replace('/','').replace('.root','')
        postfitFilename = join( os.getcwd(), 'higgsCombine_POSTFIT_{0}.MultiDimFit.mH125.root'.format( wsTag ) )

        cmd = [
            'combine',
            ws,
            '-M MultiDimFit',
            '--saveNLL',
            '--minimizerStrategy 2',
            '-v 2',
            '--floatOtherPOIs=1',
            # 
            '--saveWorkspace',
            '-m 125',
            '--setPhysicsModelParameters tau=0.0',
            '--freezeNuisances tau',
            '-n _POSTFIT_{0}'.format( wsTag ),
            # 
            '; mv {0} {1}'.format( postfitFilename, ws.replace( '.root', '_POSTFIT.root' ) )
            ]

        Tools.executeCommand( cmd, onBatch = args.onBatch )



    #____________________________________________________________________

    if args.tauscan or args.tauscanToys:

        # Define a tau scan function
        def doTauscan(
                ws,
                toyFile           = None,
                onBatch           = True,
                freezePdfIndices  = True,
                nPoints           = 39,
                nPointsPerJob     = 3,
                queue             = 'short.q',
                tauRange          = [ 0.0, 5.0 ],
                jobDirectory      = 'Scan_{0}_tauscan'.format(datestr),
                POIs              = None,
                skipInitialFit    = False,
                ):

            wsTag           = basename(ws).replace('/','').replace('.root','')
            # postfitFilename = join( os.getcwd(), 'higgsCombine_POSTFIT_{0}.MultiDimFit.mH125.root'.format( wsTag ) )

            # ======================================
            # Parse command

            if POIs is None:
                POIs = Tools.ListSet( ws, 'POI' )
                POIs.append( 'tau' )
            POIrangeStr = ':'.join([ '{0}=-1.0,4.0'.format(POI) for POI in POIs if re.match( r'r\d+', POI ) ])

            extraOptions = [
                '-P tau',
                '--setPhysicsModelParameterRanges tau={0},{1}'.format( *tauRange ) + ':' + POIrangeStr,
                # '--saveSpecifiedNuis tau',
                # 
                # '--floatOtherPOIs=1',
                '-m 125.00',
                '--snapshotName MultiDimFit',
                # '--redefineSignalPOIs {0}'.format( ','.join(POIs) ),
                ]

            if toyFile:
                extraOptions.append( '--dataset {0}:toys/toy_1 '.format( toyFile ) )

            if skipInitialFit:
                extraOptions.append('--skipInitialFit')

            if freezePdfIndices:
                pdfIndicesToFreeze = Tools.ListOfPDFIndicesToFreeze(
                    ws,
                    # verbose=args.verbose
                    )
                extraOptions.append(
                    '--freezeNuisances {0}'.format( ','.join(pdfIndicesToFreeze) )
                    )
                extraOptions.append( '--floatNuisances {0}'.format( ','.join(POIs) ) )
            else:
                extraOptions.append( '--floatOtherPOIs=1' )

            Tools.MultiDimCombineTool(
                ws,
                nPoints       = nPoints,
                nPointsPerJob = nPointsPerJob,
                queue         = queue,
                notOnBatch    = (not onBatch),
                jobDirectory  = jobDirectory,
                fastscan      = False,
                asimov        = False,
                jobPriority   = 0, # Highest; set to negative for lower
                extraOptions  = extraOptions,
                )


    if args.tauscan:

        ws = abspath(args.input)

        # ======================================
        # Flags and options

        FREEZE_PDF_INDICES = True
        # FREEZE_PDF_INDICES = False

        nPoints        = 39
        nPointsPerJob  = 3
        queue          = 'short.q'
        tauRange       = [ 0.0, 5.0 ]

        # # Debugging
        # nPoints        = 2
        # nPointsPerJob  = 2
        # queue          = 'short.q'
        # tauRange       = [ 0.0, 10.0 ]

        doTauscan(
            ws,
            onBatch           = False,
            freezePdfIndices  = True,
            nPoints           = 39,
            nPointsPerJob     = 3,
            queue             = 'short.q',
            tauRange          = [ 0.0, 5.0 ],
            POIs              = None,
            )


    #____________________________________________________________________
    def PlotTauScan(
            rootfiles,
            postfitWS = None,
            plotname  = 'tauscan',
            ):

        scanPoints = Tools.ConvertTChainToArray(
            rootfiles,
            treeName        = 'limit',
            variablePattern = '*',
            returnStyle     = 'dictPerPoint',
            verbose         = args.verbose
            )
        scanPoints.sort( key = lambda point: point['tau'] )

        mus = [ key for key in scanPoints[0].keys() if re.match( r'r\d+', key ) ]
        mus.sort( key = lambda mu: int(mu.replace('r','')) )


        if args.verbose:
            print '\nFound following scan points:'
            for scanPoint in scanPoints:
                line = []
                line.append( 'tau = {0:7.2f}'.format(scanPoint['tau']) )
                for mu in mus:
                    line.append( '{0} = {1:5.2f}'.format( mu, scanPoint[mu] ) )
                line.append( 'deltaNLL = {0:7.2f}'.format(scanPoint['deltaNLL']) )
                print ' | '.join(line)
            print

        nMus = len(mus)
        xAxis_binBoundaries = range(nMus+1)


        bestfitFound = False
        if postfitWS:
            bestfitPoint = Tools.ConvertTChainToArray(
                [ postfitWS ],
                treeName        = 'limit',
                variablePattern = '*',
                returnStyle     = 'dict',
                verbose         = args.verbose
                )
            bestfitPoint = { key : value[0] for key, value in bestfitPoint.iteritems() }
            bestfitFound = True

            if args.verbose:
                print 'Found unregularized best fit:'
                line = []
                line.append( 'tau = {0:7.2f}'.format(bestfitPoint['tau']) )
                for mu in mus:
                    line.append( '{0} = {1:5.2f}'.format( mu, bestfitPoint[mu] ) )
                line.append( 'deltaNLL = {0:7.2f}'.format(bestfitPoint['deltaNLL']) )
                print ' | '.join(line) + '\n'


        # ======================================
        # Plot

        c.Clear()
        SetCMargins()

        base = GetPlotBase(
            xMin = 0, xMax = nMus,
            yMin = -.5, yMax = 2.,
            xTitle = 'Bin_{i}', yTitle = '#mu'
            )
        base.Draw('P')


        minLineAlpha = 0.2
        maxLineAlpha = 0.9
        axisLineAlpha = [ minLineAlpha + ((float(i)/(len(scanPoints)-1))**30) * (maxLineAlpha-minLineAlpha) for i in xrange(len(scanPoints)) ]
        axisLineAlpha = list(reversed(axisLineAlpha))

        for iScanPoint, scanPoint in enumerate(scanPoints):

            H = ROOT.TH1F(
                GetUniqueRootName(), '',
                nMus,
                array( 'f', xAxis_binBoundaries )
                )
            ROOT.SetOwnership( H, False )

            for iMu, mu in enumerate(mus):
                H.SetBinContent( iMu+1, scanPoint[mu] )

            H.SetLineWidth(2)
            H.SetLineColorAlpha( 2, axisLineAlpha[iScanPoint] )
            H.Draw('HISTSAME')

            # if iScanPoint > 7: break

        if bestfitFound:

            H = ROOT.TH1F(
                GetUniqueRootName(), '',
                nMus,
                array( 'f', xAxis_binBoundaries )
                )
            ROOT.SetOwnership( H, False )

            for iMu, mu in enumerate(mus):
                H.SetBinContent( iMu+1, bestfitPoint[mu] )

            H.SetLineWidth(2)
            H.SetLineColor(4)
            H.Draw('HISTSAME')

        SaveC( plotname )


    #____________________________________________________________________
    if args.PlotTauScan:
        scandir = args.input
        PlotTauScan(
            glob( '{0}/*.root'.format(scandir) ),
            postfitWS = args.postfitWS,
            )


    #____________________________________________________________________
    if args.makeToys:

        # ======================================
        # Flags and options

        STARTSEED = 1000

        NTOYS = 1000

        postfitWS = abspath(args.postfitWS)


        # ======================================
        # 

        overDir = abspath('toys_{0}'.format(datestr))
        overDir = Tools.AppendNumberToDirNameUntilItDoesNotExistAnymore( overDir )

        currentdir = os.getcwd()
        if not Tools.IsTestMode():
            print 'Creating directory {0}'.format( overDir )
            os.makedirs( overDir )
            print 'Moving to directory {0}'.format( overDir )
            os.chdir( overDir )
        else:
            print 'Would now mkdir & cd \'{0}\''.format( overDir )


        POIs = [ p for p in Tools.ListSet( postfitWS, 'POI' ) if re.match( r'r\d+', p ) ]


        for iToy in xrange(NTOYS):
            if not Tools.IsTestMode(): os.chdir(overDir)

            seed = STARTSEED + iToy

            toyDir = 'toy_{0}_{1}'.format( iToy, seed )
            if not Tools.IsTestMode():
                print 'Creating directory {0}'.format( toyDir )
                os.makedirs( toyDir )
                print 'Moving to directory {0}'.format( toyDir )
                os.chdir( toyDir )
            else:
                print 'Would now mkdir & cd \'{0}\''.format( toyDir )

            toyFile = abspath( 'higgsCombine_toy.GenerateOnly.mH125.{0}.root'.format( seed ) )

            cmd = [
                'combine',
                postfitWS,
                '--seed {0}'.format(seed),
                '-M GenerateOnly',
                '--snapshotName MultiDimFit',
                '--saveToys',
                '-t 1',
                '--toysFrequentist --bypassFrequentistFit',
                '-m 125',
                '--setPhysicsModelParameters {0}'.format( ','.join([ '{0}=1.0'.format(p) for p in POIs ]) ),
                # 
                '-n _toy',
                ]

            Tools.executeCommand( cmd )

        os.chdir(currentdir)


    #____________________________________________________________________
    if args.fitToys:
        currentdir = os.getcwd()

        overDir   = abspath(args.input)
        postfitWS = abspath(args.postfitWS)

        POIs = [ p for p in Tools.ListSet( postfitWS, 'POI' ) if re.match( r'r\d+', p ) ]

        os.chdir( overDir )
        toyDirs = glob( 'toy_*' )

        if len(toyDirs) == 0:
            Tools.ThrowError( 'Could not find any \'toy_*\' directories in {0}'.format(overDir) )

        for toyDir in toyDirs:
            os.chdir( toyDir )

            toyFile =  glob( 'higgsCombine_toy.GenerateOnly.mH125.*.root' )[0]
            postfitToyFile = 'higgsCombine_toyPostfit.MultiDimFit.mH125.root'

            cmd = [
                'combine',
                '{0}'.format( postfitWS ),
                '-M MultiDimFit',
                '--dataset {0}:toys/toy_1 '.format( toyFile ),
                '-m 125',
                '--floatOtherPOIs=1',
                '--saveWorkspace',
                '--minimizerStrategy 2',
                '--setPhysicsModelParameters tau=0.0,{0}'.format( ','.join([ '{0}=1.0'.format(p) for p in POIs ]) ),
                '--freezeNuisances tau',
                # 
                '-n _toyPostfit',
                ]

            Tools.executeCommand( cmd, onBatch=True )                

            os.chdir( overDir )
        os.chdir(currentdir)


    #____________________________________________________________________
    if args.tauscanToys:

        # ======================================
        # Flags and options

        # REQUIRE_ALL_TOYS = True
        REQUIRE_ALL_TOYS = False

        # ONLY_TOYS_MISSING_SCAN = True
        ONLY_TOYS_MISSING_SCAN = False

        # FREEZE_PDF_INDICES = True
        FREEZE_PDF_INDICES = False

        ON_BATCH = True
        # ON_BATCH = False

        nPoints        = 120
        tauRange       = [ 0.0, 25.0 ]

        nPointsPerJob  = nPoints
        # queue          = 'short.q'
        queue          = 'all.q'

        # ======================================
        # Submit tau scan jobs

        currentdir = os.getcwd()
        overDir    = abspath(args.input)

        os.chdir( overDir )
        toyDirs = glob( 'toy_*' )

        if len(toyDirs) == 0:
            Tools.ThrowError( 'Could not find any \'toy_*\' directories in {0}'.format(overDir) )

        if REQUIRE_ALL_TOYS:
            for toyDir in toyDirs:
                if not isfile(join( toyDir, 'higgsCombine_toyPostfit.MultiDimFit.mH125.root' )):
                    Tools.ThrowError( 'File \'higgsCombine_toyPostfit.MultiDimFit.mH125.root\' does not exist in \'{0}\', but is required by REQUIRE_ALL_TOYS'.format( toyDir ) )
                POIs = Tools.ListSet( join( toyDir, 'higgsCombine_toyPostfit.MultiDimFit.mH125.root' ) )

        for toyDir in toyDirs:

            # print toyDir
            # if not 'toy_1_1001' in basename(toyDir):
            #     continue

            os.chdir( toyDir )

            if ONLY_TOYS_MISSING_SCAN and len(glob('tauscan_*')) > 0:
                print 'Skipping {0} (Only doing toys missing a tauscan directory)'.format(toyDir)
                os.chdir( overDir )
                continue

            postfitToyFile = abspath(       'higgsCombine_toyPostfit.MultiDimFit.mH125.root' )
            toyFile        = abspath( glob( 'higgsCombine_toy.GenerateOnly.mH125.*.root' )[0] )

            try:
                POIs = Tools.ListSet( postfitToyFile )
            except Tools.AnalysisError as error:
                if REQUIRE_ALL_TOYS:
                    # Should never happen; error should already be raised above
                    raise
                else:
                    print 'Skipping {0} (No workspace \'w\' found in postfit)'.format(toyDir)
                    os.chdir( overDir )
                    continue

            doTauscan(
                postfitToyFile,
                toyFile,
                onBatch           = ON_BATCH,
                freezePdfIndices  = FREEZE_PDF_INDICES,
                nPoints           = nPoints,
                nPointsPerJob     = nPointsPerJob,
                queue             = queue,
                tauRange          = tauRange,
                POIs              = None,
                jobDirectory      = 'tauscan_{0}'.format( datestr_detailed ),
                skipInitialFit    = True,
                )

            os.chdir( overDir )
        os.chdir(currentdir)
    

    #____________________________________________________________________
    if args.cleanToys:

        currentdir = os.getcwd()
        overDir    = abspath(args.input)


        CLEAN_JOBS = True
        # CLEAN_JOBS = False

        CLEAN_TAUSCANS = True
        # CLEAN_TAUSCANS = False

        # CLEAN_POSTFITS = True
        CLEAN_POSTFITS = False


        os.chdir( overDir )
        toyDirs = glob( 'toy_*' )

        if len(toyDirs) == 0:
            Tools.ThrowError( 'Could not find any \'toy_*\' directories in {0}'.format(overDir) )

        for toyDir in toyDirs:
            os.chdir( toyDir )

            filesToDelete = []

            if CLEAN_JOBS:
                filesToDelete.extend(glob( 'job*.sh*' ))

            if CLEAN_TAUSCANS:
                filesToDelete.extend(glob( 'tauscan_*' ))

            if CLEAN_POSTFITS:
                filesToDelete.extend(glob('higgsCombine_toyPostfit.MultiDimFit.mH125.root' ))


            if len(filesToDelete) == 0:
                print '{0}: {1} is already clean'.format(
                    'TESTMODE' if Tools.IsTestMode() else 'EXECUTING',
                    toyDir
                    )
                os.chdir( overDir )
                continue

            if Tools.IsTestMode():
                print 'TESTMODE: Would clean the following files from {0}:'.format( toyDir )
                print '    ' + '\n    '.join(filesToDelete)
            else:
                print 'EXECUTING: Cleaning the following files from {0}:'.format( toyDir )
                print '    ' + '\n    '.join(filesToDelete)
                for fileToDelete in filesToDelete:
                    if fileToDelete.startswith('tauscan') and CLEAN_TAUSCANS:
                        shutil.rmtree(fileToDelete)
                    else:
                        os.remove( fileToDelete )

            os.chdir( overDir )
        os.chdir(currentdir)


    #____________________________________________________________________
    if args.plotToys:

        currentdir = os.getcwd()
        overDir    = abspath(args.input)

        ONLY_LATEST_SCAN = True
        # ONLY_LATEST_SCAN = False

        os.chdir( overDir )
        toyDirs = glob( 'toy_*' )
        if len(toyDirs) == 0:
            Tools.ThrowError( 'Could not find any \'toy_*\' directories in {0}'.format(overDir) )

        for toyDir in toyDirs:
            os.chdir( toyDir )

            postfitWS   = abspath( 'higgsCombine_toyPostfit.MultiDimFit.mH125.root' )
            tauScandirs = glob( 'tauscan*' )

            if len(tauScandirs) == 0:
                print 'Skipping {0}; no tauscan available'.format(toyDir)
                os.chdir(overDir)
            elif ONLY_LATEST_SCAN and len(tauScandirs) > 1:
                # Interpret timestamps
                timestamps = [ re.search( r'\d\d\d\d\d\d_\d\d\d\d\d\d', basename(d) ).group() for d in tauScandirs ]
                # Sort tauscans simultaneously with timestamps and take first one
                tauScandirs = [ d for ts, d in sorted(zip( timestamps, tauScandirs ))]
                if args.verbose:
                    print 'Found more than 1 tauscan; take the latest one'
                    for d in tauScandirs:
                        if d == tauScandirs[-1]:
                            print ' > ' + d
                        else:
                            print '   ' + d
                tauScandirs = [ tauScandirs[-1] ]

            for tauScandir in tauScandirs:
                rootfiles = glob( '{0}/*.root'.format(tauScandir) )
                PlotTauScan( rootfiles, postfitWS, 'tauscan_{0}'.format( basename(toyDir).replace('/','') ) )

            os.chdir(overDir)
        os.chdir(currentdir)





########################################
# End of Main
########################################
if __name__ == "__main__":
    main()