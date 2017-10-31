#!/usr/bin/env python
"""
Thomas Klijnsma
"""

########################################
# Imports
########################################

import os, argparse, re, sys
from os.path import *
from glob import glob
from array import array

from time import strftime
datestr = strftime( '%b%d' )

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
    
    # ======================================
    # To add the regularization terms to the .txt datacards

    parser.add_argument( '-v', '--verbose', action='store_true' )
    parser.add_argument( '-t', '--test', action='store_true' )
    parser.add_argument( '-i', '--input', type=str )

    parser.add_argument( '--processTextDatacards', action='store_true')
    parser.add_argument( '--t2ws', action='store_true')
    parser.add_argument( '--makePostfit', action='store_true')
    parser.add_argument( '--tauscan', action='store_true')

    parser.add_argument( '--PlotTauScan', action='store_true')

    parser.add_argument( '--bestfitWS', type=str )

    # parser.add_argument( '--string', type=str, default='default', help='default string' )
    # parser.add_argument( '--list', metavar='N', type=str, nargs='+', help='list of strings' )

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
            # 
            '--saveWorkspace',
            '-m 125',
            '--setPhysicsModelParameters tau=0.0',
            '--freezeNuisances tau',
            '-n _POSTFIT_{0}'.format( wsTag ),
            # 
            '; mv {0} {1}'.format( postfitFilename, ws.replace( '.root', '_POSTFIT.root' ) )
            ]

        Tools.executeCommand( cmd )



    #____________________________________________________________________
    if args.tauscan:

        ws              = abspath(args.input)
        wsTag           = basename(ws).replace('/','').replace('.root','')
        # postfitFilename = join( os.getcwd(), 'higgsCombine_POSTFIT_{0}.MultiDimFit.mH125.root'.format( wsTag ) )


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


        # ======================================
        # Parse command

        POIs = Tools.ListSet( ws, 'POI' )
        POIs.append( 'tau' )
        POIrangeStr = ':'.join([ '{0}=-1.0,4.0'.format(POI) for POI in POIs if re.match( r'r\d+', POI ) ])

        extraOptions = [
            '-P tau',
            '--setPhysicsModelParameterRanges tau={0},{1}'.format( *tauRange ) + ':' + POIrangeStr,
            # '--saveSpecifiedNuis tau',
            # 
            '--floatOtherPOIs=1',
            '-m 125.00',
            '--snapshotName MultiDimFit',
            # '--redefineSignalPOIs {0}'.format( ','.join(POIs) ),
            ]

        if FREEZE_PDF_INDICES:
            pdfIndicesToFreeze = Tools.ListOfPDFIndicesToFreeze( ws, verbose=args.verbose )
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
            notOnBatch    = True,
            jobDirectory  = 'Scan_{0}_tauscan'.format(datestr),
            fastscan      = False,
            asimov        = False,
            jobPriority   = 0, # Highest; set to negative for lower
            extraOptions  = extraOptions,
            )



    if args.PlotTauScan:

        scandir = args.input

        rootfiles = glob( '{0}/*.root'.format(scandir) )

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
        if args.bestfitWS:
            bestfitPoint = Tools.ConvertTChainToArray(
                [ args.bestfitWS ],
                treeName        = 'limit',
                variablePattern = '*',
                returnStyle     = 'dict',
                verbose         = args.verbose
                )
            bestfitFound = True



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

        print axisLineAlpha

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
                H.SetBinContent( iMu+1, bestfitPoint[mu][0] )

            H.SetLineWidth(2)
            H.SetLineColor(4)
            H.Draw('HISTSAME')



        SaveC( 'tauscan' )






########################################
# End of Main
########################################
if __name__ == "__main__":
    main()