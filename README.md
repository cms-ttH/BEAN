# Boson Exploration Analysis Ntuple

More information and instructions are on the [Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHiggs).

## Installation
Follow These Steps:

    cmsrel CMSSW_5_3_11
    cd CMSSW_5_3_11/src/
    cmsenv
    git cms-merge-topic --unsafe muell149:merged_5_3_11

The '--unsafe' option prevents the merge from cloning dependent packages that we actually don't use

    git cms-addpkg AnalysisDataFormats/TopObjects 
    git clone https://github.com/muell149/UserCode-CMG-CMGTools-External.git CMGTools/External
    git cms-addpkg CondFormats/EgammaObjects
    git clone https://github.com/muell149/UserCode-EGamma-EGammaAnalysisTools.git EGamma/EGammaAnalysisTools
    git cms-addpkg ElectroWeakAnalysis/Utilities
    git clone https://github.com:latinos/UserCode-sixie-Muon-MuonAnalysisTools.git Muon/MuonAnalysisTools/
    git cms-addpkg RecoBTag/SecondaryVertex
    git cms-addpkg RecoVertex/AdaptiveVertexFilter
    git clone https://github.com/muell149/TopAnalysis.git

Our software
    
    git clone https://github.com/cms-ttH/BEAN.git -b charlieProjection

Build

    scram b -j32
