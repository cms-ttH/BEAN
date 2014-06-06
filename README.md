# Boson Exploration Analysis Ntuple

More information and instructions are on the [Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHiggs).

## Installation
Follow These Steps:

    cmsrel CMSSW_5_3_15_patch1
    cd CMSSW_5_3_15_patch1/src/
    cmsenv
    scram setup lhapdffull
    git cms-merge-topic --unsafe muell149:clean_merged_5_3_15p1_w_tau2014

The '--unsafe' option prevents the merge from cloning dependent packages that we actually don't use

    git cms-addpkg AnalysisDataFormats/TopObjects 
    git clone https://github.com/muell149/UserCode-CMG-CMGTools-External.git CMGTools/External
    git cms-addpkg CondFormats/EgammaObjects
    git clone https://github.com/muell149/UserCode-EGamma-EGammaAnalysisTools.git EGamma/EGammaAnalysisTools
    git clone git@github.com:latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools/
    git cms-addpkg RecoVertex/AdaptiveVertexFinder
    git clone https://github.com/muell149/TopAnalysis.git

Our software
    
    git clone https://github.com/cms-ttH/BEAN.git -b chPrj_5_3_15p1

Build

    scram b -j 64
