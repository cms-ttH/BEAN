# Boson Exploration Analysis Ntuple

More information and instructions are on the [Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHiggs).

## Installation
Follow These Steps:

    cmsrel CMSSW_5_3_11
    cd CMSSW_5_3_11/src/
    cmsenv

    git cms-merge-topic --unsafe muell149:merged_5_3_11
    #the --unsafe option prevents the merge from cloning dependent packages that we actually don't use
    git cms-addpkg AnalysisDataFormats/TopObjects 
    git clone git@github.com:cms-ttH/BEAN.git -b merged_5_3_11
    #CMGTools
    git cms-addpkg CondFormats/EgammaObjects
    #Egamma 
    git cms-addpkg ElectroWeakAnalysis/Utilities
    #Muon
    git cms-addpkg RecoBTag/SecondaryVertex
    git cms-addpkg RecoVertex/AdaptiveVertexFilter* try it
    #SusyAnalysis
    #TopAnalysis
    
    

    # Our software
    git clone https://github.com/cms-ttH/BEAN.git -b merged_5_3_11

    # Build
    scram b -j32

This might require an environment set up as described [here](http://wiki.crc.nd.edu/wiki/index.php/NDCMS_SettingUpEnvironment).
