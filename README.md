# Boson Exploration Analysis Ntuple

More information and instructions are on the [Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHiggs).

## Installation
Follow These Steps:

    cmsrel CMSSW_5_3_11
    cd CMSSW_5_3_11/src/
    cmsenv

    # Basic third-party packages
    git cms-merge-topic muell149:merged_5_3_11
    git clone https://github.com/latinos/UserCode-CMG-CMGTools-External CMGTools/External  ..notquite correct
    git clone https://github.com/latinos/UserCode-EGamma-EGammaAnalysisTools.git EGamma/EGammaAnalysisTools
    pushd EGamma/EGammaAnalysisTools/
    git checkout tags/V00-00-08
    wget -r http://nd.edu/~abrinke1/ElectronEffectiveArea.h -O interface/ElectronEffectiveArea.h
    cd data
    cat download.url | xargs wget
    popd
    git cms-addpkg ElectroWeakAnalysis/Utilities
    


    # Our software
    git clone https://github.com/cms-ttH/BEAN.git -b 

    # Build
    scram b -j32

This might require an environment set up as described [here](http://wiki.crc.nd.edu/wiki/index.php/NDCMS_SettingUpEnvironment).
