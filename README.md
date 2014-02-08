# Boson Exploration Analysis Ntuple

More information and instructions are on the [Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHiggs).

## Installation

To use this package, perform the following steps (using release 5.3.8, patch 1):

    cmsrel CMSSW_5_3_8_patch1
    cd CMSSW_5_3_8_patch1/src/
    cmsenv

    # Basic third-party packages
    cvs co -d PhysicsTools/NtupleUtils UserCode/Bicocca/PhysicsTools/NtupleUtils
    cvs co -d SusyAnalysis -r V01-02-02 SusyAnalysis/EventSelector/interface/uncorrectionTypeMET.h
    cvs co -r V07-00-01    TopQuarkAnalysis/Configuration
    cvs co -r V06-07-11-01 TopQuarkAnalysis/TopTools
    cvs co -r V06-05-06-07 DataFormats/PatCandidates
    cvs co -r V00-02-14    DataFormats/StdDictionaries
    cvs co -r V00-00-70    FWCore/GuiBrowsers
    cvs co -r V08-09-52    PhysicsTools/PatAlgos
    cvs co -r V03-09-23    PhysicsTools/PatUtils
    cvs co -r V00-03-15    CommonTools/ParticleFlow
    cvs co -r V00-00-12    CommonTools/RecoUtils
    cvs co -r V04-06-09    JetMETCorrections/Type1MET
    cvs co -r V01-08-00    RecoBTag/SecondaryVertex
    cvs co -r V00-00-08    RecoMET/METAnalyzers
    cvs co -r V00-00-07    RecoMET/METFilters
    cvs co -r V15-01-11    RecoParticleFlow/PFProducer
    cvs co -r V02-02-00    RecoVertex/AdaptiveVertexFinder

    mkdir Muon
    git clone https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools

    cvs co -r V00-00-13 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
    cd EGamma/EGammaAnalysisTools/data
    cat download.url | xargs wget
    cd -
    cd EGamma/EGammaAnalysisTools/interface/
    wget -r http://nd.edu/~abrinke1/ElectronEffectiveArea.h -O ElectronEffectiveArea.h
    cd -
    cvs co -r V01-04-23 RecoTauTag/RecoTau
    cvs co -r V01-04-10 RecoTauTag/Configuration
    cvs co -r V00-04-00 CondFormats/EgammaObjects

    mkdir CMGTools
    git clone https://github.com/latinos/UserCode-CMG-CMGTools-External CMGTools/External

    # Our software
    git clone https://github.com/cms-ttH/BEAN.git

    # Build
    scram b -j32

This might require an environment set up as described [here](http://wiki.crc.nd.edu/wiki/index.php/NDCMS_SettingUpEnvironment).

As CMS has moved from CVS to git, the recipe above may not work due to problems accessing CVS. The git-only recipe below is under
construction (currently broken for BEAN making.)

	cmsrel CMSSW_5_3_14
	cd CMSSW_5_3_14/src
	cmsenv
	# from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMETRecipe53X
	git cms-addpkg PhysicsTools/PatAlgos
	git cms-merge-topic cms-analysis-tools:5_3_14-updateSelectorUtils
	git cms-merge-topic -u TaiSakuma:53X-met-131120-01
	git clone https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools/
	git clone https://github.com/latinos/UserCode-CMG-CMGTools-External CMGTools/External
	git clone https://github.com/h2gglobe/EGammaAnalysisTools.git EGamma/EGammaAnalysisTools
	git cms-addpkg RecoJets/JetProducers
	mkdir -p RecoJets/JetProducers/data
	wget --no-check-certificate -r https://raw2.github.com/h2gglobe/External/master/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml -O RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml

	# common ttH software
	git clone https://github.com/cms-ttH/BEAN.git

	# Build
	scram b -j9
