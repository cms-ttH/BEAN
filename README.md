# Boson Exploration Analysis Ntuple

More information and instructions are on the [Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHiggs).

## Installation

Warning: Make sure to have added the appropriate public ssh key to your GitHub account before proceeding!

To use this package, perform the following steps (using release 5.3.8, patch 1):

	cmsrel CMSSW_5_3_8_patch1
	cd CMSSW_5_3_8_patch1/src
	cmsenv

	git cms-cvs-history import V07-00-01 TopQuarkAnalysis/Configuration
	git cms-cvs-history import V06-07-11-01 TopQuarkAnalysis/TopTools
	git cms-cvs-history import  V06-05-06-07 DataFormats/PatCandidates
	git cms-cvs-history import V00-02-14 DataFormats/StdDictionaries
	git cms-cvs-history import V00-00-70 FWCore/GuiBrowsers
	git cms-cvs-history import V08-09-52 PhysicsTools/PatAlgos
	git cms-cvs-history import V03-09-23 PhysicsTools/PatUtils
	git cms-cvs-history import V00-03-15 CommonTools/ParticleFlow
	git cms-cvs-history import V00-00-12 CommonTools/RecoUtils
	git cms-cvs-history import V04-06-09 JetMETCorrections/Type1MET
	git cms-cvs-history import V01-08-00 RecoBTag/SecondaryVertex
	git clone https://github.com/cms-analysis/RecoMET-METAnalyzers.git RecoMET/METAnalyzers
	cd RecoMET/METAnalyzers
	git checkout tags/RecoMET-METAnalyzers-V00-00-08
	cd -
	git cms-cvs-history import V00-00-07 RecoMET/METFilters
	git cms-cvs-history import V15-01-11 RecoParticleFlow/PFProducer
	git cms-cvs-history import V02-02-00 RecoVertex/AdaptiveVertexFinder
	git clone git@github.com:latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools/
	git clone https://github.com/latinos/UserCode-EGamma-EGammaAnalysisTools.git EGamma/EGammaAnalysisTools
	pushd EGamma/EGammaAnalysisTools/
	git checkout tags/V00-00-08
	wget -r http://nd.edu/~abrinke1/ElectronEffectiveArea.h -O interface/ElectronEffectiveArea.h
	cd data
	cat download.url | xargs wget
	popd

	git cms-cvs-history import V01-04-23 RecoTauTag/RecoTau
	git cms-cvs-history import V01-04-10 RecoTauTag/Configuration
	git cms-cvs-history import V00-04-00 CondFormats/EgammaObjects
	git clone https://github.com/latinos/UserCode-CMG-CMGTools-External CMGTools/External
	cd CMGTools/External/
	git checkout tags/V00-02-10
	cd -

	git clone https://github.com/cms-ttH/BEAN.git
	scram b -j24

This might require an environment set up as described [here](http://wiki.crc.nd.edu/wiki/index.php/NDCMS_SettingUpEnvironment).
As CMS has stopped supporting CVS, the recipe above has been adapted for git, but not yet validated against the most recent set of
BEANs, which were produced using [this recipe](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHiggs#53X_prescription).
