# Boson Exploration Analysis Ntuple

More information and instructions are on the [Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHiggs).

## Installation

To use this package, perform the following steps (using release 5.3.18 on SLC6 machines only):

	cmsrel CMSSW_5_3_18
	cd CMSSW_5_3_18/src
	cmsenv

	git cms-merge-topic --unsafe muell149:5_3_18_SL6_upgrade

	git cms-cvs-history import V07-00-01 TopQuarkAnalysis/Configuration
	git cms-cvs-history import V06-07-11-01 TopQuarkAnalysis/TopTools
	git cms-cvs-history import  V06-05-06-07 DataFormats/PatCandidates
	git cms-cvs-history import V00-02-14 DataFormats/StdDictionaries
	git cms-cvs-history import V00-00-70 FWCore/GuiBrowsers
	git cms-cvs-history import V04-06-09 JetMETCorrections/Type1MET
	git cms-cvs-history import V15-01-11 RecoParticleFlow/PFProducer
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
	
	git clone git@github.com:cms-ttH/BEAN.git -b SL6_upgrade
	scram b -j32


This might require an environment set up as described [here](http://wiki.crc.nd.edu/wiki/index.php/NDCMS_SettingUpEnvironment).
As CMS has stopped supporting CVS, the recipe above has been adapted for git, but not yet validated against the most recent set of
BEANs, which were produced using [this recipe](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTbarHiggs#53X_prescription).
