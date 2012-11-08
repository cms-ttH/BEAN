#include "../interface/BTagReshaping.h"

BTagShape::BTagShape(){}
BTagShape::BTagShape(TFile *file ,const char * name,const std::vector<std::pair<float, float> > & cutsAndSF){
	TH1F * m_h = (TH1F *) file->Get(name);

	//compute equivalents
	std::vector<std::pair<float,float> > eq;
	int lastbin =2001;
	//float integral =  m_h->Integral(-1,lastbin);
	for(unsigned int i =0;i<cutsAndSF.size(); i++){
		float oldCut=cutsAndSF[i].first;
		float sf=cutsAndSF[i].second;
		float originalIntegral = m_h->Integral(m_h->FindBin(oldCut),lastbin);
		float originalLowEdge = m_h->GetBinLowEdge(m_h->FindBin(oldCut));
		//std::cout << std::endl<<    " Scale Factor : " << sf << std::endl;
		//      float target=originalIntegral/sf;
		float target=originalIntegral*sf;
		//std::cout << " Target " << target << " orig " << originalIntegral << std::endl;
		for(int j=lastbin; j> -1; j--)
		{
			if(m_h->Integral(j,lastbin)>= target)
			{
				//equivalents.push_back(std::pair<float,float>(originalLowEdge,h->GetBinLowEdge(j))); 
				eq.push_back(std::pair<float,float>(m_h->GetBinLowEdge(j),originalLowEdge)); 
				//std::cout << "Found at " << j << " was " << m_h->FindBin(oldCut) <<  std::endl;
				//std::cout << m_h->GetBinLowEdge(j) << " was " << originalLowEdge << " cut: " << oldCut <<  std::endl;
				//std::cout << m_h->GetBinLowEdge(j) << "\t" << originalLowEdge << std::endl;
				break;
			}
		}

	}


	//Interpolator

	std::vector<double> x; x.clear();
	std::vector<double> y; y.clear();
	x.push_back(0.);
	y.push_back(0.);
	for(unsigned int i  = 0 ; i < eq.size(); i++)
	{
		x.push_back(eq[eq.size()-i-1].first);
		y.push_back(eq[eq.size()-i-1].second);
	}
	//x.push_back(10.001);
	//y.push_back(10.001);
	x.push_back(1.001);
	y.push_back(1.001);

/*	cout << "Net result: =========== " << endl;
	for(unsigned int i=0; i<x.size(); i++){
		cout << "x = " << x.at(i) << "\t\ty = " << y.at(i) << endl;
	}
	cout << "======================= " << endl; //*/

	m_i =  new ROOT::Math::Interpolator(x,y,ROOT::Math::Interpolation::kLINEAR);
	//m_i =  new ROOT::Math::Interpolator(x,y,ROOT::Math::Interpolation::kCSPLINE);

	//for( int i=0; i<int(x.size()); i++ ) printf("x.push_back(%4.4f); y.push_back(%4.4f);\n", x[i], y[i] );

}

float BTagShape::eval(float x) { return m_i->Eval(x); } 


EtaPtBin::EtaPtBin(){}
EtaPtBin::EtaPtBin(float emin,float emax,float ptmin,float ptmax) : etaMin(emin), etaMax(emax), ptMin(ptmin), ptMax(ptmax) {}
bool EtaPtBin::contains(float eta,float pt) {return eta < etaMax && eta >= etaMin && pt < ptMax && pt >= ptMin ; } 
float EtaPtBin::centerEta() { return (etaMax+etaMin)/2.;}
float EtaPtBin::centerPt() { return (ptMax+ptMin)/2.;}

BinnedBTagShape::BinnedBTagShape(){}
BinnedBTagShape::BinnedBTagShape(std::vector<EtaPtBin> & bins, std::vector< std::vector<std::pair<float, float> > > &  cutsAndSF, TFile * f, const char * name):m_bins(bins)
{
	for(unsigned int i =0; i < bins.size(); i++)
	{
		m_shapes.push_back( BTagShape(f , name, cutsAndSF[i]));
	}
}

float  BinnedBTagShape::eval(float eta,float pt,float x) 
{
	for(unsigned int i =0; i < m_bins.size(); i++)
	{
		if(m_bins[i].contains(fabs(eta),pt)) return m_shapes[i].eval(x);
	}
	//    std::cout << "Cannot reshape eta pt discr "  << eta << " " << pt << " " << x << std::endl; 
	return x;
}


BTagShapeInterface::BTagShapeInterface(){}
BTagShapeInterface::BTagShapeInterface(std::string prefix, const char * file, float scaleBC, float scaleL) : m_file(new TFile(file))
{

	std::string h_hb = std::string( prefix + "_hb" );
	std::string h_hc = std::string( prefix + "_hc" );
	std::string h_hl = std::string( prefix + "_hl" );

	std::vector<EtaPtBin> binsBC;
	std::vector< std::vector<std::pair<float, float> > > cutsAndSFB;
	std::vector< std::vector<std::pair<float, float> > > cutsAndSFC;
	float charmFactor = 2. -1.0 ; //additional uncertainty for charm
	for(unsigned int i=0;i<myBeff.bins;i++)
	{
		EtaPtBin bin(-2.5,2.5,myBeff.ptmin[i],myBeff.ptmax[i]);
		binsBC.push_back(bin);
		std::vector<std::pair< float, float > > cutsAndSFbinB;
		std::vector<std::pair< float, float > > cutsAndSFbinC;
		float sft = myBeff.CSVT_SFb(bin.centerPt());
		sft+=scaleBC * myBeff.CSVT_SFb_error[i]; // add error
		cutsAndSFbinB.push_back(std::pair<float, float>(0.898, sft));
		sft+=scaleBC * myBeff.CSVT_SFb_error[i]*charmFactor; // charm additional error
		cutsAndSFbinC.push_back(std::pair<float, float>(0.898, sft));

		float sfm = myBeff.CSVM_SFb(bin.centerPt());
		sfm+=scaleBC * myBeff.CSVM_SFb_error[i]; // add error
		cutsAndSFbinB.push_back(std::pair<float, float>(0.679,sfm));
		sfm+=scaleBC * myBeff.CSVM_SFb_error[i]*charmFactor; // charm additional error
		cutsAndSFbinC.push_back(std::pair<float, float>(0.679,sfm));

		float sfl = myBeff.CSVL_SFb(bin.centerPt());
		sfl+=scaleBC * myBeff.CSVL_SFb_error[i]; // add error
		cutsAndSFbinB.push_back(std::pair<float, float>(0.244,sfl));
		sfl+=scaleBC * myBeff.CSVL_SFb_error[i]*charmFactor; // charm additional error
		cutsAndSFbinC.push_back(std::pair<float, float>(0.244,sfl));

		//std::cout << "SFs "  << i << " " << sfl << " " << sfm << " " << sft << std::endl;
		cutsAndSFB.push_back(cutsAndSFbinB);
		cutsAndSFC.push_back(cutsAndSFbinC);
	}  
	//underflow:
	{
		std::vector<std::pair< float, float > > cutsAndSFbinC;
		std::vector<std::pair< float, float > > cutsAndSFbinB;

		binsBC.push_back(EtaPtBin(-2.5,2.5,-9e99,myBeff.ptmin[0]));
		float sft = myBeff.CSVT_SFb(myBeff.ptmin[0]);
		sft+=scaleBC * 0.12;// add error
		cutsAndSFbinB.push_back(std::pair<float, float>(0.898, sft));
		sft+=scaleBC * 0.12*charmFactor; // charm additional error
		cutsAndSFbinC.push_back(std::pair<float, float>(0.898, sft));

		float sfm = myBeff.CSVM_SFb(myBeff.ptmin[0]);
		sfm+=scaleBC * 0.12; // add error
		cutsAndSFbinB.push_back(std::pair<float, float>(0.679,sfm));
		sfm+=scaleBC * 0.12*charmFactor; // charm additional error
		cutsAndSFbinC.push_back(std::pair<float, float>(0.679,sfm));

		float sfl = myBeff.CSVL_SFb(myBeff.ptmin[0]);
		sfl+=scaleBC * 0.12; // add error
		cutsAndSFbinB.push_back(std::pair<float, float>(0.244,sfl));
		sfl+=scaleBC * 0.12*charmFactor; // charm additional error
		cutsAndSFbinC.push_back(std::pair<float, float>(0.244,sfl));

		//std::cout << "Firstbin SFs " << sfl << " " << sfm << " " << sft << std::endl;
		cutsAndSFB.push_back(cutsAndSFbinB);
		cutsAndSFC.push_back(cutsAndSFbinC);

	}
	//overflow:
	{
		std::vector<std::pair< float, float > > cutsAndSFbinC;
		std::vector<std::pair< float, float > > cutsAndSFbinB;

		binsBC.push_back(EtaPtBin(-2.5,2.5,myBeff.ptmax[myBeff.bins-1],9e99));
		float sft = myBeff.CSVT_SFb(myBeff.ptmax[myBeff.bins-1]);
		sft+=scaleBC * myBeff.CSVT_SFb_error[myBeff.bins-1]*2;// add error
		cutsAndSFbinB.push_back(std::pair<float, float>(0.898, sft));
		sft+=scaleBC * myBeff.CSVT_SFb_error[myBeff.bins-1]*charmFactor; // charm additional error
		cutsAndSFbinC.push_back(std::pair<float, float>(0.898, sft));

		float sfm = myBeff.CSVM_SFb(myBeff.ptmax[myBeff.bins-1]);
		sfm+=scaleBC * myBeff.CSVM_SFb_error[myBeff.bins-1]*2; // add error
		cutsAndSFbinB.push_back(std::pair<float, float>(0.679,sfm));
		sfm+=scaleBC * myBeff.CSVM_SFb_error[myBeff.bins-1]*charmFactor; // charm additional error
		cutsAndSFbinC.push_back(std::pair<float, float>(0.679,sfm));

		float sfl = myBeff.CSVL_SFb(myBeff.ptmax[myBeff.bins-1]);
		sfl+=scaleBC * myBeff.CSVL_SFb_error[myBeff.bins-1]*2; // add error
		cutsAndSFbinB.push_back(std::pair<float, float>(0.244,sfl));
		sfl+=scaleBC * myBeff.CSVL_SFb_error[myBeff.bins-1]*charmFactor; // charm additional error
		cutsAndSFbinC.push_back(std::pair<float, float>(0.244,sfl));

		//std::cout << "Lastbin SFs " << sfl << " " << sfm << " " << sft << std::endl;
		cutsAndSFB.push_back(cutsAndSFbinB);
		cutsAndSFC.push_back(cutsAndSFbinC);
	}


	/*     m_b = new BinnedBTagShape(binsBC,cutsAndSFB,m_file,"hb"); */
	/*     m_c = new BinnedBTagShape(binsBC,cutsAndSFC,m_file,"hc"); */
	m_b = new BinnedBTagShape(binsBC,cutsAndSFB,m_file, h_hb.c_str() );
	m_c = new BinnedBTagShape(binsBC,cutsAndSFC,m_file, h_hc.c_str() );

	std::vector<EtaPtBin> binsL;
	std::vector< std::vector<std::pair<float, float> > > cutsAndSFL;
	//20-30 is also covered for mistag
	float ptmin[] = {20,30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
	float ptmax[] = {30,40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};
	float etamin[]  = {0,0.5,1.0,1.5};
	float etamax[]  = {0.5,1.0,1.5,2.5};
	size_t bins=15;

	for(unsigned int j = 0; j < 4 ; j++)
	{

		for(unsigned int i=0;i<bins;i++)
		{
			EtaPtBin bin(etamin[j],etamax[j],ptmin[i],ptmax[i]);
			binsL.push_back(bin);
			std::vector<std::pair< float, float > > cutsAndSFbinL;
			//float sft = mistag_CSVT(bin.centerEta(),bin.centerPt(),scaleL*1.5);
			float sft = mistag_CSVT(bin.centerEta(),bin.centerPt(),scaleL);
			cutsAndSFbinL.push_back(std::pair<float, float>(0.898, sft));

			float sfm = mistag_CSVM(bin.centerEta(),bin.centerPt(),scaleL);
			cutsAndSFbinL.push_back(std::pair<float, float>(0.679,sfm));

			float sfl = mistag_CSVL(bin.centerEta(),bin.centerPt(),scaleL);
			cutsAndSFbinL.push_back(std::pair<float, float>(0.244,sfl));

			//std::cout << "SFs light " << j << " " << i << " " << sfl << " " << sfm << " " << sft << std::endl;
			cutsAndSFL.push_back(cutsAndSFbinL);
		}

		//overflow:
		{
			std::vector<std::pair< float, float > > cutsAndSFbinL;

			binsL.push_back(EtaPtBin(etamin[j],etamax[j],ptmax[bins-1],9e99));
			float sft = mistag_CSVT((etamin[j]+etamax[j])/2.,ptmax[bins-1],scaleL*2);
			cutsAndSFbinL.push_back(std::pair<float, float>(0.898, sft));
			float sfm = mistag_CSVM((etamin[j]+etamax[j])/2.,ptmax[bins-1],scaleL*2);
			cutsAndSFbinL.push_back(std::pair<float, float>(0.679, sfm));
			float sfl = mistag_CSVL((etamin[j]+etamax[j])/2.,ptmax[bins-1],scaleL*2);
			cutsAndSFbinL.push_back(std::pair<float, float>(0.244, sfl));



			//std::cout << "SFs light " << sfl << " " << sfm << " " << sft << std::endl;
			cutsAndSFL.push_back(cutsAndSFbinL);
		}

	}

	/*     m_l = new BinnedBTagShape(binsL,cutsAndSFL,m_file,"hl"); */
	m_l = new BinnedBTagShape(binsL,cutsAndSFL,m_file, h_hl.c_str() );


}

float BTagShapeInterface::reshape(float eta, float pt, float csv, int flav){
	if(csv < 0) return csv;
	if(csv >= 1.) return csv; 
	if(flav == 0) return csv;  
	if(fabs(flav) == 5) return  m_b->eval(eta,pt,csv); 
	if(fabs(flav) == 4) return  m_c->eval(eta,pt,csv);
	if(fabs(flav) != 4 && fabs(flav) != 5) return m_l->eval(eta,pt,csv);

	return -10000; 

}


