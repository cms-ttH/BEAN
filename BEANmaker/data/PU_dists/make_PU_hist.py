import os
import ROOT
ROOT.gROOT.SetBatch(True)

def main():

    hists = {}
    
    era_lumis = {}
    era_lumis['2012AB'] = 5.29
    era_lumis['2012C'] = 6.90
    era_lumis['2012D'] = 7.27
    
    root_file = ROOT.TFile('pu_distributions_SimGeneral.root', 'RECREATE')
    hists['2012'] = ROOT.TH1D('puDist_2012', 'puDist_2012', 60, 0, 60)
    
    for era in era_lumis:
        print '---------- %s ---------' % era
        #PU_file = os.popen('mix_%s_Profile_PoissonOOTPU_cfi.py' % era[0])
        hists[era] = ROOT.TH1D('puDist_%s' % era, 'puDist_%s' % era, 60, 0, 60)

        string1 = '/SDT/lxr/source/SimGeneral/MixingModule/python/mix_%s_Profile_PoissonOOTPU_cfi.py?v=CMSSW_5_3_21' % era
        string2 = '                 '
        
        all_lines = open('mix_%s_Profile_PoissonOOTPU_cfi.py' % era).read().splitlines()
        lines = []
        for line in all_lines:
            if string1 in line and string2 in line and ',' in line:
                lines.append(line)
        
        for bin in range(1,61):
            try:
                value = lines[bin].split(string2)[1]
                if ',' in value:
                    value = value.replace(',','')
            except:
                print 'In bin %d no line - filling with 0' % bin
                hists[era].SetBinContent(bin+1, 0.0)
                continue

            try:
                #print 'Filling bin %d with value %0.15f' % (bin, float(value))
                hists[era].SetBinContent(bin+1, float(value))
            except:
                print 'In bin %d value is %s - filling with 0' % (bin, value)
                hists[era].SetBinContent(bin+1, 0.0)

        #print 'Hist for era %s has integral %0.3f' % (era, hists[era].Integral())
        hists[era].Scale(1.0/hists[era].Integral())

        hists[era].Scale(era_lumis[era]/19.46)
        hists['2012'].Add( hists['2012'], hists[era] )

        hists[era].Scale(19.46/era_lumis[era])
        hists[era].SetDirectory(root_file)
        hists[era].Write()

    hists['2012'].SetDirectory(root_file)
    hists['2012'].Write()                                        
            

if __name__ == '__main__':
    main()
        
