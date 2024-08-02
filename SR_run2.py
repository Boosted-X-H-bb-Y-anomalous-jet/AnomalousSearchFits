from time import time
from TwoDAlphabet import plot
from TwoDAlphabet.twoDalphabet import MakeCard, TwoDAlphabet
from TwoDAlphabet.alphawrap import BinnedDistribution, ParametricFunction
from TwoDAlphabet.helpers import make_env_tarball, cd, execute_cmd
from TwoDAlphabet.ftest import FstatCalc
import os
import numpy as np

def _get_other_region_names(pass_reg_name):
    return pass_reg_name, pass_reg_name.replace('Fail','Pass')

def _select_signal(row, args):
    signame = args[0]
    poly_order = args[1]
    if row.process_type == 'SIGNAL':
        if signame in row.process:
            return True
        else:
            return False
    elif 'Background' in row.process:
        if row.process == 'Background_'+poly_order:
            return True
        elif row.process == 'Background':
            return True
        else:
            return False
    else:
        return True
#------------------------------------------------------------------------------------------------------------
# NEED TO ADJUST THIS FUNCTION
#------------------------------------------------------------------------------------------------------------
def _load_CR_rpf(poly_order):
    twoD_CRonly = TwoDAlphabet('XHY_CRfl','test4failtoloose.json', loadPrevious=True)
    params_to_set = twoD_CRonly.GetParamsOnMatch('rpf.*'+poly_order, 'MX_2000_MY_800_area', 'b')
    return {k:v['val'] for k,v in params_to_set.items()}
#------------------------------------------------------------------------------------------------------------

def _load_CR_rpf_as_SR(poly_order):
    params_to_set = {}
    for k,v in _load_CR_rpf(poly_order).items():
        params_to_set[k.replace('CR','SR')] = v
    return params_to_set

def _generate_constraints(nparams):
    out = {}
    for i in range(nparams):
        if i == 0:
            out[i] = {"MIN":-50,"MAX":50,"NOM":0.005}
        else:
            out[i] = {"MIN":-50,"MAX":50,"NOM":0.0}
    return out

_rpf_options = {
    '0': {
        'form': '@0',
        'constraints': _generate_constraints(1)
    },
    "1":{
        'form': '@0*(1+@1*x+@2*y)',
        'constraints': _generate_constraints(3)        
    },
    "2":{
        'form': '@0*(1+@1*x+@2*y+@3*x*x+@4*y*y+@5*x*y)',
        'constraints': _generate_constraints(6)        
    },
    "3":{
        'form': '@0*(1+@1*x+@2*y+@3*x*x+@4*y*y+@5*x*y+@6*x*x*y+@7*x*y*y+@8*x*x*x+@9*y*y*y)',
        'constraints': _generate_constraints(10)        
    }
}

def test_make(working_area='CR_run2',json_name='CR_run2.json'):
    twoD = TwoDAlphabet(working_area,json_name,loadPrevious=False)
    qcd_hists = twoD.InitQCDHists()

    for f,p in [_get_other_region_names(r) for r in twoD.ledger.GetRegions() if 'Fail' in r]:
        binning_f, _ = twoD.GetBinningFor(f)
        fail_name = 'Background_'+f
        qcd_f = BinnedDistribution(
            fail_name, qcd_hists[f],
            binning_f, constant=False
        )
        twoD.AddAlphaObj('Background',f,qcd_f,title='Multijet')

        for opt_name, opt in _rpf_options.items():

            qcd_rpf = ParametricFunction(
                        fail_name.replace('Fail','rpf')+'_'+opt_name,
                        binning_f, opt['form'],
                        constraints=opt['constraints']
                )
            qcd_p = qcd_f.Multiply(fail_name.replace('Fail','Pass')+'_'+opt_name,qcd_rpf)
            twoD.AddAlphaObj('Background'+'_'+opt_name,p,qcd_p,title='Multijet')

    twoD.Save()

def add_vae_sf_to_card(working_area,subtag,signal,SF_file):
    
    def append_to_datacard(era, signal, value, uncertainty_low, uncertainty_high, datacard_path):
        with open(datacard_path, 'r') as datacard:
            lines = datacard.readlines()
        
        for i, line in enumerate(lines):
            if line.startswith('kmax'):
                kmax_parts = line.split()
                if kmax_parts[1] != '*':
                    kmax_value = int(kmax_parts[1])
                    kmax_parts[1] = str(kmax_value + 1)
                    lines[i] = ' '.join(kmax_parts) + '\n'
                break

        vae_sf_line = f"vae_sf_{era} rateParam SR_Pass* {signal}_{era} {value}\n"
        rate_param_line = f"vae_sf_{era} param {value} -{uncertainty_low}/+{uncertainty_high}\n"
        
        lines.append(vae_sf_line)
        lines.append(rate_param_line)
        
        with open(datacard_path, 'w') as datacard:
            datacard.writelines(lines)


    card = f"{working_area}/{subtag}/card.txt"

    with open(SF_file, 'r') as file:
        for line in file:
            parts = line.strip().split(',')
            era = parts[0]
            signal_sf = parts[1]
            if signal!=signal_sf:
                continue
            value = parts[2]
            uncertainty_high = parts[3]
            uncertainty_low = parts[4]
            
            # Append the information to the datacard
            append_to_datacard(era, signal, value, uncertainty_low, uncertainty_high, card)

def test_fit(signal, tf='', working_area='CR_run2', defMinStrat=0,  extra='--robustFit 1'): #extra='--robustHesse 1'
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)

    subset = twoD.ledger.select(_select_signal, signal, tf)
    twoD.MakeCard(subset, '{}-{}_area'.format(signal, tf))
    add_vae_sf_to_card(working_area,'{}-{}_area'.format(signal, tf),signal,"SFs_vae.txt")

    twoD.MLfit('{}-{}_area'.format(signal,tf),rMin=-20,rMax=20,verbosity=1,extra=extra)


def test_plot(signal, tf='', working_area='CR_run2'):
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    subset = twoD.ledger.select(_select_signal, signal, tf)
    twoD.StdPlots('{}-{}_area'.format(signal,tf), subset)

def _gof_for_FTest(twoD, subtag, card_or_w='card.txt'):

    run_dir = twoD.tag+'/'+subtag
    
    with cd(run_dir):
        gof_data_cmd = [
            'combine -M GoodnessOfFit',
            '-d '+card_or_w,
            '--algo=saturated',
            '-n _gof_data'
        ]

        gof_data_cmd = ' '.join(gof_data_cmd)
        execute_cmd(gof_data_cmd)

def test_FTest(poly1, poly2, working_area="CR_run2",signal='MX1400_MY90'):
    '''
    Perform an F-test using existing working areas
    '''
    
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    binning = twoD.binnings['default']
    nBins = (len(binning.xbinList)-1)*(len(binning.ybinList)-1)
    
    # Get number of RPF params and run GoF for poly1
    params1 = twoD.ledger.select(_select_signal, signal, poly1).alphaParams
    rpfSet1 = params1[params1["name"].str.contains("rpf")]
    nRpfs1  = len(rpfSet1.index)
    _gof_for_FTest(twoD, f'{signal}-{poly1}_area', card_or_w='card.txt')
    gofFile1 = f'{working_area}/{signal}-{poly1}_area/higgsCombine_gof_data.GoodnessOfFit.mH120.root'

    # Get number of RPF params and run GoF for poly2
    params2 = twoD.ledger.select(_select_signal, signal, poly2).alphaParams
    rpfSet2 = params2[params2["name"].str.contains("rpf")]
    nRpfs2  = len(rpfSet2.index)
    _gof_for_FTest(twoD, f'{signal}-{poly2}_area', card_or_w='card.txt')
    gofFile2 = f'{working_area}/{signal}-{poly2}_area/higgsCombine_gof_data.GoodnessOfFit.mH120.root'

    base_fstat = FstatCalc(gofFile1,gofFile2,nRpfs1,nRpfs2,nBins)
    print(base_fstat)

    def plot_FTest(base_fstat,nRpfs1,nRpfs2,nBins):
        from ROOT import TF1, TH1F, TLegend, TPaveText, TLatex, TArrow, TCanvas, kBlue, gStyle
        gStyle.SetOptStat(0000)

        if len(base_fstat) == 0: base_fstat = [0.0]

        ftest_p1    = min(nRpfs1,nRpfs2)
        ftest_p2    = max(nRpfs1,nRpfs2)
        ftest_nbins = nBins
        fdist       = TF1("fDist", "[0]*TMath::FDist(x, [1], [2])", 0,max(10,1.3*base_fstat[0]))
        fdist.SetParameter(0,1)
        fdist.SetParameter(1,ftest_p2-ftest_p1)
        fdist.SetParameter(2,ftest_nbins-ftest_p2)

        pval = fdist.Integral(0.0,base_fstat[0])
        print('P-value: %s'%pval)

        c = TCanvas('c','c',800,600)    
        c.SetLeftMargin(0.12) 
        c.SetBottomMargin(0.12)
        c.SetRightMargin(0.1)
        c.SetTopMargin(0.1)
        ftestHist_nbins = 30
        ftestHist = TH1F("Fhist","",ftestHist_nbins,0,max(10,1.3*base_fstat[0]))
        ftestHist.GetXaxis().SetTitle("F = #frac{-2log(#lambda_{1}/#lambda_{2})/(p_{2}-p_{1})}{-2log#lambda_{2}/(n-p_{2})}")
        ftestHist.GetXaxis().SetTitleSize(0.025)
        ftestHist.GetXaxis().SetTitleOffset(2)
        ftestHist.GetYaxis().SetTitleOffset(0.85)
        
        ftestHist.Draw("pez")
        ftestobs  = TArrow(base_fstat[0],0.25,base_fstat[0],0)
        ftestobs.SetLineColor(kBlue+1)
        ftestobs.SetLineWidth(2)
        fdist.Draw('same')

        ftestobs.Draw()
        tLeg = TLegend(0.6,0.73,0.89,0.89)
        tLeg.SetLineWidth(0)
        tLeg.SetFillStyle(0)
        tLeg.SetTextFont(42)
        tLeg.SetTextSize(0.03)
        tLeg.AddEntry(ftestobs,"observed = %.3f"%base_fstat[0],"l")
        tLeg.AddEntry(fdist,"F-dist, ndf = (%.0f, %.0f) "%(fdist.GetParameter(1),fdist.GetParameter(2)),"l")
        tLeg.Draw("same")

        model_info = TPaveText(0.2,0.6,0.4,0.8,"brNDC")
        model_info.AddText('p1 = '+poly1)
        model_info.AddText('p2 = '+poly2)
        model_info.AddText("p-value = %.2f"%(1-pval))
        model_info.Draw('same')
        
        latex = TLatex()
        latex.SetTextAlign(11)
        latex.SetTextSize(0.06)
        latex.SetTextFont(62)
        latex.SetNDC()
        latex.DrawLatex(0.12,0.91,"CMS")
        latex.SetTextSize(0.05)
        latex.SetTextFont(52)
        latex.DrawLatex(0.65,0.91,"Preliminary")
        latex.SetTextFont(42)
        latex.SetTextFont(52)
        latex.SetTextSize(0.045)
        c.SaveAs(working_area+'/ftest_{0}_vs_{1}_notoys.png'.format(poly1,poly2))

    plot_FTest(base_fstat,nRpfs1,nRpfs2,nBins)

def test_GoF(signal, working_area,tf='', condor=True,extra=''):
    #assert SRorCR == 'CR'
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    signame = signal
    if not os.path.exists(twoD.tag+'/'+signame+'-{}_area/card.txt'.format(tf)):
        print('{}/{}-area/card.txt does not exist, making card'.format(twoD.tag,signame))
        subset = twoD.ledger.select(_select_signal, signame, tf)
        twoD.MakeCard(subset, signame+'_area')
        add_vae_sf_to_card(working_area,'{}-{}_area'.format(signal, tf),signal,"SFs_vae.txt")
    if condor == False:
        twoD.GoodnessOfFit(
            signame+'-{}_area'.format(tf), ntoys=100, freezeSignal=0,
            condor=False,extra=extra
        )
    else:
        twoD.GoodnessOfFit(
            signame+'-{}_area'.format(tf), ntoys=500, freezeSignal=0,
            condor=True, njobs=10,extra=extra
        )

def test_GoF_plot(signal, working_area, tf='', condor=True):
    '''Plot the GoF in ttbar<SRorCR>/TprimeB-<signal>_area (condor=True indicates that condor jobs need to be unpacked)'''
    signame = signal
    plot.plot_gof(working_area,'{}-{}_area'.format(signame,tf), condor=condor)

def load_RPF(twoD):
    params_to_set = twoD.GetParamsOnMatch('rpf.*', 'TprimeB-1800-125-_area', 'b')
    return {k:v['val'] for k,v in params_to_set.items()}

def test_SigInj(working_area, signal, rpfOrder, rpfParams, r, condor=False,njobs=10,scale_rpf=1.0,extra=''):
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    if scale_rpf!=1.0:
        print(rpfParams)
        print("Scaling rpf params for sig injection")
        for rpf_name, rpf_val in rpfParams.items():
            if("par0" in rpf_name):
            #Sometimes we want to inflate bkg to check if positive signal bias comes from low stats
            #par0 sets the overall scale
                rpfParams[rpf_name] = rpf_val*scale_rpf
        print(rpfParams)

    twoD.SignalInjection(
        '{}-{}_area'.format(signal, rpfOrder),
        injectAmount = r,       # amount of signal to inject (r=0 <- bias test)
        ntoys = 500,
        blindData = False,      # working with toy data, no need to blind
        setParams = rpfParams,     # give the toys the same RPF params
        verbosity = 0,
        extra = extra,
        condor = condor,
        njobs=njobs)

def test_SigInj_plot(working_area,signal,rpfOrder, r, condor=False):
    signal_area='{}-{}_area'.format(signal, rpfOrder)
    plot.plot_signalInjection(working_area,signal_area, injectedAmount=r, condor=condor)

def test_Impacts(working_area,signal,rpfOrder,extra='-t -1'):
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    signal_area='{}-{}_area'.format(signal, rpfOrder)
    twoD.Impacts(signal_area, cardOrW='card.txt', extra=extra)

def _load_fit_rpf(working_area,signal_name,rpfOrder,json_file):
    twoD_blindFit = TwoDAlphabet(working_area,json_file, loadPrevious=True)
    params_to_set = twoD_blindFit.GetParamsOnMatch('rpf.*', f'{signal_name}-{rpfOrder}_area', 's')
    return {k:v['val'] for k,v in params_to_set.items()}


def test_limit(working_area,signal,rpfOrder,json_file,blind=True,rpf_params={}):
    '''Perform a blinded limit. To be blinded, the Combine algorithm (via option `--run blind`)
    will create an Asimov toy dataset from the pre-fit model. Since the TF parameters are meaningless
    in our true "pre-fit", we need to load in the parameter values from a different fit so we have
    something reasonable to create the Asimov toy. 
    '''
    twoD = TwoDAlphabet(working_area, json_file, loadPrevious=True)

    #Make a subset and card as in test_fit()
    subset = twoD.ledger.select(_select_signal, signal, rpfOrder)
    twoD.MakeCard(subset, '{}-{}_area'.format(signal, rpfOrder))
    add_vae_sf_to_card(working_area,'{}-{}_area'.format(signal, rpfOrder),signal,"SFs_vae.txt")
    #Run the blinded limit with our dictionary of TF parameters
    twoD.Limit(
        subtag='{}-{}_area'.format(signal, rpfOrder),
        blindData=blind,
        verbosity=0,
        setParams=rpf_params,
        condor=False
    )

if __name__ == "__main__":
    #make_env_tarball()
    working_area_SR='SR_run2'
    test_make(working_area=working_area_SR,json_name='SR_run2.json')
    rpf_order="0"
    load_rpf_from_signal_name = "MX1400_MY90"
    params_to_set = _load_fit_rpf("CR_run2","MX1400_MY90",rpf_order,"CR_run2.json")
    params_to_set = {key.replace('CR', 'SR'): value for key, value in params_to_set.items()}#Replace CR parameter names with SR
    print(params_to_set)
    # #Btw. signal is normalized to xsec=5fb!
    # #Signal injection part
    for MX in [1400,2000]:
        signal=f"MX{MX}_MY90"
        #test_fit(signal=signal,working_area=working_area_SR,tf=rpf_order)#Uncomment if needed to create a card
        test_SigInj(working_area_SR, signal, rpf_order, params_to_set, r=0.0, condor=True,scale_rpf=1.0)
        if MX==2000:
            r_inj = 0.2
        else:
            r_inj = 1.0
        test_SigInj(working_area_SR, signal, rpf_order, params_to_set, r=r_inj, condor=True,scale_rpf=1.0)
        #test_SigInj_plot(working_area_SR,signal,rpf_order, r=0.0, condor=True)
        #test_SigInj_plot(working_area_SR,signal,rpf_order, r=0.2, condor=True)
    
    #Limit part
    for MX in [1200, 1400, 1600, 2000, 2500, 3000]:
        signal=f"MX{MX}_MY90"
        test_limit(working_area_SR,signal,rpf_order,'SR_run2.json',blind=True,rpf_params=params_to_set)

    signal=f"MX2000_MY90"
    par0 = params_to_set["Background_SR_rpf_0_par0"]
    test_Impacts(working_area_SR,signal,rpf_order,extra=f'-t -1 --setParameters Background_SR_rpf_0_par0={par0} --expectSignal=0.2')