import uproot
import pandas as pd
import time
import math
import numpy as np

import infofile


lumi = 10.0643 # 10 fb-1
fraction = 0.05
MC_to_data_ratio = 1
tuple_path = "/eos/project/a/atlas-outreach/projects/open-data/OpenDataTuples/renamedLargeRJets/2lep/" # local
#tuple_path = "https://atlas-opendata.web.cern.ch/atlas-opendata/samples/2020/2lep/" # web


samples = {

    'data': {
        'list' : ['data_A','data_B','data_C','data_D']
    },

    'HWW' : {
        'list' : ['VBFH125_WW2lep','ggH125_WW2lep'],#,'WpH125J_qqWW2lep','WpH125J_lvWW2lep','ZH125J_qqWW2lep','ZH125J_llWW2lep','ZH125J_vvWW2lep'],
        'color' : "#0074bf"
    },

    'WW' : {
        'list' : ['llvv'],#,'WpqqWmlv','WplvWmqq','WqqZll','WlvZqq','lllv','lvvv','ZqqZll','llll'],
        'color' : "#ff7400"
    },

    'ttbar' : {
        'list' : ['ttbar_lep'],#,'single_top_wtchan','single_antitop_wtchan','single_top_tchan','single_antitop_tchan','single_top_schan','single_antitop_schan'],
        'color' : "#00ab00"
    },

    'Z' : {
        'list' : ['Zee','Zmumu','Ztautau'],#,'Wplusenu','Wplusmunu','Wplustaunu','Wminusenu','Wminusmunu','Wminustaunu'],
        #'Zmumu_PTV0_70_CVetoBVeto','Zmumu_PTV0_70_CFilterBVeto','Zmumu_PTV0_70_BFilter','Zmumu_PTV70_140_CVetoBVeto','Zmumu_PTV70_140_CFilterBVeto','Zmumu_PTV70_140_BFilter','Zmumu_PTV140_280_CVetoBVeto','Zmumu_PTV140_280_CFilterBVeto','Zmumu_PTV140_280_BFilter','Zmumu_PTV280_500_CVetoBVeto','Zmumu_PTV280_500_CFilterBVeto','Zmumu_PTV280_500_BFilter','Zmumu_PTV500_1000','Zmumu_PTV1000_E_CMS','Zee_PTV0_70_CFilterBVeto','Zee_PTV0_70_BFilter','Zee_PTV70_140_CVetoBVeto','Zee_PTV70_140_CFilterBVeto','Zee_PTV70_140_BFilter','Zee_PTV140_280_CVetoBVeto','Zee_PTV140_280_CFilterBVeto','Zee_PTV140_280_BFilter','Zee_PTV280_500_CVetoBVeto','Zee_PTV280_500_CFilterBVeto','Zee_PTV280_500_BFilter','Zee_PTV500_1000','Zee_PTV1000_E_CMS','Ztautau_PTV0_70_CVetoBVeto','Ztautau_PTV0_70_CFilterBVeto','Ztautau_PTV0_70_BFilter','Ztautau_PTV70_140_CVetoBVeto','Ztautau_PTV70_140_CFilterBVeto','Ztautau_PTV70_140_BFilter','Ztautau_PTV140_280_CVetoBVeto','Ztautau_PTV140_280_CFilterBVeto','Ztautau_PTV140_280_BFilter','Ztautau_PTV280_500_CVetoBVeto','Ztautau_PTV280_500_CFilterBVeto','Ztautau_PTV280_500_BFilter','Ztautau_PTV500_1000'],
        'color' : "#ed0000"
    },

    #'Wt' : {
    #    'list' : ['single_top_wtchan','single_antitop_wtchan'],#,'single_top_tchan','single_antitop_tchan','single_top_schan','single_antitop_schan'],
    #    'color' : "#00c7f7"
    #}
    
}


def read_sample(s):
    print('Processing '+s+' samples')
    frames = []
    for val in samples[s]['list']:
        prefix = "MC/mc_"
        if s == 'data':
            prefix = "Data/"
        else: prefix += str(infofile.infos[val]["DSID"])+"."
        fileString = tuple_path+prefix+val+".2lep.root" 
        if fileString != "":
            temp = read_file(fileString,val)
            print('\t\tsum of weights',sum(temp['weight']))
            frames.append(temp)
        else:
            print("Error: "+val+" not found!")
    data_s = pd.concat(frames)
    data_s['LepDeltaPhi'] = data_s['LepDeltaPhi'].round(2)
    data_s.to_csv('13TeV_TransData.csv', mode='a', index=False, header=False)
    return data_s


def get_data_from_files():
    data = {}
    df=pd.DataFrame(columns=["type","NJets","MET","Mll","TransMass","LepDeltaPhi","METLLDeltaPhi","SumLepPt","BTags","weight"])
    df.to_csv('13TeV_TransData.csv',index=False)
    for s in samples:
        data[s] = read_sample(s)
    return data

# multiply event weights and scale factors
def calc_weight(mcWeight,scaleFactor_PILEUP,scaleFactor_ELE,
                scaleFactor_MUON, scaleFactor_LepTRIGGER):
    return mcWeight*scaleFactor_PILEUP*scaleFactor_ELE*scaleFactor_MUON*scaleFactor_LepTRIGGER


# multiply totalWeight by cross-section weight
def get_xsec_weight(totalWeight,sample):
    info = infofile.infos[sample]
    weight = (lumi*1000*info["xsec"])/(info["sumw"]*info["red_eff"]) #*1000 to go from fb-1 to pb-1
    weight *= totalWeight
    return round(weight/MC_to_data_ratio,5)

def find_good_lep_0_index(lep_n,lep_type,lep_pt,lep_eta,lep_ptcone,lep_etcone,lep_isTightID,lep_z0,lep_d0,lep_sigd0):
    for i in range(lep_n):
        if lep_pt[i]>15000 and abs(lep_eta[i])<2.5 and lep_ptcone[i]/lep_pt[i]<0.1 and lep_etcone[i]/lep_pt[i]<0.1 and abs(lep_d0[i])/lep_sigd0[i]<5 and lep_isTightID[i]:
            if lep_type[i]==13 and abs(lep_d0[i])/lep_sigd0[i]>3: continue
            theta_i = 2*math.atan(math.exp(-lep_eta[i]))
            if abs(lep_z0[i]*math.sin(theta_i))<0.5:
                return i
    return -1

def find_good_lep_1_index(lep_n,lep_type,lep_pt,lep_eta,lep_ptcone,lep_etcone,lep_isTightID,lep_z0,lep_d0,lep_sigd0,
                          good_lep_0_index):
    if good_lep_0_index!=-1:
        for i in range(good_lep_0_index+1,lep_n):
            if lep_pt[i]>15000 and abs(lep_eta[i])<2.5 and  lep_ptcone[i]/lep_pt[i]<0.1 and lep_etcone[i]/lep_pt[i]<0.1 and abs(lep_d0[i])/lep_sigd0[i]<5 and lep_isTightID[i]:
                if lep_type[i]==13 and abs(lep_d0[i])/lep_sigd0[i]>3: continue
                theta_i = 2*math.atan(math.exp(-lep_eta[i]))
                if abs(lep_z0[i]*math.sin(theta_i))<0.5:
                    return i
    return -1

def find_good_lep_2_index(lep_n,lep_type,lep_pt,lep_eta,lep_ptcone,lep_etcone,lep_isTightID,lep_z0,lep_d0,lep_sigd0,
                          good_lep_1_index):
    if good_lep_1_index!=-1:
        for i in range(good_lep_1_index+1,lep_n):
            if lep_pt[i]>15000 and abs(lep_eta[i])<2.5 and  lep_ptcone[i]/lep_pt[i]<0.1 and lep_etcone[i]/lep_pt[i]<0.1 and abs(lep_d0[i])/lep_sigd0[i]<5 and lep_isTightID[i]:
                if lep_type[i]==13 and abs(lep_d0[i])/lep_sigd0[i]>3: continue
                theta_i = 2*math.atan(math.exp(-lep_eta[i]))
                if abs(lep_z0[i]*math.sin(theta_i))<0.5:
                    return i
    return -1

def find_good_jet_0_index(jet_n,jet_pt,jet_eta,jet_jvt):
    for i in range(jet_n):
        if jet_pt[i]>25000:
            if jet_pt[i]<60000 and abs(jet_eta[i])<2.4:
                if jet_jvt[i]<0.59: continue
            return i
    return -1

def find_good_jet_1_index(jet_n,jet_pt,jet_eta,jet_jvt,good_jet_0_index):
    if good_jet_0_index!=-1:
        for i in range(good_jet_0_index+1,jet_n):
            if jet_pt[i]>25000:
                if jet_pt[i]<60000 and abs(jet_eta[i])<2.4:
                    if jet_jvt[i]<0.59: continue
                return i
    return -1

def find_good_jet_2_index(jet_n,jet_pt,jet_eta,jet_jvt,good_jet_1_index):
    if good_jet_1_index!=-1:
        for i in range(good_jet_1_index+1,jet_n):
            if jet_pt[i]>25000:
                if jet_pt[i]<60000 and abs(jet_eta[i])<2.4:
                    if jet_jvt[i]<0.59: continue
                return i
    return -1

def find_good_jet_3_index(jet_n,jet_pt,jet_eta,jet_jvt,good_jet_2_index):
    if good_jet_2_index!=-1:
        for i in range(good_jet_2_index+1,jet_n):
            if jet_pt[i]>25000:
                if jet_pt[i]<60000 and abs(jet_eta[i])<2.4:
                    if jet_jvt[i]<0.59: continue
                return i
    return -1

def find_good_jet_4_index(jet_n,jet_pt,jet_eta,jet_jvt,good_jet_3_index):
    if good_jet_3_index!=-1:
        for i in range(good_jet_3_index+1,jet_n):
            if jet_pt[i]>25000:
                if jet_pt[i]<60000 and abs(jet_eta[i])<2.4:
                    if jet_jvt[i]<0.59: continue
                return i
    return -1

def find_good_jet_5_index(jet_n,jet_pt,jet_eta,jet_jvt,good_jet_4_index):
    if good_jet_4_index!=-1:
        for i in range(good_jet_4_index+1,jet_n):
            if jet_pt[i]>25000:
                if jet_pt[i]<60000 and abs(jet_eta[i])<2.4:
                    if jet_jvt[i]<0.59: continue
                return i
    return -1

def find_good_jet_6_index(jet_n,jet_pt,jet_eta,jet_jvt,good_jet_5_index):
    if good_jet_5_index!=-1:
        for i in range(good_jet_5_index+1,jet_n):
            if jet_pt[i]>25000:
                if jet_pt[i]<60000 and abs(jet_eta[i])<2.4:
                    if jet_jvt[i]<0.59: continue
                return i
    return -1

def find_good_jet_7_index(jet_n,jet_pt,jet_eta,jet_jvt,good_jet_6_index):
    if good_jet_6_index!=-1:
        for i in range(good_jet_6_index+1,jet_n):
            if jet_pt[i]>25000:
                if jet_pt[i]<60000 and abs(jet_eta[i])<2.4:
                    if jet_jvt[i]<0.59: continue
                return i
    return -1

def find_good_jet_8_index(jet_n,jet_pt,jet_eta,jet_jvt,good_jet_7_index):
    if good_jet_7_index!=-1:
        for i in range(good_jet_7_index+1,jet_n):
            if jet_pt[i]>25000:
                if jet_pt[i]<60000 and abs(jet_eta[i])<2.4:
                    if jet_jvt[i]<0.59: continue
                return i
    return -1



# return number to represent which process
def mc_type(sample):
    if sample in samples['HWW']['list']: return 0
    elif sample in samples['WW']['list']: return 1
    elif sample in samples['ttbar']['list']: return 2
    elif sample in samples['Z']['list']: return 3
    #elif sample in samples['Wt']['list']: return 5
    else: return 4 #data

# return number to represent which channel
def channel(lep_type,good_lep_0_index,good_lep_1_index):
    if lep_type[good_lep_0_index]*lep_type[good_lep_1_index]==121: return 0 #ee
    elif lep_type[good_lep_0_index]*lep_type[good_lep_1_index]==169: return 1 #mm
    else: return 2 #em 

# calculate invariant mass of dilepton pair
def calc_mll(lep_pt,lep_eta,lep_phi,lep_E,good_lep_0_index,good_lep_1_index):
    px_0 = lep_pt[good_lep_0_index]*math.cos(lep_phi[good_lep_0_index])
    py_0 = lep_pt[good_lep_0_index]*math.sin(lep_phi[good_lep_0_index])
    pz_0 = lep_pt[good_lep_0_index]*math.sinh(lep_eta[good_lep_0_index])
    px_1 = lep_pt[good_lep_1_index]*math.cos(lep_phi[good_lep_1_index])
    py_1 = lep_pt[good_lep_1_index]*math.sin(lep_phi[good_lep_1_index])
    pz_1 = lep_pt[good_lep_1_index]*math.sinh(lep_eta[good_lep_1_index])
    sumpx = px_0 + px_1
    sumpy = py_0 + py_1
    sumpz = pz_0 + pz_1
    sump = math.sqrt(sumpx**2 + sumpy**2 + sumpz**2)
    sumE = lep_E[good_lep_0_index] + lep_E[good_lep_1_index]
    return round(math.sqrt(sumE**2 - sump**2)/1000,2) #/1000 to go from MeV to GeV

# calculate azimuthal angle difference between the 2 leptons
def calc_dPhiLL(lep_phi,good_lep_0_index,good_lep_1_index):
    dPhi = lep_phi[good_lep_0_index]-lep_phi[good_lep_1_index]
    if dPhi >= math.pi: dPhi -= 2*math.pi
    elif dPhi < -math.pi: dPhi += 2*math.pi
    return round(abs(dPhi),2)

# calculate azimuthal angle difference between the MET and the vector sum of the 2 leptons
def calc_dPhiLLmet(lep_pt,lep_phi,met_phi,good_lep_0_index,good_lep_1_index):
    px_0 = lep_pt[good_lep_0_index]*math.cos(lep_phi[good_lep_0_index])
    py_0 = lep_pt[good_lep_0_index]*math.sin(lep_phi[good_lep_0_index])
    px_1 = lep_pt[good_lep_1_index]*math.cos(lep_phi[good_lep_1_index])
    py_1 = lep_pt[good_lep_1_index]*math.sin(lep_phi[good_lep_1_index])
    sumpx = px_0 + px_1
    sumpy = py_0 + py_1
    phi_LL = math.atan2(sumpy,sumpx) # normal arctan but returns angle from -pi to +pi
    dPhi = phi_LL - met_phi
    if dPhi >= math.pi: dPhi -= 2*math.pi
    elif dPhi < -math.pi: dPhi += 2*math.pi
    return round(abs(dPhi),2)

# calculate the pt of the vector sum of the 2 leptons
def calc_ptLL(lep_pt,lep_phi,good_lep_0_index,good_lep_1_index):
    px_0 = lep_pt[good_lep_0_index]*math.cos(lep_phi[good_lep_0_index])
    py_0 = lep_pt[good_lep_0_index]*math.sin(lep_phi[good_lep_0_index])
    px_1 = lep_pt[good_lep_1_index]*math.cos(lep_phi[good_lep_1_index])
    py_1 = lep_pt[good_lep_1_index]*math.sin(lep_phi[good_lep_1_index])
    sumpx = px_0 + px_1
    sumpy = py_0 + py_1
    return round(math.sqrt(sumpx**2 + sumpy**2)/1000,2) #/1000 to go from MeV to GeV 
    
# calculate transverse mass
def calc_Mt(lep_pt,lep_eta,lep_E,met_et,good_lep_0_index,good_lep_1_index):
    E = (lep_E[good_lep_0_index]+lep_E[good_lep_1_index]+met_et)
    pz_0 = lep_pt[good_lep_0_index]*math.sinh(lep_eta[good_lep_0_index])
    pz_1 = lep_pt[good_lep_1_index]*math.sinh(lep_eta[good_lep_1_index])
    sumpz = pz_0 + pz_1
    Mt2 = E**2 - sumpz**2
    Mt = math.sqrt(Mt2)/1000
    return round(Mt,2)

# determine whether any jets are bjets
def bjets(jet_MV2c10,
              good_jet_0_index,good_jet_1_index,good_jet_2_index,good_jet_3_index,good_jet_4_index,
                      good_jet_5_index,good_jet_6_index,good_jet_7_index,good_jet_8_index):
    bjets = 0
    all_jets_indices = [good_jet_0_index,good_jet_1_index,good_jet_2_index,good_jet_3_index,good_jet_4_index,good_jet_5_index,good_jet_6_index,good_jet_7_index,good_jet_8_index]
    good_jets_indices = [jet_i for jet_i in all_jets_indices if jet_i!=-1]
    for i in good_jets_indices:
        if jet_MV2c10[i]>0.1758475:
            bjets = 1
    return bjets


# throw away events that don't have 2 leptons
def cut_good_lep_n(good_lep_1_index,good_lep_2_index):
    # return when number of good leptons is not equal to 2                                                               
    # good_lep_index_1==-1 means there's no 2nd good lepton                                                              
    # lep_index_2!=-1 means there's a 3rd lepton                                                                         
    return good_lep_1_index==-1 or good_lep_2_index!=-1

# throw away events that don't have opposite-charge leptons
def cut_lep_charge(lep_charge,good_lep_0_index,good_lep_1_index):
    # return when sum of lepton charges is not equal to 0                                                                
    # first lepton is [0], 2nd lepton is [1]                                                                             
    return lep_charge[good_lep_0_index] + lep_charge[good_lep_1_index] != 0

def calc_good_jet_n(jet_pt,good_jet_0_index,good_jet_1_index,good_jet_2_index,good_jet_3_index,good_jet_4_index,
                      good_jet_5_index,good_jet_6_index,good_jet_7_index,good_jet_8_index):
    all_jets_indices = [good_jet_0_index,good_jet_1_index,good_jet_2_index,good_jet_3_index,good_jet_4_index,good_jet_5_index,good_jet_6_index,good_jet_7_index,good_jet_8_index]
    good_jets_indices = [jet_i for jet_i in all_jets_indices if jet_i!=-1]
    pt30_jets_indices = [jet_i for jet_i in good_jets_indices if jet_pt[jet_i]>30000]
    return len(pt30_jets_indices)

def cut_Mt_lower(TransMass):
    return TransMass<50

def cut_Mt_upper(TransMass):
    return TransMass>300

# throw away events where Mll < 10 GeV
def Mll_cut_lower(Mll):
    return Mll<10

def Mll_cut_upper(Mll):
    return Mll>105

def cut_SumLepPt(SumLepPt):
    return SumLepPt>200

def cut_met_et(met_et):
    return met_et>200*1000

def cut_LepDeltaPhi(LepDeltaPhi):
    return LepDeltaPhi==0

def cut_METLLDeltaPhi(METLLDeltaPhi):
    return METLLDeltaPhi==0

# throw away events where weight is 0
def cut_weight(weight):
    return weight<0.00005

def cut_weight_upper(weight):
    return weight>=1

# throw away events which aren't emu channel
def channel_cut(lep_type,good_lep_0_index,good_lep_1_index):
    return lep_type[good_lep_0_index]*lep_type[good_lep_1_index]!=143

def read_file(path,sample):
    start = time.time()
    print("\tProcessing: "+sample+" file")
    data_all = pd.DataFrame()
    mc = uproot.open(path)["mini"]
    numevents = uproot.numentries(path, "mini")
    if 'data' in sample: fraction_MC=fraction
    else: fraction_MC=fraction*MC_to_data_ratio
    entrystart=0
    entrystop=numevents*fraction_MC
    for data in mc.iterate(["lep_n","lep_pt","lep_eta","lep_phi","lep_E","lep_charge","lep_type",
                            "lep_isTightID","lep_ptcone30","lep_etcone20","lep_z0","lep_trackd0pvunbiased","lep_tracksigd0pvunbiased",
                            "jet_n","jet_pt","jet_eta","jet_jvt","jet_MV2c10","met_et","met_phi",
                         "mcWeight","scaleFactor_PILEUP","scaleFactor_ELE","scaleFactor_MUON", # add more variables here if you make cuts on them ,              
                            "scaleFactor_LepTRIGGER"], flatten=False, entrysteps=2459370, outputtype=pd.DataFrame, entrystart=entrystart, entrystop=entrystop):

        nIn = len(data.index)

        data['good_lep_0_index'] = np.vectorize(find_good_lep_0_index)(data.lep_n,data.lep_type,data.lep_pt,data.lep_eta,
                                                                       data.lep_ptcone30,data.lep_etcone20,data.lep_isTightID,
                                                                       data.lep_z0,data.lep_trackd0pvunbiased,data.lep_tracksigd0pvunbiased)
        data['good_lep_1_index'] = np.vectorize(find_good_lep_1_index)(data.lep_n,data.lep_type,data.lep_pt,data.lep_eta,
                                                                       data.lep_ptcone30,data.lep_etcone20,data.lep_isTightID,
                                                                       data.lep_z0,data.lep_trackd0pvunbiased,data.lep_tracksigd0pvunbiased,data.good_lep_0_index)
        data['good_lep_2_index'] = np.vectorize(find_good_lep_2_index)(data.lep_n,data.lep_type,data.lep_pt,data.lep_eta,
                                                                       data.lep_ptcone30,data.lep_etcone20,data.lep_isTightID,
                                                                       data.lep_z0,data.lep_trackd0pvunbiased,data.lep_tracksigd0pvunbiased,data.good_lep_1_index)
        data['good_jet_0_index'] = np.vectorize(find_good_jet_0_index)(data.jet_n,data.jet_pt,
                                                                       data.jet_eta,data.jet_jvt)
        data['good_jet_1_index'] = np.vectorize(find_good_jet_1_index)(data.jet_n,data.jet_pt,
                                                                       data.jet_eta,data.jet_jvt,
                                                                      data.good_jet_0_index)
        data['good_jet_2_index'] = np.vectorize(find_good_jet_2_index)(data.jet_n,data.jet_pt,
                                                                       data.jet_eta,data.jet_jvt,
                                                                      data.good_jet_1_index)
        data['good_jet_3_index'] = np.vectorize(find_good_jet_3_index)(data.jet_n,data.jet_pt,
                                                                       data.jet_eta,data.jet_jvt,
                                                                      data.good_jet_2_index)
        data['good_jet_4_index'] = np.vectorize(find_good_jet_4_index)(data.jet_n,data.jet_pt,
                                                                       data.jet_eta,data.jet_jvt,
                                                                      data.good_jet_3_index)
        data['good_jet_5_index'] = np.vectorize(find_good_jet_5_index)(data.jet_n,data.jet_pt,
                                                                       data.jet_eta,data.jet_jvt,
                                                                      data.good_jet_4_index)
        data['good_jet_6_index'] = np.vectorize(find_good_jet_6_index)(data.jet_n,data.jet_pt,
                                                                       data.jet_eta,data.jet_jvt,
                                                                      data.good_jet_5_index)
        data['good_jet_7_index'] = np.vectorize(find_good_jet_7_index)(data.jet_n,data.jet_pt,
                                                                       data.jet_eta,data.jet_jvt,
                                                                      data.good_jet_6_index)
        data['good_jet_8_index'] = np.vectorize(find_good_jet_8_index)(data.jet_n,data.jet_pt,data.jet_eta,data.jet_jvt,
                                                                      data.good_jet_7_index)

        # throw away events where number of leptons isn't 2
        fail = data[ np.vectorize(cut_good_lep_n)(data.good_lep_1_index,data.good_lep_2_index) ].index
        data.drop(fail, inplace=True)

        # throw away events where leptons aren't oppositely charged
        fail = data[ np.vectorize(cut_lep_charge)(data.lep_charge,data.good_lep_0_index,data.good_lep_1_index) ].index
        data.drop(fail, inplace=True)

        # cut on channel
        fail = data[ np.vectorize(channel_cut)(data.lep_type,data.good_lep_0_index,data.good_lep_1_index) ].index
        data.drop(fail, inplace=True)

        fail = data[ np.vectorize(cut_met_et)(data.met_et) ].index
        data.drop(fail, inplace=True)

        # label for each mc type
        data['type'] = np.vectorize(mc_type)(sample)

        # label for channel (ee, mm or em)
        #data['Channel'] = np.vectorize(channel)(data.lep_type,data.good_lep_0_index,data.good_lep_1_index)

        # number of jets
        data['NJets'] = np.vectorize(calc_good_jet_n)(data.jet_pt,data.good_jet_0_index,data.good_jet_1_index,data.good_jet_2_index,
                                          data.good_jet_3_index,data.good_jet_4_index,data.good_jet_5_index,
                                          data.good_jet_6_index,data.good_jet_7_index,data.good_jet_8_index)

        # MET
        data['MET'] = round(data['met_et']/1000,2)

        # calculation of 2-lepton invariant mass
        data['Mll'] = np.vectorize(calc_mll)(data.lep_pt,data.lep_eta,data.lep_phi,data.lep_E,data.good_lep_0_index,
                                            data.good_lep_1_index)
        fail = data[ np.vectorize(Mll_cut_lower)(data.Mll) ].index
        data.drop(fail, inplace=True)
        fail = data[ np.vectorize(Mll_cut_upper)(data.Mll) ].index
        data.drop(fail, inplace=True)

        # transverse mass
        data['TransMass'] = np.vectorize(calc_Mt)(data.lep_pt,data.lep_eta,data.lep_E,data.met_et,data.good_lep_0_index,data.good_lep_1_index)
        fail = data[ np.vectorize(cut_Mt_lower)(data.TransMass) ].index
        data.drop(fail, inplace=True)
        fail = data[ np.vectorize(cut_Mt_upper)(data.TransMass) ].index
        data.drop(fail, inplace=True)

        # Angular separation between leptons
        data['LepDeltaPhi'] = np.vectorize(calc_dPhiLL)(data.lep_phi,data.good_lep_0_index,data.good_lep_1_index)
        fail = data[ np.vectorize(cut_LepDeltaPhi)(data.LepDeltaPhi) ].index
        data.drop(fail, inplace=True)

        # Angular separation between leptons and MET dPhi(MET,ll)
        data['METLLDeltaPhi'] = np.vectorize(calc_dPhiLLmet)(data.lep_pt,data.lep_phi,data.met_phi,data.good_lep_0_index,data.good_lep_1_index)
        fail = data[ np.vectorize(cut_METLLDeltaPhi)(data.METLLDeltaPhi) ].index
        data.drop(fail, inplace=True)

        # Sum of lepton pt
        data['SumLepPt'] = np.vectorize(calc_ptLL)(data.lep_pt,data.lep_phi,data.good_lep_0_index,
                                             data.good_lep_1_index        )
        fail = data[ np.vectorize(cut_SumLepPt)(data.SumLepPt) ].index
        data.drop(fail, inplace=True)

        # whether at least 1 jet is btagged
        data['BTags'] = np.vectorize(bjets)(data.jet_MV2c10,
                                        data.good_jet_0_index,data.good_jet_1_index,data.good_jet_2_index,
                                          data.good_jet_3_index,data.good_jet_4_index,data.good_jet_5_index,
                                          data.good_jet_6_index,data.good_jet_7_index,data.good_jet_8_index        )

        if 'data' not in sample:
            data['weight'] = np.vectorize(calc_weight)(data.mcWeight,data.scaleFactor_PILEUP,data.scaleFactor_ELE,data.scaleFactor_MUON,data.scaleFactor_LepTRIGGER)
            data['weight'] = np.vectorize(get_xsec_weight)(data.weight,sample)
            # throw away events with weight 0
            fail = data[ np.vectorize(cut_weight)(data.weight) ].index
            data.drop(fail, inplace=True)
            fail = data[ np.vectorize(cut_weight_upper)(data.weight) ].index
            data.drop(fail, inplace=True)
        else:
            data['weight'] = 1

        data.drop(["lep_n","lep_pt","lep_eta","lep_phi","lep_E","lep_charge","lep_type","lep_isTightID","lep_ptcone30","lep_etcone20","lep_z0","lep_trackd0pvunbiased","lep_tracksigd0pvunbiased","jet_n","jet_pt","jet_eta","jet_jvt","jet_MV2c10","met_et","met_phi","mcWeight","scaleFactor_PILEUP","scaleFactor_ELE","scaleFactor_MUON","scaleFactor_LepTRIGGER","good_lep_0_index","good_lep_1_index","good_lep_2_index","good_jet_0_index","good_jet_1_index","good_jet_2_index","good_jet_3_index","good_jet_4_index","good_jet_5_index","good_jet_6_index","good_jet_7_index","good_jet_8_index"], axis=1, inplace=True)

        #print(data[['LepDeltaPhi','METLLDeltaPhi','TransMass','weight']])
        #print(data['LepDeltaPhi'])

        nOut = len(data.index)
        data_all = data_all.append(data)
        elapsed = time.time() - start
        print("\t\t"+sample+" time taken: "+str(elapsed)+"s, nIn: "+str(nIn)+", nOut: "+str(nOut))

    return data_all


start = time.time()
data = get_data_from_files()
elapsed = time.time() - start
print("Time taken: "+str(elapsed))

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator,LogLocator,LogFormatterSciNotation # for minor ticks
import HWWHistograms
import labelfile
stack_order = ['WW','ttbar','Z']

def plot_data(data):

    signal_format = 'hist' # 'line' for line above SM stack
                           # 'hist' for bar above SM stack
                           # None for signal as part of SM stack
    Total_SM_label = False # for Total SM black line in plot and legend
    plot_label = r'$H \rightarrow WW \rightarrow e\nu\mu\nu$'
    signal_label = r'Signal ($m_H=125$ GeV)' # r''

    # *******************
    # general definitions (shouldn't need to change)
    lumi_used = str(lumi*fraction)    
    signal = None
    for s in samples.keys():
        if s not in stack_order and s!='data': signal = s

    for x_variable,hist in HWWHistograms.hist_dict.items():

        h_bin_width = hist['bin_width']
        h_num_bins = hist['num_bins']
        h_xrange_min = hist['xrange_min']
        h_log_y = hist['log_y']
        h_y_label_x_position = hist['y_label_x_position']
        h_legend_loc = hist['legend_loc']
        h_log_top_margin = hist['log_top_margin'] # to decrease the separation between data and the top of the figure, remove a 0
        h_linear_top_margin = hist['linear_top_margin'] # to decrease the separation between data and the top of the figure, pick a number closer to 1
    
        bins = [h_xrange_min + x*h_bin_width for x in range(h_num_bins+1) ]
        bin_centres = [h_xrange_min+h_bin_width/2 + x*h_bin_width for x in range(h_num_bins) ]

        data_x,_ = np.histogram(data['data'][x_variable].values, bins=bins)
        data_x_errors = np.sqrt(data_x)

        signal_x = None
        if signal_format=='line':
            signal_x,_ = np.histogram(data[signal][x_variable].values,bins=bins,weights=data[signal].weight.values)
        elif signal_format=='hist':
            signal_x = data[signal][x_variable].values
            signal_weights = data[signal].weight.values
            signal_color = samples[signal]['color']
    
        mc_x = []
        mc_weights = []
        mc_colors = []
        mc_labels = []
        mc_x_tot = np.zeros(len(bin_centres))

        for s in stack_order:
            mc_labels.append(s)
            mc_x.append(data[s][x_variable].values)
            mc_colors.append(samples[s]['color'])
            mc_weights.append(data[s].weight.values)
            mc_x_heights,_ = np.histogram(data[s][x_variable].values,bins=bins,weights=data[s].weight.values)
            mc_x_tot = np.add(mc_x_tot, mc_x_heights)
    
        mc_x_err = np.sqrt(mc_x_tot)
    
    
        # *************
        # Main plot 
        # *************
        plt.clf()
        plt.axes([0.1,0.3,0.85,0.65]) #(left, bottom, width, height)
        main_axes = plt.gca()
        main_axes.errorbar( x=bin_centres, y=data_x, yerr=data_x_errors, fmt='ko', label='Data')
        mc_heights = main_axes.hist(mc_x,bins=bins,weights=mc_weights,stacked=True,color=mc_colors, label=mc_labels)
        if Total_SM_label:
            totalSM_handle, = main_axes.step(bins,np.insert(mc_x_tot,0,mc_x_tot[0]),color='black')
        if signal_format=='line':
            main_axes.step(bins,np.insert(signal_x,0,signal_x[0]),color=samples[signal]['color'], linestyle='--',
                       label=signal)
        elif signal_format=='hist':
            main_axes.hist(signal_x,bins=bins,bottom=mc_x_tot,weights=signal_weights,color=signal_color,label=signal)
        main_axes.bar(bin_centres,2*mc_x_err,bottom=mc_x_tot-mc_x_err,alpha=0.5,color='none',hatch="////",
                  width=h_bin_width, label='Stat. Unc.')
        
        main_axes.set_xlim(left=h_xrange_min,right=bins[-1])
        main_axes.xaxis.set_minor_locator(AutoMinorLocator()) # separation of x axis minor ticks
        main_axes.tick_params(which='both',direction='in',top=True,labeltop=False,labelbottom=False,right=True,labelright=False)
        main_axes.set_ylabel(r'Events / '+str(h_bin_width)+r' GeV',fontname='sans-serif',horizontalalignment='right',y=1.0,fontsize=11)
        if h_log_y:
            main_axes.set_yscale('log')
            smallest_contribution = mc_heights[0][0]
            smallest_contribution.sort()
            bottom = smallest_contribution[-2]
            top = np.amax(data_x)*h_log_top_margin
            main_axes.set_ylim(bottom=bottom,top=top)
            main_axes.yaxis.set_major_formatter(CustomTicker())
            locmin = LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=12)
            main_axes.yaxis.set_minor_locator(locmin)
        else: 
            main_axes.set_ylim(bottom=0,top=(np.amax(data_x)+math.sqrt(np.amax(data_x)))*h_linear_top_margin)
            main_axes.yaxis.set_minor_locator(AutoMinorLocator())
        
        plt.text(0.05,0.97,'ATLAS',ha="left",va="top",family='sans-serif',transform=main_axes.transAxes,style='italic',weight='bold',fontsize=13)
        plt.text(0.19,0.97,'Open Data',ha="left",va="top",family='sans-serif',transform=main_axes.transAxes,fontsize=13)
        plt.text(0.05,0.9,'for education only',ha="left",va="top",family='sans-serif',transform=main_axes.transAxes,style='italic',fontsize=8)
        plt.text(0.05,0.86,r'$\sqrt{s}=13\,\mathrm{TeV},\;\int L\,dt=$'+lumi_used+'$\,\mathrm{fb}^{-1}$',ha="left",va="top",family='sans-serif',transform=main_axes.transAxes)
        plt.text(0.05,0.78,plot_label,ha="left",va="top",family='sans-serif',transform=main_axes.transAxes)
    
        # Create new legend handles but use the colors from the existing ones 
        handles, labels = main_axes.get_legend_handles_labels()
        if signal_format=='line':
            handles[labels.index(signal)] = Line2D([], [], c=samples[signal]['color'], linestyle='dashed')
        if Total_SM_label:
            uncertainty_handle = mpatches.Patch(facecolor='none',hatch='////')
            handles.append((totalSM_handle,uncertainty_handle))
            labels.append('Total SM')
    
        # specify order within legend
        new_handles = [handles[labels.index('Data')]]
        new_labels = ['Data']
        for s in reversed(stack_order):
            new_handles.append(handles[labels.index(s)])
            new_labels.append(s)
        if Total_SM_label:
            new_handles.append(handles[labels.index('Total SM')])
            new_labels.append('Total SM')
        else: 
            new_handles.append(handles[labels.index('Stat. Unc.')])
            new_labels.append('Stat. Unc.')
        if signal is not None:
            new_handles.append(handles[labels.index(signal)])
            new_labels.append(signal_label)
        main_axes.legend(handles=new_handles, labels=new_labels, frameon=False, loc=h_legend_loc)
    
    
        # *************
        # Data/MC ratio 
        # *************
        plt.axes([0.1,0.1,0.85,0.2]) #(left, bottom, width, height)
        ratio_axes = plt.gca()
        ratio_axes.errorbar( x=bin_centres, y=data_x/mc_x_tot, yerr=data_x_errors/mc_x_tot, fmt='ko')
        ratio_axes.bar(bin_centres,2*mc_x_err/mc_x_tot,bottom=1-mc_x_err/mc_x_tot,alpha=0.5,color='none',
            hatch="////",width=h_bin_width)
        ratio_axes.plot(bins,np.ones(len(bins)),color='k')
        ratio_axes.set_xlim(left=h_xrange_min,right=bins[-1])
        ratio_axes.xaxis.set_minor_locator(AutoMinorLocator()) # separation of x axis minor ticks
        ratio_axes.xaxis.set_label_coords(0.9,-0.2) # (x,y) of x axis label # 0.2 down from x axis
        ratio_axes.set_xlabel(labelfile.variable_labels[x_variable],fontname='sans-serif',fontsize=11)
        ratio_axes.tick_params(which='both',direction='in',top=True,labeltop=False,right=True,labelright=False)
        ratio_axes.set_ylim(bottom=0,top=2.5)
        ratio_axes.set_yticks([0,1,2])
        ratio_axes.yaxis.set_minor_locator(AutoMinorLocator())
        if signal is not None:
            ratio_axes.set_ylabel(r'Data/SM',fontname='sans-serif',x=1,fontsize=11)
        else:
            ratio_axes.set_ylabel(r'Data/MC',fontname='sans-serif',fontsize=11)
        
        
        # Generic features for both plots
        main_axes.yaxis.set_label_coords(h_y_label_x_position,1)
        ratio_axes.yaxis.set_label_coords(h_y_label_x_position,0.5)
    
        plt.savefig("HWW_"+x_variable+".pdf")
    
    return signal_x,mc_x_tot

#signal_yields,background_yields = plot_data(data)
