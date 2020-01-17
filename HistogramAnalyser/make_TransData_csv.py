import uproot
import pandas as pd
import time
import math
import numpy as np

import infofile


lumi = 10 # 10 fb-1
fraction = 1
tuple_path = "/eos/project/a/atlas-outreach/projects/open-data/OpenDataTuples/renamedLargeRJets/exactly2lep/" # local


samples = {

    'data': {
        'list' : ['data_A','data_B','data_C','data_D']
    },

    'HWW' : {
        'list' : ['VBFH125_WW2lep','ggH125_WW2lep','WpH125J_qqWW2lep','ZH125J_qqWW2lep','ZH125J_vvWW2lep'],
    },

    'WW' : {
        'list' : ['WpqqWmlv','WplvWmqq','llvv'],
    },

    'ttbar' : {
        'list' : ['ttbar_lep'],
    },

    'Z' : {
        'list' : ['Zee','Zmumu','Ztautau'],
    }

}


def read_sample(s):
    print('Processing '+s+' samples')
    frames = []
    for val in samples[s]['list']:
        prefix = "MC/mc_"
        if s == 'data':
            prefix = "Data/"
        else: prefix += str(infofile.infos[val]["DSID"])+"."
        fileString = tuple_path+prefix+val+".exactly2lep.root" 
        if fileString != "":
            temp = read_file(fileString,val)
            frames.append(temp)
        else:
            print("Error: "+val+" not found!")
    data_s = pd.concat(frames)
    data_s.to_csv('13TeVTransData.csv', mode='a', index=False, header=False)
    return data_s


def get_data_from_files():

    data = {}
    df=pd.DataFrame(columns=["type","Channel","NJets","MET","Mll","TransverseMass","TransMass","LepDeltaPhi","METLLDeltaPhi","SumLepPt","BTags","weight"])
    df.to_csv('13TeVTransData.csv',index=False)
    for s in samples:
        data[s] = read_sample(s)
    
    return data


def calc_weight(mcWeight,scaleFactor_PILEUP,scaleFactor_ELE,
                scaleFactor_MUON, scaleFactor_LepTRIGGER):
    return mcWeight*scaleFactor_PILEUP*scaleFactor_ELE*scaleFactor_MUON*scaleFactor_LepTRIGGER


def get_xsec_weight(totalWeight,sample):
    info = infofile.infos[sample]
    weight = (lumi*1000*info["xsec"])/(info["sumw"]*info["red_eff"]) #*1000 to go from fb-1 to pb-1
    weight *= totalWeight
    return weight


def mc_type(sample):
    if sample in samples['HWW']['list']: return 0
    elif sample in samples['WW']['list']: return 1
    elif sample in samples['ttbar']['list']: return 2
    elif sample in samples['Z']['list']: return 3
    else: return 4 # data

def channel_cut(lep_type):
    return lep_type[0]*lep_type[1]!=143

def channel(lep_type):
    if lep_type[0]*lep_type[1]==121: return 0 #ee
    elif lep_type[0]*lep_type[1]==169: return 1 #mm
    else: return 2 #em

def NJets(jet_n):
    return jet_n

def MET(met_et):
    return met_et/1000

def calc_dPhiLL(lep_phi):
    return abs(lep_phi[0]-lep_phi[1])

def calc_dPhiLLmet(lep_pts,lep_etas,lep_phis,met_phi):
    theta_0 = 2*math.atan(math.exp(-lep_etas[0]))
    theta_1 = 2*math.atan(math.exp(-lep_etas[1]))
    p_0 = lep_pts[0]/math.sin(theta_0)
    p_1 = lep_pts[1]/math.sin(theta_1)
    px_0 = p_0*math.sin(theta_0)*math.cos(lep_phis[0])
    px_1 = p_1*math.sin(theta_1)*math.cos(lep_phis[1])
    py_0 = p_0*math.sin(theta_0)*math.sin(lep_phis[0])
    py_1 = p_1*math.sin(theta_1)*math.sin(lep_phis[1])
    sumpx = px_0 + px_1
    sumpy = py_0 + py_1
    sumpt = math.sqrt(sumpx**2 + sumpy**2)
    phi_LL = np.sign(sumpy)*math.acos(sumpx/sumpt)
    return abs(phi_LL-met_phi)

def calc_ptLL(lep_pts,lep_etas,lep_phis):
    theta_0 = 2*math.atan(math.exp(-lep_etas[0]))
    theta_1 = 2*math.atan(math.exp(-lep_etas[1]))
    p_0 = lep_pts[0]/math.sin(theta_0)
    p_1 = lep_pts[1]/math.sin(theta_1)
    px_0 = p_0*math.sin(theta_0)*math.cos(lep_phis[0])
    px_1 = p_1*math.sin(theta_1)*math.cos(lep_phis[1])
    py_0 = p_0*math.sin(theta_0)*math.sin(lep_phis[0])
    py_1 = p_1*math.sin(theta_1)*math.sin(lep_phis[1])
    sumpx = px_0 + px_1
    sumpy = py_0 + py_1
    return math.sqrt(sumpx**2 + sumpy**2)/1000
    
def calc_mT(ptLL,met_et,dPhiLLmet):
    return math.sqrt(2*ptLL*met_et*(1-math.cos(dPhiLLmet)))/1000

def bjet_n(jet_n,jet_MV2c10):
    bjet_n = 0
    for i in range(jet_n):
        if jet_MV2c10[i]>0.1758475: 
            bjet_n+=1
    return bjet_n


def calc_mll(lep_pts,lep_etas,lep_phis):
    mll = 2*lep_pts[0]*lep_pts[1]
    cosh = math.cosh(lep_etas[0]-lep_etas[1])
    cos = math.cos(lep_phis[0]-lep_phis[1])
    mll *= ( cosh - cos )
    return math.sqrt(mll)/1000


def read_file(path,sample):
    start = time.time()
    print("\tProcessing: "+sample+" file")
    data_all = pd.DataFrame()
    mc = uproot.open(path)["mini"]
    numevents = uproot.numentries(path, "mini")
    for data in mc.iterate(["lep_pt","lep_eta","lep_phi","lep_type",
                            "jet_n","jet_MV2c10","met_et","met_phi",
                         "mcWeight","scaleFactor_PILEUP","scaleFactor_ELE","scaleFactor_MUON", # add more variables here if you make cuts on them ,  
                            "scaleFactor_LepTRIGGER"], flatten=False, entrysteps=2500000, outputtype=pd.DataFrame, entrystop=numevents*fraction):

        nIn = len(data.index)

        # label for each mc type
        data['type'] = np.vectorize(mc_type)(sample)

        # cut on lepton type
        fail = data[ np.vectorize(channel_cut)(data.lep_type) ].index
        data.drop(fail, inplace=True)

        # label for channel (ee, mm or em)                                                                            
        data['Channel'] = np.vectorize(channel)(data.lep_type)

        # number of jets
        data['NJets'] = data['jet_n']

        # MET
        data['MET'] = data['met_et']/1000

        # calculation of 2-lepton invariant mass
        data['Mll'] = np.vectorize(calc_mll)(data.lep_pt,data.lep_eta,data.lep_phi)

        # Angular separation between leptons
        data['LepDeltaPhi'] = np.vectorize(calc_dPhiLL)(data.lep_phi)

        # Angular separation between leptons and MET dPhi(MET,ll)
        data['METLLDeltaPhi'] = np.vectorize(calc_dPhiLLmet)(data.lep_pt,data.lep_eta,data.lep_phi,data.met_phi)

        # Sum of lepton pt
        data['SumLepPt'] = np.vectorize(calc_ptLL)(data.lep_pt,data.lep_eta,data.lep_phi)

        # transverse mass
        data['TransMass'] = np.vectorize(calc_mT)(data.SumLepPt,data.met_et,data.METLLDeltaPhi)

        # transverse mass
        data['TransverseMass'] = np.vectorize(calc_mT)(data.SumLepPt,data.met_et,data.METLLDeltaPhi)

        # calculation of bjet_n
        data['BTags'] = np.vectorize(bjet_n)(data.jet_n,data.jet_MV2c10)

        if 'data' not in sample:
            data['weight'] = np.vectorize(calc_weight)(data.mcWeight,data.scaleFactor_PILEUP,data.scaleFactor_ELE,data.scaleFactor_MUON,data.scaleFactor_LepTRIGGER)
            data['weight'] = np.vectorize(get_xsec_weight)(data.weight,sample)
        else:
            data['weight'] = 1

        data.drop(["lep_pt","lep_eta","lep_phi","lep_type","jet_n","jet_MV2c10","met_et","met_phi","mcWeight","scaleFactor_PILEUP","scaleFactor_ELE","scaleFactor_MUON","scaleFactor_LepTRIGGER"], axis=1, inplace=True)

        data = data[data.weight != 0]

        #print(data[['lep_eta']])
        #print(data)

        nOut = len(data.index)
        data_all = data_all.append(data)
        elapsed = time.time() - start
        print("\t\t"+sample+" time taken: "+str(elapsed)+"s, nIn: "+str(nIn)+", nOut: "+str(nOut))
    
    return data_all


start = time.time()
data = get_data_from_files()
elapsed = time.time() - start
print("Time taken: "+str(elapsed))
