import uproot
import pandas as pd
import time
import math
import numpy as np

import infofile


lumi = 10 # 10 fb-1
fraction = 1
tuple_path = "/eos/project/a/atlas-outreach/projects/open-data/OpenDataTuples/renamedLargeRJets/4lep/" # local


samples = {

    #'data': {
    #    'list' : ['data_A','data_B','data_C','data_D']
    #},

    'HZZ' : {
        'list' : ['ZH125_ZZ4lep','WH125_ZZ4lep','VBFH125_ZZ4lep','ggH125_ZZ4lep'],
    },

    'ZZ' : {
        'list' : ['llll'],
    },

    'ttbar' : {
        'list' : ['ttbar_lep'],
    },

    'Z' : {
        'list' : ['Zee','Zmumu'],
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
        fileString = tuple_path+prefix+val+".4lep.root" 
        if fileString != "":
            temp = read_file(fileString,val)
            frames.append(temp)
        else:
            print("Error: "+val+" not found!")
    data_s = pd.concat(frames)
    data_s.to_csv('13TeV_HZZ.csv', mode='a', index=False, header=False)
    return data_s


def get_data_from_files():

    data = {}
    df=pd.DataFrame(columns=["type","Channel","SumCharges","MinmOCST","LepPt0","LepPt1","LepPt2","MZ1","MZ2","Mllll","weight"])
    df.to_csv('13TeV_HZZ.csv',index=False)
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
    if sample in samples['HZZ']['list']: return 0
    elif sample in samples['ZZ']['list']: return 1
    elif sample in samples['ttbar']['list']: return 2
    elif sample in samples['Z']['list']: return 3
    else: return 4 # data

def channel(lep_type):
    if lep_type[0]+lep_type[1]+lep_type[2]+lep_type[3]==44: return 0 #eeee
    elif lep_type[0]+lep_type[1]+lep_type[2]+lep_type[3]==52: return 1 #mmmm
    elif lep_type[0]+lep_type[1]+lep_type[2]+lep_type[3]==48: return 2 #eemm
    elif lep_type[0]+lep_type[1]+lep_type[2]+lep_type[3]==46: return 3 #eeem
    else: return 4 #emmm

def sum_charge(lep_charge):
    return lep_charge[0]+lep_charge[1]+lep_charge[2]+lep_charge[3]

def lep_pt_0(lep_pt):
    return lep_pt[0]/1000

def lep_pt_1(lep_pt):
    return lep_pt[1]/1000

def lep_pt_2(lep_pt):
    return lep_pt[2]/1000

def calc_min_mOCST(lep_pts,lep_etas,lep_phis,lep_charges,lep_types):
    # same-type-opposite-flavour
    mOCST = []
    for i in range(len(lep_pts)-1):
        for j in range(len(lep_pts)):
            if j>i and lep_charges[i]!=lep_charges[j] and lep_types[i]==lep_types[j]:
                mll = 2*lep_pts[i]*lep_pts[j]
                cosh = math.cosh(lep_etas[i]-lep_etas[j])
                cos = math.cos(lep_phis[i]-lep_phis[j])
                mll *= ( cosh - cos )
                mOCST.append(math.sqrt(mll)/1000)
    if len(mOCST)>0: return min(mOCST)
    else: return 0

def calc_m_Z1(lep_pts,lep_etas,lep_phis,lep_charges,lep_types):
    # mass of Z boson candidate 1
    mOCST = []
    for i in range(len(lep_pts)-1):
        for j in range(len(lep_pts)):
            if j>i and lep_charges[i]!=lep_charges[j] and lep_types[i]==lep_types[j]:
                mll = 2*lep_pts[i]*lep_pts[j]
                cosh = math.cosh(lep_etas[i]-lep_etas[j])
                cos = math.cos(lep_phis[i]-lep_phis[j])
                mll *= ( cosh - cos )
                mOCST.append(math.sqrt(mll)/1000.)      
    if len(mOCST)>0: return min([abs(mOCST-90) for mOCST in mOCST]) + 90 
    else: return 0

def find_Z1_pair(lep_pts,lep_etas,lep_phis,lep_charges,lep_types):
    # mass of Z boson candidate 1
    mOCST = []
    for i in range(len(lep_pts)-1):
        for j in range(len(lep_pts)):
            if j>i:
                if lep_charges[i]!=lep_charges[j] and lep_types[i]==lep_types[j]:
                    mll = 2*lep_pts[i]*lep_pts[j]
                    cosh = math.cosh(lep_etas[i]-lep_etas[j])
                    cos = math.cos(lep_phis[i]-lep_phis[j])
                    mll *= ( cosh - cos )
                    mOCST.append(math.sqrt(mll)/1000.)      
                else: mOCST.append(0)
    mOCST_minus90 = [abs(mOCST-90) for mOCST in mOCST]
    if mOCST_minus90.index(min(mOCST_minus90))==0: return '01'
    elif mOCST_minus90.index(min(mOCST_minus90))==1: return '02'
    elif mOCST_minus90.index(min(mOCST_minus90))==2: return '03'
    elif mOCST_minus90.index(min(mOCST_minus90))==3: return '12'
    elif mOCST_minus90.index(min(mOCST_minus90))==4: return '13'
    elif mOCST_minus90.index(min(mOCST_minus90))==5: return '23'

def calc_m_Z2(lep_pts,lep_etas,lep_phis,lep_charges,lep_types,Z1_pair):
    # mass of Z boson candidate 2
    mOCST = []
    for i in range(len(lep_pts)-1):
        for j in range(len(lep_pts)):
            if j>i and str(j) not in Z1_pair and str(i) not in Z1_pair and lep_charges[i]!=lep_charges[j] and lep_types[i]==lep_types[j]:
                mll = 2*lep_pts[i]*lep_pts[j]
                cosh = math.cosh(lep_etas[i]-lep_etas[j])
                cos = math.cos(lep_phis[i]-lep_phis[j])
                mll *= ( cosh - cos )
                mOCST.append(math.sqrt(mll)/1000.)      
    if len(mOCST)>0: return mOCST[0]
    else: return 0

def calc_mllll(lep_pts,lep_etas,lep_phis):
    theta_0 = 2*math.atan(math.exp(-lep_etas[0]))
    theta_1 = 2*math.atan(math.exp(-lep_etas[1]))
    theta_2 = 2*math.atan(math.exp(-lep_etas[2]))
    theta_3 = 2*math.atan(math.exp(-lep_etas[3]))
    p_0 = lep_pts[0]/math.sin(theta_0)
    p_1 = lep_pts[1]/math.sin(theta_1)
    p_2 = lep_pts[2]/math.sin(theta_2)
    p_3 = lep_pts[3]/math.sin(theta_3)
    pz_0 = p_0*math.cos(theta_0)
    pz_1 = p_1*math.cos(theta_1)
    pz_2 = p_2*math.cos(theta_2)
    pz_3 = p_3*math.cos(theta_3)
    px_0 = p_0*math.sin(theta_0)*math.cos(lep_phis[0])
    px_1 = p_1*math.sin(theta_1)*math.cos(lep_phis[1])
    px_2 = p_2*math.sin(theta_2)*math.cos(lep_phis[2])
    px_3 = p_3*math.sin(theta_3)*math.cos(lep_phis[3])
    py_0 = p_0*math.sin(theta_0)*math.sin(lep_phis[0])
    py_1 = p_1*math.sin(theta_1)*math.sin(lep_phis[1])
    py_2 = p_2*math.sin(theta_2)*math.sin(lep_phis[2])
    py_3 = p_3*math.sin(theta_3)*math.sin(lep_phis[3])
    sumpz = pz_0 + pz_1 + pz_2 + pz_3
    sumpx = px_0 + px_1 + px_2 + px_3
    sumpy = py_0 + py_1 + py_2 + py_3
    sumE = p_0 + p_1 + p_2 + p_3
    mllll = sumE**2 - sumpz**2 - sumpx**2 - sumpy**2
    return math.sqrt(mllll)/1000.

def read_file(path,sample):
    start = time.time()
    print("\tProcessing: "+sample+" file")
    data_all = pd.DataFrame()
    mc = uproot.open(path)["mini"]
    numevents = uproot.numentries(path, "mini")
    for data in mc.iterate(["lep_pt","lep_eta","lep_phi","lep_type","lep_charge",
                         "mcWeight","scaleFactor_PILEUP","scaleFactor_ELE","scaleFactor_MUON", # add more variables here if you make cuts on them ,  
                            "scaleFactor_LepTRIGGER"], flatten=False, entrysteps=2500000, outputtype=pd.DataFrame, entrystop=numevents*fraction):

        nIn = len(data.index)

        # label for each mc type
        data['type'] = np.vectorize(mc_type)(sample)

        # label for channel (ee, mm or em)
        data['Channel'] = np.vectorize(channel)(data.lep_type)

        data['SumCharges'] = np.vectorize(sum_charge)(data.lep_charge)
        data['MinmOCST'] = np.vectorize(calc_min_mOCST)(data.lep_pt,data.lep_eta,data.lep_phi,data.lep_charge,data.lep_type)
        data['LepPt0'] = np.vectorize(lep_pt_0)(data.lep_pt)
        data['LepPt1'] = np.vectorize(lep_pt_1)(data.lep_pt)
        data['LepPt2'] = np.vectorize(lep_pt_2)(data.lep_pt)
        data['MZ1'] = np.vectorize(calc_m_Z1)(data.lep_pt,data.lep_eta,data.lep_phi,data.lep_charge,data.lep_type)
        data['Z1_pair'] = np.vectorize(find_Z1_pair)(data.lep_pt,data.lep_eta,data.lep_phi,data.lep_charge,data.lep_type)
        data['MZ2'] = np.vectorize(calc_m_Z2)(data.lep_pt,data.lep_eta,data.lep_phi,data.lep_charge,data.lep_type,data.Z1_pair)
        data['Mllll'] = np.vectorize(calc_mllll)(data.lep_pt,data.lep_eta,data.lep_phi)

        if 'data' not in sample:
            data['weight'] = np.vectorize(calc_weight)(data.mcWeight,data.scaleFactor_PILEUP,data.scaleFactor_ELE,data.scaleFactor_MUON,data.scaleFactor_LepTRIGGER)
            data['weight'] = np.vectorize(get_xsec_weight)(data.weight,sample)
        else:
            data['weight'] = 1

        data.drop(["Z1_pair","lep_pt","lep_eta","lep_phi","lep_type","lep_charge","mcWeight","scaleFactor_PILEUP","scaleFactor_ELE","scaleFactor_MUON","scaleFactor_LepTRIGGER"], axis=1, inplace=True)

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
