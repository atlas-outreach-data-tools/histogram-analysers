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

    'data': {
        'list' : ['data_A','data_B','data_C','data_D']
    },

    'HZZ' : {
        'list' : ['ZH125_ZZ4lep','WH125_ZZ4lep','VBFH125_ZZ4lep','ggH125_ZZ4lep'],
        'color' : "#CCCCCC"
    },

    'ZZ' : {
        'list' : ['llll'],
        'color' : "#2A9000"
    },

    'ttbar' : {
        'list' : ['ttbar_lep'],
        'color' : "#1C5EA8"
    },

    'Z' : {
        'list' : ['Zee','Zmumu'],
        'color' : "#F66A00"
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
    return round(weight,5)


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
    return round(lep_pt[0]/1000,2)

def lep_pt_1(lep_pt):
    return round(lep_pt[1]/1000,2)

def lep_pt_2(lep_pt):
    return round(lep_pt[2]/1000,2)

def lep_pt_3(lep_pt):
    return round(lep_pt[3]/1000,2)

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
    if len(mOCST)>0: return round(min(mOCST),2)
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
    if len(mOCST)>0: return round(min([abs(mOCST-90) for mOCST in mOCST]) + 90,2)
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
    if len(mOCST)>0: return round(mOCST[0],2)
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
    return round(math.sqrt(mllll)/1000,2)

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
        data['LepPt3'] = np.vectorize(lep_pt_3)(data.lep_pt)
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

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator,LogLocator,LogFormatterSciNotation # for minor ticks
import HZZHistograms
import labelfile
stack_order = ['Z','ttbar','ZZ']

def plot_data(data):

    signal_format = 'hist' # 'line' for line above SM stack
                           # 'hist' for bar above SM stack
                           # None for signal as part of SM stack
    Total_SM_label = False # for Total SM black line in plot and legend
    plot_label = r'$H \rightarrow ZZ^* \rightarrow \ell\ell\ell\ell$'
    signal_label = r'Signal ($m_H=125$ GeV)' # r''

    # *******************
    # general definitions (shouldn't need to change)
    lumi_used = str(lumi*fraction)    
    signal = None
    for s in samples.keys():
        if s not in stack_order and s!='data': signal = s

    for x_variable,hist in HZZHistograms.hist_dict.items():

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
    
        plt.savefig("HZZ_"+x_variable+".pdf")
    
    return signal_x,mc_x_tot

signal_yields,background_yields = plot_data(data)
