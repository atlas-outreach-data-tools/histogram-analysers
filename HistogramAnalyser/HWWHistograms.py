Channel = {
    # change plotting parameters
    'bin_width':1,
    'num_bins':3,
    'xrange_min':-0.5,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.1 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

Mll = {
    # change plotting parameters
    'bin_width':5,
    'num_bins':22,
    'xrange_min':0,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.1 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

NJets = {
    # change plotting parameters
    'bin_width':1,
    'num_bins':11,
    'xrange_min':-0.5,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.4 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

BTags = {
    # change plotting parameters
    'bin_width':1,
    'num_bins':2,
    'xrange_min':-0.5,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.1 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

SumLepPt = {
    # change plotting parameters
    'bin_width':10,
    'num_bins':21,
    'xrange_min':0,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.1 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

MET = {
    # change plotting parameters
    'bin_width':10,
    'num_bins':21,
    'xrange_min':0,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.4 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

METLLDeltaPhi = {
    # change plotting parameters
    'bin_width':0.2,
    'num_bins':32,
    'xrange_min':-3.2,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.4 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

LepDeltaPhi = {
    # change plotting parameters
    'bin_width':0.2,
    'num_bins':32,
    'xrange_min':-3.2,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.4 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

TransMass = {
    # change plotting parameters
    'bin_width':10,
    'num_bins':31,
    'xrange_min':0,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.4 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

hist_dict = {"Mll":Mll,"NJets":NJets,"BTags":BTags,"SumLepPt":SumLepPt,"MET":MET,"METLLDeltaPhi":METLLDeltaPhi,"LepDeltaPhi":LepDeltaPhi,"TransMass":TransMass}#,"Channel":Channel}
