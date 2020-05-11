Mllll = {
    # change plotting parameters
    'bin_width':5,
    'num_bins':34,
    'xrange_min':80,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.4 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

MZ1 = {
    # change plotting parameters
    'bin_width':5,
    'num_bins':24,
    'xrange_min':0,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.1 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

MZ2 = {
    # change plotting parameters
    'bin_width':5,
    'num_bins':24,
    'xrange_min':0,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.1 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

LepPt3 = {
    # change plotting parameters
    'bin_width':20,
    'num_bins':20,
    'xrange_min':0,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.4 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

LepPt2 = {
    # change plotting parameters
    'bin_width':20,
    'num_bins':20,
    'xrange_min':0,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.4 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

LepPt1 = {
    # change plotting parameters                                                                                                   
    'bin_width':20,
    'num_bins':20,
    'xrange_min':0,
    'log_y':False,

    # change aesthetic parameters if you want                                                                                      
    'y_label_x_position':-0.09, # 0.09 to the left of y axis                                                                       
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0                        
    'linear_top_margin':1.4 # to decrease the separation between data and the top of the figure, pick a number closer to 1         
}

LepPt0 = {
    # change plotting parameters                                                                                                   
    'bin_width':20,
    'num_bins':20,
    'xrange_min':0,
    'log_y':False,

    # change aesthetic parameters if you want                                                                                      
    'y_label_x_position':-0.09, # 0.09 to the left of y axis                                                                       
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0                        
    'linear_top_margin':1.4 # to decrease the separation between data and the top of the figure, pick a number closer to 1         
}

MinmOCST = {
    # change plotting parameters
    'bin_width':5,
    'num_bins':24,
    'xrange_min':0,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.1 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

SumCharges = {
    # change plotting parameters
    'bin_width':1,
    'num_bins':9,
    'xrange_min':-4.5,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.1 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

Channel = {
    # change plotting parameters
    'bin_width':1,
    'num_bins':5,
    'xrange_min':-0.5,
    'log_y':False,

    # change aesthetic parameters if you want
    'y_label_x_position':-0.09, # 0.09 to the left of y axis
    'legend_loc':'best',
    'log_top_margin':10000, # to decrease the separation between data and the top of the figure, remove a 0
    'linear_top_margin':1.1 # to decrease the separation between data and the top of the figure, pick a number closer to 1
}

hist_dict = {"Channel":Channel,"SumCharges":SumCharges,"MinmOCST":MinmOCST,"LepPt0":LepPt0,"LepPt1":LepPt1,"LepPt2":LepPt2,"LepPt3":LepPt3,"MZ1":MZ1,"MZ2":MZ2,'Mllll':Mllll}
