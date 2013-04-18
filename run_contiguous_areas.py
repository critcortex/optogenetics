# -*- coding: utf-8 -*-

import run_stimulation

experiment_type = 'opsinonly' #'BAP'
experiment_name = '130204_sepNpHR_noisyn'
experiment_name = '130204_testChR_noisyn'
experiment_name = "130204_bothconcurrent_stimdur500"

# for expressions, array is [expression, pd, del_light, t_on, t_off, distance*]
# *NB that distance is only relevant for apical_distal or apical_proximal
sampleexpression = [5e-4,10,600,200,350,600]
samplelocations = ['basal','soma','axon','apical_proximal','apical_distal']
samplevalues = [[x,sampleexpression] for x in samplelocations]


ChRexpression= [5e-4,10,400,400,350,600]
ChRlocations = [None] + samplelocations[:2]
NpHRexpression = [5e-4,10,600,400,350,600]
NpHRlocations = samplelocations + [None]

ChR_expression = [5e-4,10,500,1000,350,600]
NpHR_expression = [5e-4,10,1500,2000,350,600]
#both_expression = [5e-4,10,2500,3000,350,600]
both_expression = [5e-4,10,500,500,350,600]

        
for i in range(1,len(NpHRlocations)+1): # loop for NpHR expression

    for j in range(1,len(ChRlocations)+1):

        ChRvalues= [] #[[x,ChR_expression] for x in ChRlocations] 
        ChRvalues = ChRvalues + [[x,both_expression] for x in ChRlocations[:j]] 
        if ChRlocations[j-1] is None:
            ChR_descript = 'noChR'
        else:
            ChR_descript = "ChR_"+ChRlocations[j-1]+'_downwards'
        
        NpHRvalues = [] #[[x,NpHR_expression] for x in NpHRlocations[-i:]]
        NpHRvalues = NpHRvalues + [[x,both_expression] for x in NpHRlocations[-i:]]

        if NpHRlocations[-i] is None:
            NpHR_descript = "noNpHR_"
        else:
            NpHR_descript = "NpHR_"+NpHRlocations[-i]+"upwards_"
        
        
        
        opdict = {  'ChR': ChRvalues , 
                    'NpHR': NpHRvalues} 
        print opdict

        expdict = {}
        expdict['experiment_type']=experiment_type
        expdict['expname']=experiment_name+"_%g_"%(i-1)+NpHR_descript+ChR_descript
        expdict['opdict']=opdict
        expdict['savedata']=True
        expdict['tstart'] = 0
        expdict['tstop'] = 1500
        print expdict['expname']

        run_stimulation.main(expdict)


