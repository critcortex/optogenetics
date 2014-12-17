import matplotlib
import run_analysis
import numpy as np

abstract = {'filebase': '',}

expL5PC_vivo = {'filebase': '140723_L5PC_basal_irr_ChRonly',
                'light_on': 1050,
                'light_dur': 1000,
                'tstop': 2100,
                'freqs': range(2,15,2)+range(15,151,5),
                'Js':[2.],
                'nsite_range': [80],
                'factors': [0.001,0.125,0.25,0.5,0.75,1.] ,
                'irrs': [0.002]}


expL5PC_vitro = {'filebase': '130710_increasing_uneq_irradiance_factor_somacurr',
                'light_on': 700,
                'light_dur': 1000,
                'tstop': 2500,
                'current_amps': np.arange(-2.,2.01,0.1),
                'factors': [0.125,0.25,0.5,1.,2.,4.,8.] ,
                'irrs': [1.0]}

expStellate = {'filebase': '140408_scan_SHstell',
               'light_on': 700,
                'light_dur': 1000,
                'tstop': 2700,
                'current_amps': np.arange(-1.,3.1,0.2),
                'freqs': range(10,201,10),
                'Js':[2.],
                'nsite_range': [40],
                'stimlocs': ['stim2','stim3'],
                'factors': [ 0.0625, 0.125,0.25,0.375,0.5,0.75,1.,1.5,2.],
                'irrs': [0.007,0.012]}


def plot_Fig1():
    """ Abstract cells """
    pass









def plot_L5PC_vivo():
    """ 
    L5PC cells:
        in vivo
        modulation index
    """
    data = expL5PC_vivo
    description = ['whole','whole']
    
    for irr in data['irrs']:
        af = run_analysis.AnalyticFrame(cmap_index=1)
    
        af.update_params({'tstart':data['light_on'],'tstop':data['light_on']+data['light_dur'],
                              'tstart_bg': 50,'tstop_bg':data['light_on'],
                              'tstart_post':data['light_on']+data['light_dur'],'tstop_post':data['tstop']})
        
        exp_comp_list = [['irr%.3f_'%irr+'factor%.2f_'%f+'freq%g'+'_J%g_nsites%g'%(data['Js'][0],data['nsite_range'][0])+'_NpHR_%s_ChR_%s'%(description[1],description[0]),'%.2f'%f] for f in data['factors']]
    
        print exp_comp_list
  
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        print explabels
        
        af.populate_expset(data['filebase'],expss,explabels, [data['freqs']])
        
        af.submenu_extractSpikes()
        af.submenu_runFI()
        for exp in af.experimentset:
            exp.calculate_responses('FI')
            exp.calculate_responses('FI_bg')
            exp.calculate_responses('FI_post')
        
        af.submenu_save()
        af.submenu_print()
        
        af.perform_analysis()
        mi = af._calculateModulationIndex('MI_bg')
        #mi2 = af._calculateModulationIndex('MI_bg',5)
        
        
        af.submenu_plot(10, 'Fig2_L5PC_vivo_'+data['filebase']+'FI_gain_varyFactor_invivo_irr%g'%irr+'_NpHR_%s_ChR_%s'%(description[1],description[0]))
        
        print 'Modulation index for L5PC_vivo = ',mi
        

def plot_L5PC_vitro():
    """ 
    L5PC cells:
        in vitro
        modulation index
    """
    data = expL5PC_vitro
    description = ['whole','whole']
    
    for irr in data['irrs']:
        af = run_analysis.AnalyticFrame(cmap_index=1)
    
        af.update_params({'tstart':data['light_on'],'tstop':data['light_on']+data['light_dur']})
            
        exp_comp_list = [['_irr%.1f'%(irr)+'_factor%.2f'%(factor)+'_Idist%.1f_'+'NpHR_%s_ChR_%s'%(description[1],description[0]),'%g'%factor] for factor in data['factors']]
        
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        
        af.populate_expset(data['filebase'],expss,explabels, [data['current_amps']])
        
        af.submenu_load()
        af.submenu_print()
        
        af.submenu_plot(5, 'Fig2_L5PC_vitro_'+data['filebase']+'FI_gain_irr%.1f_varyFactor_exp%s%s_'%(irr,description[0],description[1]))
        
        af.perform_analysis()
        mi = af._calculateModulationIndex()
        print 'Modulation index for L5PC_vitro = ',mi
        """
        mi = af._calculateModulationIndex(eval_at_x=0.)
        print 'Modulation index for L5PC_vitro = ',mi
        """
        
        
           
def plot_Stellate_vivo():
    
    data = expStellate
    description = ['whole','whole']
    
    for irr in data['irrs']:
        af = run_analysis.AnalyticFrame(cmap_index=2)
    
        af.update_params({'tstart':data['light_on'],'tstop':data['light_on']+data['light_dur'],
                              'tstart_bg': 0,'tstop_bg':data['light_on'],
                              'tstart_post':data['light_on']+data['light_dur'],'tstop_post':data['tstop']})
        
        exp_comp_list = [['irr%.3f_'%irr+'factor%.2f_'%f+'freq%g'+'_J%g'%(2.)+'_NpHR_%s_ChR_%s'%(description[1],description[0]),'=%.3f'%f] for f in data['factors']]
    
        print exp_comp_list
  
        expss = [ec[0] for ec in exp_comp_list]
        explabels = [ec[1] for ec in exp_comp_list]
        print explabels
        
        af.populate_expset(data['filebase'],expss,explabels, [data['freqs']])
                
        af.submenu_extractSpikes()
        
        af.submenu_runFI()
        for exp in af.experimentset:
            exp.calculate_responses('FI')
            exp.calculate_responses('FI_bg')
            exp.calculate_responses('FI_post')
    
        af.submenu_save()
        af.submenu_print()
        af.submenu_plot(10, 'Fig3_stellate_vivo_'+data['filebase']+'FI_bggain_invivo_varyFactor_irr%.3f'%irr)
        
    
def plot_Stellate_vitro():
    data = expStellate
    description = ['whole','whole']
    
    for irr in data['irrs']:
        for stimloc in data['stimlocs']:
            af = run_analysis.AnalyticFrame(cmap_index=2)
        
            af.update_params({'tstart':data['light_on'],'tstop':data['light_on']+data['light_dur'],
                                  'tstart_bg': 0,'tstop_bg':data['light_on'],
                                  'tstart_post':data['light_on']+data['light_dur'],'tstop_post':data['tstop']})
            
            
            exp_comp_list = [['irr%.3f_'%irr+'factor%.2f_'%f+'I%.2f'+'_stimloc_%s'%(stimloc)+'_NpHR_%s_ChR_%s'%(description[1],description[0]),'=%.3f'%f] for f in data['factors']]
            print exp_comp_list
      
            expss = [ec[0] for ec in exp_comp_list]
            explabels = [ec[1] for ec in exp_comp_list]
            print explabels
            
            af.populate_expset(data['filebase'],expss,explabels, [data['current_amps']])
            af.submenu_extractSpikes()
        
            af.submenu_runFI()
            for exp in af.experimentset:
                exp.calculate_responses('FI')
                exp.calculate_responses('FI_bg')
                exp.calculate_responses('FI_post')
       
            af.submenu_save()
        
            af.submenu_print()
            af.submenu_plot(0, 'Fig3_stellate_vitro_'+data['filebase']+'FI_gain_invitro_varyFactor_irr%.3f_stimloc%s_TEST'%(irr,stimloc))
            #af.submenu_plot(5, 'Fig3_stellate_vitro_'+data['filebase']+'FI_gain_invitro_varyFactor_irr%.3f_stimloc%s'%(irr,stimloc))
            
        
def plot_Stellate():
    #plot_Stellate_vivo()
    plot_Stellate_vitro()
    
#plot_L5PC_vivo()    
#plot_L5PC_vitro()    
plot_Stellate()    
    