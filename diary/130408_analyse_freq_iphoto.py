import run_analysis

af = run_analysis.AnalyticFrame()
basename = '130328_FIcurve'
exp_comp_list = [['_freq%g_pw10*NpHR_whole*ChR_whole','ChR whole,NpHR whole'],
                 ['_freq%g_pw10*NpHR_soma*ChR_whole','ChR whole,NpHR soma'],
                 ['_freq%g_pw10*NpHR_apical*ChR_whole','ChR whole,NpHR apical'],
                 ['_freq%g_pw10*NpHR_none*ChR_whole','ChR whole,NpHR none'],
                 ['_freq%g_pw10*NpHR_whole*ChR_apical','ChR apical,NpHR whole'],
                 ['_freq%g_pw10*NpHR_apical*ChR_apical','ChR apical,NpHR apical'],
                 ['_freq%g_pw10*NpHR_soma*ChR_apical','ChR apical,NpHR soma'],
                 ['_freq%g_pw10*NpHR_whole*ChR_apical','ChR apical,NpHR whole']]
expss = [ec[0] for ec in exp_comp_list]
explabels = [ec[1] for ec in exp_comp_list]
af.populate_expset(basename,expss,explabels,run_analysis.DEFAULT_FREQ)
af.run_analysis_menu()