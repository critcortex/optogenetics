import run_experiments as rs

DEFAULT_FREQ = [5]

es = rs.ExpSetup()
es.run_test_exp('sampletest',params={'num_threads':1},runon='local')

