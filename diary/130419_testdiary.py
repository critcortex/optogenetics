import run_experiments as rs

DEFAULT_FREQ = [5]

es = rs.ExpSetup()
es.run_test_exp('130417_test',params={'num_threads':2})

