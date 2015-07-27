
#import os, sys
from subprocess import call
"""

Sets up directory structures
Runs nrnivmodl to generate objects from mod files
"""

# create folders
folders = ['experiments/','diary/']
for f in folders:
    try:
        call(['mkdir',f])
        print 'Created directory', f
    except:
        pass
    
    
# run nrnivmodl
#call(['nrnivmodl'])
call(['nrnivmodl nrn/mod/'])
# will create executable folder


# EXPERIMENT COLLECTION:
"""
sjarvis1@bg-sjarvis1.bg.ic.ac.uk:/home/sjarvis/git/exp_data_optogen/experiments/
"""