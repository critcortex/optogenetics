
#import os, sys
from subprocess import call
"""

Sets up directory structures
Runs nrnivmodl to generate objects from mod files
"""

# create folders
folders = ['experiments/']
for f in folders:
    try:
        call(['mkdir',f])
        print 'Created directory', f
    except:
        pass
    
    
# run nrnivmodl
call(['nrnivmodl'])
# will create executable folder