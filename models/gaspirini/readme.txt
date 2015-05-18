NEURON mod files from the paper:
Sonia Gasparini, Michele Migliore, and Jeffrey C. Magee
On the initiation and propagation of dendritic spikes in CA1 pyramidal neurons,
J.Neurosci. 24:11046-11056 (2004).

The paper investigates the dendritic spike probability as a function 
of the number of activated synaptic input and their spatial and temporal distribution.
The results suggest that highly synchronized and moderately clustered inputs 
have the highest probability to generate a dendritic spike.

The main effect is demonstrated here by reproducing the simulations
in Fig.5A of the paper, in which 5 synapses were activated 
in the same compartment with different time delays.
Simulations can be run for different kind of synaptic input:
- current pulses 
- AMPA 
- AMPA+NMDA

Under unix systems:
to compile the mod files use the command 
nrnivmodl 
and run the simulation hoc file with the command 
nrngui forfig5A.hoc

Under Windows systems:
to compile the mod files use the "mknrndll" command.
A double click on the simulation file
forfig5A.hoc 
will open the simulation window.

Questions on how to use this model
should be directed to michele.migliore@pa.ibf.cnr.it


