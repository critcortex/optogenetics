from neuron import h, run, init

# this file is from the TREES toolbox, so make sure you have the right value for PATH_TO_FILE  
h.load_file("PATH_TO_FILE/neu_tree.hoc")


# create tree in h. session here
# TODO

# then when that's done, save it as a .neu file ready for the TREES toolbox
h.neu_tree('mytree.neu')



# and then change over to matlab and load it in using 
"""
fname = sprintf('mytree.neu');
tt = load_tree(fname);
xdend_tree(tt,'-s');
saveas(gcf,sprintf('%s_shape.png',fname),'png');
"""