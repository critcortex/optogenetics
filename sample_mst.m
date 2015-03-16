start_trees;

%% First: create tree and save to hoc file
num_pts = 10;

X = rand(num_pts,1)*200 - 100;
Y = rand(num_pts,1)*200 - 100;
Z = zeros(num_pts,1);

msttree = MST_tree(1,[0;X],[0;Y],[0;Z]);
msttree = resample_tree(msttree,5);
msttree = soma_tree(msttree);

figure;
plot_tree(msttree,[1 0 0])

neuron_tree(msttree,'sample.hoc');


%% Second: load back into TREES

load_tree('sample_out.neu')

% compare
figure;
plot_tree(1,[1 0 0])




