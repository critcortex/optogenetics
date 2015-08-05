
fname = sprintf('tree_L5PC.neu');
L5PCtree = load_tree(fname);
L5PCtree.Gm = 1./30000.0;
L5PCtree.Ri = 100.0;

fname = sprintf('tree_SHstellate.neu');
SHtree = load_tree(fname);
SHtree.Gm = 1./30000.0;
SHtree.Ri = 139.09;


%a = sse_tree(tree,[]);
%figure;
%dA_tree(L5PCtree);
%figure;
%dA_tree(SHtree);


dstats_tree(stats_tree({L5PCtree,SHtree},{'L5PC','SHStellate'},[],'-x -s'))