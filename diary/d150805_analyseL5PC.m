function [BI_values,OT_values,BC_values,figh] =d150805_analyseL5PC(injection_level,inj_curr,background_curr,options)
    
    LEN_DEND_SEG = 50;
    LEN_EPSILON = 1;
fname = sprintf('tree_L5PC.neu');
L5PCtree = load_tree(fname);
L5PCtree.Gm = 1./30000.0;
L5PCtree.Ri = 100.0;

%fname = sprintf('tree_SHstellate.neu');
%SHtree = load_tree(fname);
%SHtree.Gm = 1./30000.0;
%SHtree.Ri = 139.09;


%a = sse_tree(tree,[]);
%figure;
%dA_tree(L5PCtree);
%figure;
%dA_tree(SHtree);


%dstats_tree(stats_tree({L5PCtree,SHtree},{'L5PC','SHStellate'},[],'-x -s'))

tree = L5PCtree;

 % find terminating branches
    T = T_tree (tree);
    term_nodes = find(T);
    % work out distances to root
    dists = Pvec_tree(tree);
    % work out paths
    ipar = ipar_tree (tree);

    % BI: Branch injection
    BI_name  = 'L5PCtemplate[]';
    % find sections that are units lambda from the soma
    BI_dists = find(dist_tree(tree,[150]));
    BI_node = find(tree.R==find(ismember(tree.rnames,BI_name)));
    BI_node = intersect(BI_dists,BI_node);
    
    % MB: main branch
    MB_name = 'L5PCtemplate[]';
    MB_node = find(tree.R==find(ismember(tree.rnames,MB_name)));
    MB_node = intersect(term_nodes,MB_node);
    
    % soma / root
    % note that the connecting node for soma is always the even number
    % (also note that this usually is 2, but not guaranteed)
    soma_node = find(tree.R==find(ismember(tree.rnames,'L5PCtemplate[]')));
    soma_node = soma_node(end);
    
    % Set of Vm for injecting at BI
    input = background_curr*ones(size (dists, 1), 1);
    input(BI_node) = inj_curr+background_curr;
    sse = sse_tree(tree,input);
    
    % 1. Obtain SSE response for injecting current at BT, for points along
    % path
    sect = [soma_node MB_node];
    % following line: taken from TREES toolbox, plotsect_tree.m
    indy = ipar  (sect (1, 2), 1 : find (ipar (sect (1, 2), :) == sect (1, 1)));
    indy = fliplr(indy);

    BI_values.name = MB_name;
    BI_values.indy = indy;
    BI_values.vm = sse(indy);
    BI_values.dists = dists(indy);
    
   
    
    % 3. SSE for injection at BT for points in same domain but at other
    % termination points to root.
    % Note that we don't find BC for a split at the first level 
    if nc > 1,
        % Generate lists for other terminating branches
        for i = 1:1:nl-1,
            
            BC_name = 'L5PCtemplate[]';
            BC_node = find(tree.R==find(ismember(tree.rnames,BC_name)));
            % find the section that is the terminating section
            BC_node = intersect(BC_node,term_nodes);
            
            sect = [soma_node BC_node];
            
            % following line: taken from TREES toolbox, plotsect_tree.m
            indy = ipar  (sect (1, 2), 1 : find (ipar (sect (1, 2), :) == sect (1, 1)));
            indy = fliplr(indy);
            
            % save Vm and distances 
            BC_values(i).name = BC_name;
            BC_values(i).indy = indy;
            BC_values(i).vm = sse(indy);
            BC_values(i).dists = dists(indy);
            
            
        end
    end
    
    if strfind (options, '-s'),
        % Plot the figure
        % Plot OT, if OT
        figh = figure;
        hold on;
        if nb > 1,
            plot(OT_values.dists,OT_values.vm,'LineWidth',1);
        end
        % Plot BC
        if nc > 1,
            for i = 1:size(BC_values,2),
                plot(BC_values(i).dists,BC_values(i).vm,'LineWidth',2);
            end
        end
        % plot BI
        plot(BI_values.dists,BI_values.vm,'LineWidth',3);
        
        
    end