% PLOT_ELECTROTONIC
%
% Finds proximal location to inject current in
% - inject current at level = injection_level
% - 
function [BI_values,OT_values,BC_values,figh] = plot_electrotonic_scan(tree,nb,nc,nl,injection_level,inj_curr,background_curr,options)
    
    LEN_DEND_SEG = 50;
    LEN_EPSILON = 1;

    if (nargin<5)||isempty(injection_level),
        injection_level = nl-1;
    end
    
    if (nargin<6)||isempty(inj_curr),
        inj_curr = 1;
    end

    if (nargin<7)||isempty(background_curr),
        background_curr = 0;
    end
    
    if (nargin<8)||isempty(options),
        options = '';
    end

    figh = 0;
    BI_values = {};
    OT_values = {};
    BC_values = {};
    
    tree.Ri = 150;
    tree.Gm = 5e-4;
    tree.Cm = 1;
    
    
    % find terminating branches
    T = T_tree (tree);
    term_nodes = find(T);
    % work out distances to root
    dists = Pvec_tree(tree);
    % work out paths
    ipar = ipar_tree (tree);

    % BI: Branch injection
    BI_name  = strcat('dend0', repmat('_0',1,injection_level));%'dend0';
    % find sections that are units lambda from the soma
    %BI_dists = find(dist_tree(tree,[LEN_DEND_SEG-LEN_EPSILON]));
    %BI_node = find(tree.R==find(ismember(tree.rnames,BI_name)));
    %BI_node = intersect(BI_dists,BI_node);
    % new version
    BI_node = find(tree.R==find(ismember(tree.rnames,BI_name)));
    % grab the middle one from this set of nodes (as there should be more 
    % than one, as we resampled the tree)
    BI_node=BI_node(ceil(end/2));
    
    % MB: main branch
    MB_name = strcat('dend0', repmat('_0',1,nl-1));
    MB_node = find(tree.R==find(ismember(tree.rnames,MB_name)));
    MB_node = intersect(term_nodes,MB_node);
    
    % soma / root
    % note that the connecting node for soma is always the even number
    % (also note that this usually is 2, but not guaranteed)
    soma_node = find(tree.R==find(ismember(tree.rnames,'soma')));
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
    
    
    % 2. SSE for injection at BT, for points along OT i.e. Other tree
    if nb > 1,
        % Generate the name for the other branches where current is not
        % injected
        OT_name = strcat('dend1',repmat('_0',1,nl-1));
        OT_node = find(tree.R==find(ismember(tree.rnames,OT_name)));
        % find the section that is the terminating section
        OT_node = intersect(OT_node,term_nodes);
        
        sect = [soma_node OT_node];
        % following line: taken from TREES toolbox, plotsect_tree.m
        indy = ipar  (sect (1, 2), 1 : find (ipar (sect (1, 2), :) == sect (1, 1)));
        indy = fliplr(indy);
        
        OT_values.name = OT_name;
        OT_values.indy = indy;
        OT_values.vm = sse(indy);
        OT_values.dists = dists(indy);
    end
    
    % 3. SSE for injection at BT for points in same domain but at other
    % termination points to root.
    % Note that we don't find BC for a split at the first level 
    if nc > 1,
        % Generate lists for other terminating branches
        for i = 1:1:nl-1,
            BC_name = strcat(repmat('_0',1,i),'_1',repmat('_0',1,nl-1-i));
            BC_name = strcat('dend', BC_name(2:end));
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
        % Plot OT (other tree) in grey
        color = [0.45 0.45 0.45];
        if nb > 1,
            plot(-1*OT_values.dists,OT_values.vm,'LineWidth',2,'Color',color);
        end
        % Plot BC (branch cousins?) in black
        if nc > 1,
            for i = 1:size(BC_values,2),
                plot(BC_values(i).dists,BC_values(i).vm,'LineWidth',3,'Color','k');
            end
        end
        % plot BI in red
        color = [0.66 0.12 0];
        plot(BI_values.dists,BI_values.vm,'LineWidth',4,'Color',color);
        
        
    end
    
    
    
    
    
        
       



    