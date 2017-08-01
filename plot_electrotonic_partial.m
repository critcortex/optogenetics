% PLOT_ELECTROTONIC
%
% Finds distal location to inject current in
%
function [BI_values,OT_values,BC_values,figh] = plot_electrotonic_partial(tree,nb,nc,nl,inj_curr,background_curr,partial_branches,options)
    


    if (nargin<5)||isempty(inj_curr),
        inj_curr = 1;
    end

    if (nargin<6)||isempty(background_curr),
        background_curr = 0;
    end
    
    if (nargin<7)||isempty(partial_branches),
        partial_branches = 1:nc;
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
    BI_name  = strcat('dend0', repmat('_0',1,nl-1));
    BI_node = find(tree.R==find(ismember(tree.rnames,BI_name)));
    BI_node = intersect(term_nodes,BI_node);
    
    % soma / root
    % note that the connecting node for soma is always the even number
    % (also note that this usually is 2, but not guaranteed)
    soma_node = find(tree.R==find(ismember(tree.rnames,'soma')));
    soma_node_end = soma_node(end);
    
    % Set of Vm for injecting at BI
    input = zeros(size (dists, 1), 1);
    for i = 1:size(partial_branches,2),
        refString = sprintf('dend%u',partial_branches(i)-1); %as we start with dend0 rather than dend1
        index = strfind(tree.rnames,refString);
        index = find(~cellfun(@isempty,index));
        index = ismember(tree.R,index);
        input(index) = 1;
    end
    % and set the soma to receive photocurrents too
    % input(soma_node) = 1;
    % now set the photocurrent to be it's chosen value
    input = background_curr*input;
    % add extra current for injection site
    input(BI_node) = inj_curr+background_curr;
    % and finally, calculate the value in response
    sse = sse_tree(tree,input);
    
    % 1. Obtain SSE response for injecting current at BT, for points along
    % path
    sect = [soma_node_end BI_node];
    % following line: taken from TREES toolbox, plotsect_tree.m
    indy = ipar  (sect (1, 2), 1 : find (ipar (sect (1, 2), :) == sect (1, 1)));
    indy = fliplr(indy);

    BI_values.name = BI_name;
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
        
        sect = [soma_node_end OT_node];
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
            
            sect = [soma_node_end BC_node];
            
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
        
        % find max length, so that we can normalize 
        maxlength = max(BI_values.dists);
        figh = figure('position', [0, 0, 600, 350]) ;
        hold on;
        color = [0.75 0.75 0.75];
        if nb > 1,
            plot(OT_values.dists/maxlength,OT_values.vm,'LineWidth',3,'Color',color);
        end
        % Plot BC
        color = [0.45 0.45 0.45];
        if nc > 1,
            for i = 1:size(BC_values,2),
                plot(BC_values(i).dists/maxlength,BC_values(i).vm,'LineWidth',3,'Color',color);
            end
        end
        % plot BI
        plot(BI_values.dists/maxlength,BI_values.vm,'LineWidth',6,'Color','k');
        
        set(gca, 'FontSize', 18);
        xlim([0 1.1]);
        xlabel('Distance (normalized)');
        ylabel('V_{m} (mV)');
        set(gca,'XTick',0:0.25:1);   
        %set(gcf,'units','centimeters','position',[0,0,50,30])
        
    end
    
    
    
    
    
        
       



    