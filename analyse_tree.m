function analyse_tree(filename,tree_index) 

    global trees

    % set a few values that weren't set when we imported the tree from NEURON
    
    %Rm = 150;
    % input resistance
    %trees{tree_index}.Ri = Rm/(4*pi*trees{tree_index}.D^2); % Ri = Rm/4.pi.a^2 (a in cm^2)
    trees{tree_index}.Ri = 25e7; % Ri = Rm/4.pi.a^2 (a in cm^2)
    %trees{tree_index}.Gm = 13e-6; % units: S. Taken from "Princ Neural Sci", p137
    trees{tree_index}.Gm = 2000; % units: S. Taken from "Princ Neural Sci", p137
    trees{tree_index}.Gi = 100;
    %trees{tree_index}.Cm = 1;
    trees{tree_index}.D(:)=trees{tree_index}.D(:)*100;
    
    
    % values from Gidon & Segev 2009
    
    
    
    return

    % compensate for different units of 
    %diam_factor = 80;
    %trees{1}.D = trees{1}.D*diam_factor;

    figure(1);
    sse = sse_tree(tree_index);
    imagesc(sse); axis equal; axis off; colorbar;
    saveas(1,sprintf('%s_fig1.png',filename),'png')
    
    % Inject at branches
    T = T_tree (tree_index);
    term = find(T,1);
    figure(2);
    sse = sse_tree(tree_index,term,'-s');
    saveas(2,sprintf('%s_fig2.png',filename),'png')

    %trees{1}.D = trees{1}.D*10
    figure(3);
    plot_tree(tree_index,[1 0 0]);
    saveas(3,sprintf('%s_fig3.png',filename),'png')

    figure(4);
    elen_tree(tree_index,'-s');
    saveas(4,sprintf('%s_fig4.png',filename),'png')

    

    %imagesc(sse); axis equal; axis off; colorbar;
    
    



