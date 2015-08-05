clear all
close all

EXPNAME = '_d150209_analyseTreesProximalXX';

start_trees

INJ_AMOUNT = 100;

PLOT_INDIV = 0;
PLOT_SUMMARY = 0;

trees_config = [ [1 2 7]; [2 1 62]; ... 
          [2 2 6]; [2 7 3]; [2 61 2]; ...
          [4 5 3]; [4 2 5]; [4 30 2]; [4 1 31];  [11 1 11]; [11 10 2]; ...
          [18 1 7]; [18 2 3]; 
          [18 6 2] ; [31 1 4]; [31 3 2]; [62 1 2]; 
          [124 1 1];[1 1 124] ;];

trees_config = [[1 2 7]]; % [2 2 6];]; % [18 1 7]; [1 2 7]]; % [2 7 3];[18 1 7]; [18 6 2]];

plot_show = '';
if PLOT_INDIV,
    plot_show = '-s';
end

global trees
global resampled_trees

data = {};

photo_max = 0.0005;

for i = 1:size(trees_config,1),
    close all;
    %i
    t = trees_config(i,:);
    fname = sprintf('../experiments/electrotonic/tree_nb%u_nc%u_nl%u.neu',t(1),t(2),t(3));
    %fname
    tt = load_tree(fname);
    tt = resample_tree(tt,1,'-d');
    trees{i} = tt;
    analyse_tree(fname,i);
    
    for injSite = 0:t(3)-1, 
        [a1,b1,c1,figh1] = plot_electrotonic_scan(trees{i},t(1),t(2),t(3),injSite,1.,0,'-s');
        drawnow; pause(0.05);
        % get gca / gcf and rescale so that max=1
        
        print(sprintf('%s_injSite%d.png',fname,injSite),'-dpng','-r300');
        close all
        
        [a1,b1,c1,figh1] = plot_electrotonic_scan(trees{i},t(1),t(2),t(3),injSite,1.,photo_max,'-s');
        drawnow; pause(0.05);
        % get gca / gcf and rescale so that max=1
        print(sprintf('%s_injSite%d_ChR2.png',fname,injSite),'-dpng','-r300');
        close all
        
        [a1,b1,c1,figh1] = plot_electrotonic_scan(trees{i},t(1),t(2),t(3),injSite,1.,-1*photo_max,'-s');
        drawnow; pause(0.05);
        % get gca / gcf and rescale so that max=1
        print(sprintf('%s_injSite%d_NpHR.png',fname,injSite),'-dpng','-r300');
        close all
        
    end
end
        



        