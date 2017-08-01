clear all
close all

EXPNAME = '_d150129_analyseTrees';

% load style
cssname = 'default_paper';


start_trees

PLOT_INDIV = 0;
PLOT_SUMMARY = 1;

ORANGE = [1 0.7 0.02];

INJ_AMOUNT = 1;

trees_config = [ [1 1 124] ;[1 2 7]; [1 123 2]; [2 2 6]; [2 7 3]; [2 1 62]; [2 61 2]; ...
          [4 5 3]; [4 2 5]; [4 30 2]; [4 1 31];  [11 1 11]; [11 10 2]; ...
          [18 1 7]; [18 2 3]; [18 6 2] ; [31 1 4]; [31 3 2]; [62 1 2]; [124 1 1]];

%trees_config = [[2 61 2];];% [1 2 7]; [2 2 6];  [2 7 3]; [4 5 3] ]; %[2 7 3];[18 1 7]; [18 6 2]];
% [(2,2,6),(2,7,3),(4,5,3)]
%trees_config = [ [1 1 124]; [11 1 11]; [2 2 6]; [2 61 2]]; 
%trees_config = [[1 123 2]];


% for 124+-1segments
trees_config = [ [1 1 124];  [1 123 2]; [2 1 62]; [2 61 2]; [3 1 41]; [3 40 2]; ...
    [4 1 31]; [4 2 5]; [4 5 3]; [4 30 2]; [5 1 25]; [5 24 2]; [25 1 5]; [25 4 2]; ...
    [31 1 4]; [31 3 2]; [41 1 3]; [41 2 2]; [62 1 2]; [124 1 1] ];


% for 124+-3segments
trees_config = [[1 1 124]; [1 2 7]; [1 3 5]; [1 123 2]; [2 1 62]; [2 2 6]; ...
    [2 61 2]; [3 1 41]; [3 1 42]; [3 40 2]; [3 41 2]; [4 1 31]; [4 2 5]; ... 
    [4 5 3]; [4 30 2]; [5 1 25]; [5 24 2]; [6 1 21]; [6 4 3]; [6 20 2]; ...
    [7 1 18]; [7 17 2]; [9 1 14]; [9 13 2]; [11 1 11]; [11 10 2]; ...
    [14 1 9]; [14 8 2]; [18 1 7]; [18 2 3]; [18 6 2]; [21 1 6]; [21 5 2]; ...
    [25 1 5]; [25 4 2]; [31 1 4]; [31 3 2]; [41 1 3]; [41 2 2]; [42 1 3]; ...
    [42 2 2]; [62 1 2]];

trees_config = [[1 1 124]; [1 2 7]; [2 2 6] ; [4 1 31];];

plot_show = '';
if PLOT_INDIV,
    plot_show = '-s';
end

global trees
global resampled_trees

data = {};



for i = 1:size(trees_config,1),
    close all;
    
    t = trees_config(i,:);
    fname = sprintf('../experiments/electrotonic/tree_nb%u_nc%u_nl%u.neu',t(1),t(2),t(3));
    tt = load_tree(fname);
    trees{i} = tt;
    resample_tree(i,1,'-d');
    analyse_tree(fname,i);
      
    
%     [a,b,c,figh] = plot_electrotonic(trees{i},t(1),t(2),t(3),INJ_AMOUNT,0,plot_show);
%     if PLOT_INDIV,
%         title(gca,sprintf('Inj=1 for I_{photo} = 0 in tree (%g,%g,%g)',t(1),t(2),t(3)));
%         %saveas(figh,sprintf('%s%s_electrotonic.png',fname,EXPNAME),'png');
%         print(sprintf('%s%s_electrotonic.png',fname,EXPNAME),'-dpng','-r300');
%     end
    
    if PLOT_SUMMARY,
        fig_bg = figure('position', [0, 0, 600, 350]) ;
        fig_inj = figure('position', [0, 0, 600, 350]) ;
    end
    
    % photocurrent values to be tested
    for curr =0.01:0.01:0.05, % 0.005:0.005, % 
        
%         % when we have interaction between injected (evoked) amount and
%         % background (photocurrent) for photoactive and photoinhibition
        [a1,b1,c1,figh1] = plot_electrotonic(trees{i},t(1),t(2),t(3),INJ_AMOUNT,-1*curr,plot_show);
        [a2,b2,c2,figh2] = plot_electrotonic(trees{i},t(1),t(2),t(3),INJ_AMOUNT,curr,plot_show);
        % when we have intrinsic background (photocurrent) for photoactive and photoinhibition
        [a3,b3,c3,figh3] = plot_electrotonic(trees{i},t(1),t(2),t(3),0,-1*curr,plot_show);
        [a4,b4,c4,figh4] = plot_electrotonic(trees{i},t(1),t(2),t(3),0,curr,plot_show);
        % when only a single distal injection
        %[a5,b5,c5,figh5] = plot_electrotonic(trees{i},t(1),t(2),t(3),INJ_AMOUNT,0,plot_show);
        
                
%         if curr == 0.005,
%            data{i}.name = t;
%            data{i}.vm_max_evoked = a1;
%            data{i}.vm_min_evoked = a2;
%            data{i}.vm_max_intrinsic = a3;
%            data{i}.vm_min_intrinsic = a4;
%            data{i}.vm_inj_evoked = a5;
%         end
        
        
        if PLOT_INDIV,
            
%             figure(figh1);
%             title(gca,sprintf('Inj=%g for I_{photo} = %g in tree (%g,%g,%g)',INJ_AMOUNT,-1*curr,t(1),t(2),t(3)));
%             %saveas(figh1,sprintf('%s%s_inj%g_electrotonic_i-%g.png',fname,EXPNAME,INJ_AMOUNT,curr),'png');
%             print(sprintf('%s%s_inj%g_electrotonic_i-%g.png',fname,EXPNAME,INJ_AMOUNT,curr),'-dpng','-r300');
%         
%             figure(figh2);
%             title(gca,sprintf('Inj=%g for I_{photo} = %g in tree (%g,%g,%g)',INJ_AMOUNT,curr,t(1),t(2),t(3)));
%             %saveas(figh2,sprintf('%s%s_inj%g_electrotonic_i%g.png',fname,EXPNAME,INJ_AMOUNT,curr),'png');
%             print(sprintf('%s%s_inj%g_electrotonic_i%g.png',fname,EXPNAME,INJ_AMOUNT,curr),'-dpng','-r300');
%             
%             figure(figh3);
%             title(gca,sprintf('Inj=0 for I_{photo} = %g in tree (%g,%g,%g)',-1*curr,t(1),t(2),t(3)));
%             %saveas(figh3,sprintf('%s%s_noInj_electrotonic_i-%g.png',fname,EXPNAME,curr),'png');
%             print(sprintf('%s%s_noInj_electrotonic_i-%g.png',fname,EXPNAME,curr),'-dpng','-r300');
%             
%             figure(figh4);
%             title(gca,sprintf('Inj=0 for I_{photo} = %g in tree (%g,%g,%g)',curr,t(1),t(2),t(3)));
%             %saveas(figh4,sprintf('%s%s_noInj_electrotonic_i%g.png',fname,EXPNAME,curr),'png');
%             print(sprintf('%s%s_noInj_electrotonic_i%g.png',fname,EXPNAME,curr),'-dpng','-r300');
            
%             figure(figh5);
%             figValues(figh5);
%             name = sprintf('%s%s_Inj_noElectrotonic_i%g',fname,EXPNAME,curr);
%             %title(gca,sprintf('Injected current I=1 in tree (%g,%g,%g)',t(1),t(2),t(3)));
%             saveas(figh5,sprintf('%s_2.png',name),'png');
%             export_fig(sprintf('%s.png',name));
%             export_fig(sprintf('%s.eps',name));
  
%             close(figh1);
%             close(figh2);
%             close(figh3);
%             close(figh4);
%             close(figh5);
        end

        
        if PLOT_SUMMARY,
            line_factor = 100;
            curr*line_factor
            
            figure(fig_inj);
            hold on;
            plot(a1.dists/max(a1.dists),a1.vm,'LineWidth',curr*line_factor);
            plot(a2.dists/max(a2.dists),a2.vm,'LineWidth',curr*line_factor,'Color','r');
%            if t(1) > 1,
%                plot(b1.dists,b1.vm,'LineWidth',curr*line_factor,'LineStyle','--');
%                plot(b2.dists,b2.vm,'LineWidth',curr*line_factor,'Color','r','LineStyle','--');
%            end
            hold off
            
            figure(fig_bg);
            hold on;
            plot(a3.dists/max(a3.dists),a3.vm,'LineWidth',curr*line_factor,'Color','b');
            plot(a4.dists/max(a4.dists),a4.vm,'LineWidth',curr*line_factor,'Color',ORANGE);
            hold off
            
        end
    end
    
    if PLOT_SUMMARY,
        
        figValues(fig_inj);
        xlabel('Distance (normalized)');
        ylabel('I_{photo} (a.u.)');
        ylim([-1000 1000]);
        title(gca,sprintf('Evoked response for tree (%g,%g,%g)',t(1),t(2),t(3)));
        
        name = sprintf('%s%s_electrotonicRangeWithCurrent.png',fname,EXPNAME);
        saveas(fig_inj,sprintf('%s_2.png',name),'png')
        export_fig(sprintf('%s.png',name));
        export_fig(sprintf('%s.eps',name));
        
        
        
        
        figValues(fig_bg);
        xlabel('Distance (normalized)');
        ylabel('I_{photo} (a.u.)');
        ylim([-1000 1000]);
        title(gca,sprintf('Background response for tree (%g,%g,%g)',t(1),t(2),t(3)));
            
        name = sprintf('%s%s_electrotonicRange.png',fname,EXPNAME);
        saveas(fig_bg,sprintf('%s_2.png',name),'png')
        export_fig(sprintf('%s.png',name));
        export_fig(sprintf('%s.eps',name));
    end
    
    
    close all;
     
end


close all




%% Spread of gain
% Calculate the max - min (bandwidth) at the soma for effect of photocurrents
% and check that it's the same for evoked as for intrinsic i.e. there's no
% interaction between injected current and photocurrents
% NOTE THAT THIS CONSIDERATION IS IMPORTANT <-- PASSIVE ELEMENTS
% (what happens when we have non-linear elements???)


% Compare shift of comparative gains (at soma)
% THis is performed by looking at effective change of evoked - intrinsic

for i = 1:size(trees_config,1),
    
    trees_config(i,:)
    
    % injection + excitatory photocurrent
    dd_max_evoked = data{i}.vm_max_evoked.dists;
    vm_max_evoked = data{i}.vm_max_evoked.vm;
    vm_at_soma_evoked = vm_max_evoked(dd_max_evoked==0);
    data{i}.vm_at_soma_evoked = vm_at_soma_evoked;
    
    % injection + inhibitory photocurrent
    dd_min_evoked = data{i}.vm_min_evoked.dists;
    vm_min_evoked = data{i}.vm_min_evoked.vm;
    vm_min_at_soma_evoked = vm_min_evoked(dd_min_evoked==0);
    
    % excitatory photocurrent
    dd_max_intrinsic = data{i}.vm_max_intrinsic.dists;
    vm_max_intrinsic = data{i}.vm_max_intrinsic.vm;
    vm_at_soma_intrinsic = vm_max_intrinsic(dd_max_intrinsic==0);
    % inhibitory photocurrent
    dd_min_intrinsic = data{i}.vm_min_intrinsic.dists;
    vm_min_intrinsic = data{i}.vm_min_intrinsic.vm;
    vm_min_at_soma_intrinsic = vm_min_intrinsic(dd_min_intrinsic==0);
    
    % check that the values for using inh and exc photocurrents aren't affected by current injection    
    data{i}.diff_max_min_evoked = abs(vm_at_soma_evoked - vm_min_at_soma_evoked);
    data{i}.diff_max_min_intrinsic = abs(vm_at_soma_intrinsic - vm_min_at_soma_intrinsic);
    
        
    % should give the evoked (current injection) response only
    data{i}.diff_evoked_intrins = vm_at_soma_evoked - vm_at_soma_intrinsic;
    data{i}.diff_evoked_intrins;
    
    % and finally, get shift of evoked at soma
    dd_inj_evoked = data{i}.vm_inj_evoked.dists;
    vm_inj_evoked = data{i}.vm_inj_evoked.vm;
    vm_at_soma_inj = vm_inj_evoked(dd_inj_evoked==0);
    data{i}.vm_at_soma_inj  = vm_at_soma_inj;
    
    % note that this should be the same as the data{i}.diff_evoked_intrins
    
    
end


%save('d150129_distalInjection_124_test.mat','data');
























% 
% % Plot all on same axis
% fig = figure();
% hold on;
% for i = 1:size(trees_config,1),
%     i
%     %t = trees_config(i,:);
%     %fname = sprintf('../tree_nb%u_nc%u_nl%u.neu',t(1),t(2),t(3));
%     %tt = load_tree(fname);
%     %trees{i} = tt;
%     %analyse_tree(fname,i);
%     %resample_tree(i,10,'-d');
%     
%     if t(1) == 1,
%         col = 'k';
%     elseif t(2) == 1,
%         col = 'b';
%     else
%         col = 'r';
%     end
%     
%     [a,b,c,figh] = plot_electrotonic(trees{i},t(1),t(2),t(3),1,0);
%     
%     plot(a.dists/max(a.dists),a.vm,'LineWidth',4,'Color',col);
%     
%     
%     for curr = 0.01:0.01:0.05,
%         [a1,b1,c1,figh2] = plot_electrotonic(trees{i},t(1),t(2),t(3),1,-1*curr);
%         [a2,b2,c2,figh3] = plot_electrotonic(trees{i},t(1),t(2),t(3),1,curr);
%         plot(a1.dists/max(a.dists),a1.vm,'LineWidth',curr*100,'Color',col);
%         plot(a2.dists/max(a.dists),a2.vm,'LineWidth',curr*100,'Color',col);
%         
%     end
%     
%     
%     
% end
% saveas(fig,sprintf('%s_electrotonicRangeCompare.png',fname),'png');
