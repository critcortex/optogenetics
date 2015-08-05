clear all
close all

EXPNAME = '_d150209_analyseTreesProximalXX';

start_trees

INJ_AMOUNT = 100;

PLOT_INDIV = 0;
PLOT_SUMMARY = 0;

trees_config = [ [1 1 124] ;[1 2 7]; [2 1 62]; ... 
          [2 2 6]; [2 7 3]; [2 61 2]; ...
          [4 5 3]; [4 2 5]; [4 30 2]; [4 1 31];  [11 1 11]; [11 10 2]; ...
          [18 1 7]; [18 2 3]; 
          [18 6 2] ; [31 1 4]; [31 3 2]; [62 1 2]; 
          [124 1 1]];

trees_config = [ [2 2 6];] % [2 7 3];[18 1 7]; [18 6 2]];

plot_show = '';
if PLOT_INDIV,
    plot_show = '-s';
end

global trees
global resampled_trees

data = {};

for i = 1:size(trees_config,1),
    close all;
    i
    t = trees_config(i,:);
    fname = sprintf('../experiments/electrotonic/tree_nb%u_nc%u_nl%u.neu',t(1),t(2),t(3));
    fname
    tt = load_tree(fname);
    tt = resample_tree(tt,1,'-d');
    trees{i} = tt;
    analyse_tree(fname,i);
    
    
    
    %[xd t2] = xdend_tree(trees{i});
    %resampled_trees{i} = t2;
    %xdend_tree(t2,'-s');
    %saveas(gcf,sprintf('%s_shape.png',fname),'png');
    %close all;
    
    
    [a,b,c,figh] = plot_electrotonic_proximal(trees{i},t(1),t(2),t(3),1,0,plot_show);
    if PLOT_INDIV,
        title(gca,sprintf('Inj=%g for I_{photo} = 0 in tree (%g,%g,%g) - proximal',INJ_AMOUNT,t(1),t(2),t(3)));
        saveas(figh,sprintf('%s%s_electrotonic_prox.png',fname,EXPNAME),'png');
    end
    
    if PLOT_SUMMARY,
        fig_bg = figure();
        fig_inj = figure();
    end
    
    for curr = 0.01:0.01:0.05,
        
        [a1,b1,c1,figh1] = plot_electrotonic_proximal(trees{i},t(1),t(2),t(3),INJ_AMOUNT,-1*curr,plot_show);
        [a2,b2,c2,figh2] = plot_electrotonic_proximal(trees{i},t(1),t(2),t(3),INJ_AMOUNT,curr,plot_show);
        [a3,b3,c3,figh3] = plot_electrotonic_proximal(trees{i},t(1),t(2),t(3),0,-1*curr,plot_show);
        [a4,b4,c4,figh4] = plot_electrotonic_proximal(trees{i},t(1),t(2),t(3),0,curr,plot_show);

        
        if curr == 0.05,
           data{i}.name = t;
           data{i}.vm_max_evoked = a1;
           data{i}.vm_min_evoked = a2;
           data{i}.vm_max_intrinsic = a3;
           data{i}.vm_min_intrinsic = a4;
        end
        
        
        if PLOT_INDIV,
            
            figure(figh1);
            title(gca,sprintf('Inj=%g for I_{photo} = %g in tree (%g,%g,%g) - proximal',INJ_AMOUNT,-1*curr,t(1),t(2),t(3)));
            saveas(figh1,sprintf('%s%s_inj%g_electrotonic_i-%g_prox.png',fname,EXPNAME,INJ_AMOUNT,curr),'png');
        
            figure(figh2);
            title(gca,sprintf('Inj=%g for I_{photo} = %g in tree (%g,%g,%g) - proximal',INJ_AMOUNT,curr,t(1),t(2),t(3)));
            saveas(figh2,sprintf('%s%s_inj%g_electrotonic_i%g_prox.png',fname,EXPNAME,INJ_AMOUNT,curr),'png');
            
            figure(figh3);
            title(gca,sprintf('Inj=0 for I_{photo} = %g in tree (%g,%g,%g) - proximal',-1*curr,t(1),t(2),t(3)));
            saveas(figh3,sprintf('%s%s_noInj_electrotonic_i-%g_prox.png',fname,EXPNAME,curr),'png');

            figure(figh4);
            title(gca,sprintf('Inj=0 for I_{photo} = %g in tree (%g,%g,%g) - proximal',curr,t(1),t(2),t(3)));
            saveas(figh4,sprintf('%s%s_noInj_electrotonic_i%g_prox.png',fname,EXPNAME,curr),'png');

            close(figh1);
            close(figh2);
            close(figh3);
            close(figh4);
        end
        
        if PLOT_SUMMARY,
            line_factor = 100;
            figure(fig_bg);
            hold on;
            plot(a1.dists,a1.vm,'LineWidth',curr*line_factor);
            plot(a2.dists,a2.vm,'LineWidth',curr*line_factor,'Color','r');
            if t(1) > 1,
                plot(b1.dists,b1.vm,'LineWidth',curr*line_factor,'LineStyle','--');
                plot(b2.dists,b2.vm,'LineWidth',curr*line_factor,'Color','r','LineStyle','--');
            end
            hold off
            ylim([-12 12]);
            title(gca,sprintf('Evoked response for tree (%g,%g,%g) - proximal',t(1),t(2),t(3)));
            
            figure(fig_inj);
            hold on;
            plot(a3.dists,a3.vm,'LineWidth',curr*line_factor);
            plot(a4.dists,a4.vm,'LineWidth',curr*line_factor,'Color','r');
            hold off
            ylim([-12 12]);
            title(gca,sprintf('Background response for tree (%g,%g,%g) - proximal',t(1),t(2),t(3)));
            
        end
    end
    
    if PLOT_SUMMARY,
        
        saveas(fig_bg,sprintf('%s%s_electrotonicRangeWithCurrent.png',fname,EXPNAME),'png')
        saveas(fig_inj,sprintf('%s%s_electrotonicRange.png',fname,EXPNAME),'png')
    end
    close all;
    
    
    
    
    close all;
     
end


close all




%% Spread of gain
% Max - min optogenetics
% and check that it's the same for evoked as for intrinsic
for i = 1:size(trees_config,1),
    
    dd_max_evoked = data{i}.vm_max_evoked.dists;
    vm_max_evoked = data{i}.vm_max_evoked.vm;
    vm_at_soma_evoked = vm_max_evoked(dd_max_evoked==0);
    dd_min_evoked = data{i}.vm_min_evoked.dists;
    vm_min_evoked = data{i}.vm_min_evoked.vm;
    vm_min_at_soma_evoked = vm_min_evoked(dd_min_evoked==0);
    
    dd_max_intrinsic = data{i}.vm_max_intrinsic.dists;
    vm_max_intrinsic = data{i}.vm_max_intrinsic.vm;
    vm_at_soma_intrinsic = vm_max_intrinsic(dd_max_intrinsic==0);
    dd_min_intrinsic = data{i}.vm_min_intrinsic.dists;
    vm_min_intrinsic = data{i}.vm_min_intrinsic.vm;
    vm_min_at_soma_intrinsic = vm_min_intrinsic(dd_min_intrinsic==0);
    
        
    trees_config(i,:);
    data{i}.diff_max_min_evoked = abs(vm_at_soma_evoked - vm_min_at_soma_evoked);
    data{i}.diff_max_min_intrinsic = abs(vm_at_soma_intrinsic - vm_min_at_soma_intrinsic);
    
    
    data{i}.diff_max_min_evoked;
    data{i}.diff_max_min_intrinsic;
end





%% Compare shift of comparative gains (at soma)
% THis is performed by looking at effective change of evoked - intrinsic

for i = 1:size(trees_config,1),
    
    dd_max_evoked = data{i}.vm_max_evoked.dists;
    vm_max_evoked = data{i}.vm_max_evoked.vm;
    vm_at_soma_evoked = vm_max_evoked(dd_max_evoked==0);
    
    dd_max_intrinsic = data{i}.vm_max_intrinsic.dists;
    vm_max_intrinsic = data{i}.vm_max_intrinsic.vm;
    vm_at_soma_intrinsic = vm_max_intrinsic(dd_max_intrinsic==0);
        
    trees_config(i,:);
    data{i}.diff_evoked_intrins = vm_at_soma_evoked - vm_at_soma_intrinsic;
    data{i}.diff_evoked_intrins;
    
    
end


%save('d150209_proximalInjection.mat','-struct',data);
save('d150209_proximalInjection.mat','data');






