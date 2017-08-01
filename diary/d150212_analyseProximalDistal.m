
function d150212_analyseProximalDistal()
    clear all
    %dataProx = load('d150209_proximalInjection_test.mat');
    %processData(dataProx.data,'d150209_proximalInjection_test');
    
    clear dataClean
    
    dataDist = load('d150129_distalInjection_124_3.mat');
    processData(dataDist.data,'d150129_distalInjection');
    
end

function processData(data,savename)
    
    for i = 1:size(data,2)
        dataClean(i,:) = [data{i}.name data{i}.diff_evoked_intrins data{i}.diff_max_min_intrinsic data{i}.vm_at_soma_inj];   
    end

    % Plot evoked - intrinsic
    %plotForTrees(dataClean(:,1:4),'Evoked - intrinsic',strcat(savename,'_evokedIntrinsic'));

    % Plot max-min intrinsic (i.e. photocurrent only)  
    plotForTrees([dataClean(:,1:3) dataClean(:,5)],'Max. photocurrent (normalized)',strcat(savename,'_photoBW'));
    
    % Plot max-min evoked (i.e. point where it is single distal current injection)  
    plotForTrees([dataClean(:,1:3) dataClean(:,6)],'Evoked bandwidth',strcat(savename,'_evokedBW'));
    
end

function plotForTrees(data,clabel,savename)
    
    plotPolarityBranch([data(:,1:2) data(:,4)],clabel,savename);
    plotBranchLevels(data(:,2:4),clabel,savename);
    plotPolarityLevels([data(:,1) data(:,3:4)],clabel,savename);
end

function plotPolarityBranch(data,clabel,savename)
    plotXY(data,'Polarity','Branching',clabel,strcat(savename,'_polarityBranch'));
end

function plotBranchLevels(data,clabel,savename)

    plotXY(data,'Branching','#Levels',clabel,strcat(savename,'_branchLevels'));
end

function plotPolarityLevels(data,clabel,savename)

    plotXY(data,'Polarity','#Levels',clabel,strcat(savename,'_polarityLevels'));
end

function plotXY(data,xlbl,ylbl,clabel,savename)

    xx = log(data(:,1));
    yy = log(data(:,2));
    ww = data(:,3);
    % normalize ww
    ww = ww/max(ww);
    
    incr = 0.01;
    gx = [0:incr:max(xx)]; %6.22];
    gy = [0.05:incr:max(yy)];

    xxmod = xx*((1/incr)-1);
    yymod = yy*((1/incr)-1);
    [XX,YY] = meshgrid(gx,gy);
    WI = griddata(xx,yy,ww,XX,YY,'cubic');
    
    clims = [min(ww) max(ww)];
    subplot(1,3,3);
    fig = plotScanImage(WI,[],xxmod,yymod,clims,xlbl,ylbl,clabel);
    figure(fig);
 
    print(strcat(savename,'.png'),'-dpng','-r300');
    saveas(fig,strcat(savename,'.png'),'png');
    %print(strcat(savename,'.svg'),'-dpng','-r300');
    saveas(fig,strcat(savename,'.eps'),'epsc');
    %plot2svg(strcat(savename,'.svg'),fig);
    close(fig);
    
    close('all');
end

% 
% function plotXY(data,xlbl,ylbl,clabel,savename)
% % Data is in format of [x,y,c/s]
%     ms = 100;
%     fig = figure();
%     set(gca, 'FontSize', 18);
%     scatter(data(:,1),data(:,2),ms,data(:,3),'filled');
%     xlabel(xlbl);
%     ylabel(ylbl);
%     set(gca,'xscale','log')
%     set(gca,'yscale','log')
%     h = colorbar;
%     ylabel(h,clabel);
%     title(savename);
%     
%     %saveas(fig,strcat(savename,'.png'),'png');
%     print(strcat(savename,'.png'),'-dpng','-r300');
%     close(fig)
% end


function plotXY_2Dinterpolate(data,xlbl,ylbl,clabel,savename)
    ms = 100;
    fig = figure();
    scatter(data(:,1),data(:,2),ms,'o');
    hold on
    
    x = data(:,1);
    y = data(:,2);
    v = data(:,3);
    
    % define grid
    [xq,yq] = meshgrid(0:1:1000, 0:1:1000);
    vq = griddata(x,y,v,xq,yq);
    mesh(xq,yq,vq);
    
    
    xlabel(xlbl);
    ylabel(ylbl);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    set(gca, 'FontSize', 180);
    h = colorbar;
    ylabel(h,clabel);
    title(savename);
    %saveas(fig,strcat(savename,'.png'),'png');
    print(strcat(savename,'.png'),'-dpng','-r300');
    close(fig)

end
