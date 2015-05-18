
function d150212_analyseProximalDistal()
    clear all
    dataProx = load('../d150209_proximalInjection.mat');
    processData(dataProx.data,'d150209_proximalInjection');
    
    clear dataClean
    
    dataDist = load('../d150129_distalInjection.mat');
    processData(dataDist.data,'d150129_distalInjection');
    
end

function processData(data,savename)
    
    for i = 1:size(data,2)
        dataClean(i,:) = [data{i}.name data{i}.diff_evoked_intrins data{i}.diff_max_min_intrinsic];   
    end

    % Plot evoked - intrinsic
    plotForTrees(dataClean(:,1:4),'Evoked - intrinsic',strcat(savename,'_evokedIntrinsic'));

    % Plot max-min    
    plotForTrees([dataClean(:,1:3) dataClean(:,5)],'Gain bandwidth',strcat(savename,'_gainBW'));
    
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
% Data is in format of [x,y,c/s]
    ms = 100;
    fig = figure();
    scatter(data(:,1),data(:,2),ms,data(:,3),'filled');
    xlabel(xlbl);
    ylabel(ylbl);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    h = colorbar;
    ylabel(h,clabel);
    title(savename);
    saveas(fig,strcat(savename,'.png'),'png');
    close(fig)
end

