function f = plotScanImage(data,lines,xxmod,yymod,clims,xlab,ylab,clab)
cLow = 0;
cHigh = 1;
f= figure;
pos = [0.15 0.2 0.6 0.7];
axes1 = axes('Parent',f,'FontName','Arial','FontSize',12,...
    'YScale','log','YMinorTick','on',...
    'XScale','log','XMinorTick','on',...
    'TickDir','out','TickLength',[0.04 0.075],...
    'Position', pos,...
    'LineWidth',1.5,'clim',clims);
%xlim(axes1,[0 500]);
%ylim(axes1,[0 1000]);
hold(axes1,'all');
incr = 0.005;
pos(1) = pos(1)+incr;
pos(2) = pos(2)+incr;
pos(3) = pos(3)-incr;
pos(4) = pos(4)-incr;
caxis([0 1]);
axes2 = axes('Parent',f,'Layer','top','position',pos,'clim',clims);
set(axes1, 'Color', [1 0.5 1])
imagescnan(data,'Parent',axes2,'NanColor',[0.99999 0.999999 0.999999]);
set(gca, 'CLim', [cLow, cHigh]);
colormap(CMRmap(20));
%imagescnan(ZI,'Parent',axes2)
set(axes2,'Visible','off');
axis xy;
%title('Mean dwell times for reproducible stimulus','Fontsize',12)
xlabel(axes1,xlab,'Fontsize',12,'FontName','Arial');
ylabel(axes1,ylab,'Fontsize',12,'FontName','Arial');
set(gca, 'CLim', [cLow, cHigh]);
c =colorbar('peer',axes2,'LineWidth',2,'fontsize',12,'fontname','Arial','TickDir','out');
cpos = get(c,'position');
cpos(1) = pos(1)+pos(3)+0.03;
set(c,'position',cpos,'Clim',clims);
ylabel(c,clab);
xlim([-1 500]);
hold on

plot(xxmod,yymod,'o','MarkerSize',6,'Linewidth',2,'MarkerEdgeColor',[0.7 0.7 0.7]);
set(gca, 'Clim', clims)

end