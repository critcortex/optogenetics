%%
% Plotting function from ESN Matlab project (2009)
% xx = xvalues
% yy = yvalues
% ww = observations at x and y coordinates

%xx = [xx bbones];
%yy = [yy reservoir./bbones];
%ww = [ww; lhat'];
%f = plotFunction(xx,yy,[],ww);
function g = plotFunction(xx,yy,ww)
%origxx = xx;
%origyy = yy;

xx=log(xx);
yy =log(yy);

incr = 0.01;
gx = [0:incr:6.22];
gy = [0.05:incr:6.91];

xxmod = xx*((1/incr)-1);
yymod = yy*((1/incr)-1);

spoint = find(xxmod==min(xxmod));
for i = 1:size(spoint,2)
    lines(i,1) = spoint(i);
    if(i ==size(spoint,2))
        lines(i,2) = size(xxmod,2);
    else
        lines(i,2) = spoint(i+1)-1;
    end
end

[XX,YY] = meshgrid(gx,gy);

WI = griddata(xx,yy,ww,XX,YY);
%clims = [0 20];
%clims = [-0.25 -0.06];
subplot(1,3,3);
g = plotScanImage(WI,lines,xxmod,yymod,clims,...
    'Number of clusters (log)','Size of clusters (log)','\lambda');
figure(g);
[C,h] = contour(WI,[0 0],'color',[0.6 0 0],'LineWidth',2);

max(ww)
min(ww)


end

