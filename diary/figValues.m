% Print values for figures
function figValues(figh)
    
    figure(figh);
    
    % Defaults 
    width = 1;     % Width in inches
    height = 3;    % Height in inches
    alw = 3;    % AxesLineWidth
    fsz = 11;      % Fontsize
    lw = 1.5;      % LineWidth
    msz = 8;       % MarkerSize

    set(figh, 'Color', 'w');
    % The properties we've been using in the figures
    set(figh,'defaultLineLineWidth',lw);   % set the default line width to lw
    set(figh,'defaultLineMarkerSize',msz); % set the default line marker size to msz
    set(figh,'defaultLineLineWidth',lw);   % set the default line width to lw
    set(figh,'defaultLineMarkerSize',msz); % set the default line marker size to msz
    set(figh,'defaultLineMarkerSize',msz); % set the default line marker size to msz
    
    a=findobj(figh); % get the handles associated with the current figure
    
    allaxes=findall(a,'Type','axes');
    %alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    
    set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',3,'FontSize',20);
    %set(alllines,'Linewidth',2);
    set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',20);

end