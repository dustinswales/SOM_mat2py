function plot_SOMcomposites(varIN,lonINf,latINf,pat_freq,minP,maxP,dP,lonmin,latmin,lonmax,latmax,cmap,num_rows,num_cols,cbarTitle)
                 
% Scale colormap to data range provided
x1 = minP:(maxP-minP)/(size(cmap,1)-1):maxP;
x2 = minP:dP:maxP;
cmap2 = zeros(length(x2),3);
cmap2(:,1) = interp1(x1,cmap(:,1),x2);
cmap2(:,2) = interp1(x1,cmap(:,2),x2);
cmap2(:,3) = interp1(x1,cmap(:,3),x2);
cmap = cmap2;

% Rows and columns for plot title
if (num_rows == 2 && num_cols == 2)
    cols = [1,2,1,2];
    rows = [1,1,2,2];
end
if (num_rows == 3 && num_cols == 3)
    cols = [1,2,3,1,2,3,1,2,3];
    rows = [1,1,1,2,2,2,3,3,3];
end
if (num_rows == 4 && num_cols == 4)
    cols = [1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4];
    rows = [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4];
end
% Load in coastal ans state boundary data (for plotting).
coast  = load('coast');
states = shaperead('usastatelo.shp','usegeocoords',true);

figure('color','w','units','inches','position',[1,1,7,8],'PaperPositionMode','Auto')
colormap(cmap)
for iplot=1:num_rows*num_cols
    % Set up plot
    subaxis(num_rows,num_cols,iplot,'SpacingVert',0.0,'SpacingHoriz',0.0,...
        'PaddingBottom',0.03,'PaddingTop',0.01,'PaddingLeft',0.01,'PaddingRight',0.01,...
        'MarginBottom', 0.1, 'MarginTop',0.05, 'MarginLeft',0.05, 'MarginRight',0.0)
    
    axesm('eqdcylin','MapLatLimit',[latmin latmax],'MapLonLimit',[lonmin lonmax],...
        'grid','on','PLineLocation',[0:5:90],'MLineLocation',[-180:5:180],...
        'MLineFill',50,'PlineFill',50,'glinewidth',0.3)
    axis off; framem on;
    setm(gca,'MLabelLocation',40)
    plotm(coast.lat,coast.long,'color','k')
    
    % Contour data
    contourfm(latINf,lonINf,varIN(:,:,iplot),[minP:dP:maxP],'LineStyle','none');
    caxis([minP maxP]);
    %Ptitle = ['(',num2str(rows(iplot),'%.i'),',',num2str(cols(iplot),'%.i'),')'];        
    %t = title(Ptitle,'FontSize',12);
    hold on
    plotm(coast.lat,coast.long,'color','k')
    plotm([states.Lat],[states.Lon],'k')
    axesm('eqdcylin','MapLatLimit',[latmin latmax],'MapLonLimit',[lonmin lonmax],...
        'grid','on','PLineLocation',[0:5:90],'MLineLocation',[-180:5:180],...
        'MLineFill',50,'PlineFill',50,'glinewidth',0.3)
    if (cols(iplot) == 1) 
         axesm('eqdcylin','MapLatLimit',[latmin latmax],'MapLonLimit',[lonmin lonmax],...
        'grid','on','PLineLocation',[0:5:90],'MLineLocation',[-180:5:180],...
        'ParallelLabel','on','MeridianLabel','off','LabelFormat','none',...
        'PlabelMeridian',-140,'FontSize',8,'PLabelLocation',[20,30,40,50],...
        'MLineFill',50,'PlineFill',50,'glinewidth',0.3) 
    end
    if (rows(iplot) == num_rows)
         axesm('eqdcylin','MapLatLimit',[latmin latmax],'MapLonLimit',[lonmin lonmax],...
        'grid','on','PLineLocation',[0:5:90],'MLineLocation',[-180:5:180],...
        'ParallelLabel','off','MeridianLabel','on','LabelFormat','none',...
        'PlabelMeridian',-140,'FontSize',8,'MlabelParallel',20,...
        'MlabelLocation',[-135,-125,-115,-105],...
        'MLineFill',50,'PlineFill',50,'glinewidth',0.3)
    end    
    tightmap;

    % Add text box with frequency and mean pattern correlation
    Ptitle1 = ['F = ',   num2str(100*pat_freq(iplot),'%.2f'),'%'];
    %Ptitle2 = ['\rho = ',num2str(pat_corr(iplot)   ,'%.2f')    ];
    ax = gca;
    dx = ax.XLim(2)-ax.XLim(1);
    dy = ax.YLim(2)-ax.YLim(1);
    text(ax.XLim(1)+dx*0.05,ax.YLim(1)+dy*0.15,Ptitle1)
    %text(ax.XLim(1)+dx*0.05,ax.YLim(1)+dy*0.05,Ptitle2)
    
end
axes('Position', [0.20 0.00 0.65 0.50], 'Visible', 'off','CLim',[minP maxP]);
ct=colorbar('horizontal');
title(ct,cbarTitle,'FontSize',10);

end