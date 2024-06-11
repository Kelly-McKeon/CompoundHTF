clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data
load data.mat %observations
load bootstrap_data.mat %raw output from bootstrapping (1000 iterations)
load bstats_data.mat %summarized bootstrap data

%% Sample Event (Figure 1)

%get dates to plot
storm_dn = datenum(2021,10,24);
daterange = storm_dn-2 : storm_dn+20;
ind = find(datenum(data(1).dn) >= min(daterange) & datenum(data(1).dn) <= max(daterange));

%plot
sites = [3,6,8];
for i = 1:numel(sites);
    nn=sites(i);
    dn = data(nn).dn(ind);
    sl = data(nn).sl_daily(ind);
    td = data(nn).td_daily(ind);
    NTOOE = data(1).NTOOE(ind) * data(nn).a;
    q = data(1).q(ind) * data(nn).b;
    
    if i==1;
        figure
    end
    
    subplot(3,1,i);
    hold on
    yyaxis left
    hd1 = plot(dn, sl, 'LineWidth',2);
    hd2 = plot(dn, td, 'LineWidth',2);
    hd4 = plot(dn, NTOOE, 'LineWidth',2);
    yline(data(nn).thresh, 'LineWidth',2);
    ylabel('SL (m, MHHW)','FontSize',12)
    ylim([-0.2 1.5])
    xlim([min(daterange) max(daterange)]);
    datetick

    yyaxis right
    hd3=plot(dn, q,'LineWidth',2);
    ylim([-0.2 1.5])
    xlim([min(daterange) max(daterange)]);
    datetick
    ylabel('SL from Q (m, MHHW)','FontSize',12)
    hold off
    box on
    grid on
 

    legend([hd1,hd2,hd3,hd4],{'Total SL','Tide','Q','NTOOE'})
    title([convertCharsToStrings(data(nn).na)])
    
    sgtitle('October 21, 2021 Storm')
    
end

%% Map (Figure 2)

%load GEBCO Bathymetry
bathy = openNetCDF('gebco.nc');

%plot whole thing
lat=bathy.lat;
lon=bathy.lon;
h=double(bathy.elevation);
z1=h';
z1(z1 >= 0) = nan;

[lagrid, logrid] = meshgrid(lon,lat);

latlim = [data(5).la-1 data(5).la+1];
lonlim = [data(5).lo-0.5 data(5).lo+1];

%get latitudes and longitudes of points
longs = [];
lats = [];
for nn = 2:numel(data);
    longs = [longs; data(nn).lo];
    lats = [lats; data(nn).la];
end

b_longs = [-75.1193, -74.7780556];
b_lats = [38.7828, 40.22166667];

%plot the map
figure
m1 = worldmap([data(5).la-1 data(5).la+1], [data(5).lo-0.5 data(5).lo+1])
geoshow(logrid,lagrid,z1,'DisplayType','surface')
states = shaperead('usastatelo', 'UseGeoCoords', true);
geoshow('usastatelo.shp','facecolor',[0.9 0.9 0.9],'edgecolor',[.5 .5 .5])

hold on
setm(gca,'FontSize',14);
cb=colorbar
cmap=cbrewer2('seq',['Blues'],25);
colormap(flipud(cmap));
clim([-25 -1.5])
cb.FontSize = 16
cb.Label.String = "depth (m)"
northarrow("latitude", 40.25, "longitude", -74.8)
scaleruler('Units','km','RulerStyle','patches')

%plot gauges
s1=scatterm(lats,longs,100,'filled','k');
s2=scatterm(b_lats,b_longs, 100, 'diamond','filled');

legend([s1,s2],'NOAA Gauges','Boundary Gauges')

%inset map
m2 = axes('Position',[0.30 0.7 0.2 0.2],'Visible','off');
usamap({'VA','VT'})
plabel off
mlabel off
setm(m2,'FFaceColor','w')
geoshow(states,'FaceColor',[0.9 0.9 0.9],'Parent',m2)
plotm(latlim([1 2 2 1 1]),lonlim([2 2 1 1 2]), ...
    'Color','red','LineWidth',1)

f=gcf;
f.Renderer='painters'
saveas(f,'map.eps')
%% Regression Coefficients (Figure 3)

%create vectors by site to plot with

%observations
rsq = []; %multiple regression
a = [];
b = [];
sl_rsq = []; %NTOOE single regression rsq
q_rsq = []; %q single regression rsq
sum_r = []; % sum of both

%bootstrapped
bootr_sl = []; %NTOOE single regression
bootr_sl_std =[];

bootr_q = []; %Q single regression
bootr_q_std =[];

bootr_sum = []; % single regression sum
bootr_sum_std =[];

bootr = []; %multiple regression
bootr_std =[];

boota = [];
boota_std =[];

bootb = [];
bootb_std = [];

names = [];

for nn = 2:numel(data);
    %%%%% multiple regression
    %observations
    rsq = [rsq; data(nn).r]; %rsq from observations
    a = [a; data(nn).a]; %a from observations
    b = [b; data(nn).b]; %b from observations    

    %bootstrapped
    bootr = [bootr; bstats(nn).meanr]; %r from bootstrapping
    bootr_std = [bootr_std; bstats(nn).stdr]; %rsq std

    boota = [boota; bstats(nn).meana]; % a from bootstrapping
    boota_std = [boota_std; bstats(nn).stda]; % a std

    bootb = [bootb; bstats(nn).meanb]; % b from bootstrapping
    bootb_std = [bootb_std; bstats(nn).stdb]; %b std

    %%%%% single regressions
    %observations
    sl_rsq = [sl_rsq; data(nn).sl_r]; %NTOOE single regression rsq
    q_rsq = [q_rsq; data(nn).q_r]; %q single regression rsq
    sum_r = [sum_r; data(nn).r_sum]; % sum of both
    
    %bootstrapped
    bootr_sl = [bootr_sl; bstats(nn).meanr_sl]; %NTOOE single regression boostrapped
    bootr_sl_std =[bootr_sl_std; bstats(nn).stdr_sl];
    
    bootr_q = [bootr_q; bstats(nn).meanr_q]; %Q single regression bootstrapped
    bootr_q_std =[bootr_q_std; bstats(nn).stdr_q];

    bootr_sum = [bootr_sum; bstats(nn).meanr_sum]; %single regression sums bootstrapped
    bootr_sum_std =[bootr_sum_std; bstats(nn).stdr_sum];
    
    %just get the names in plotting condition
    names=[names;convertCharsToStrings(data(nn).na)];
end

%make names in categories to use later 
plot_names=categorical(names);
plot_names=reordercats(plot_names,names);

%%%%% now actually plot it
% Scatter plot of regression coefficients with errors
figure

%%%%% r squares
subplot(1,3,1);
hold on
%multiple r square
s1 = scatter(plot_names, rsq, 90,'filled','k')
e1 =errorbar(plot_names, rsq, bootr_std*-2, bootr_std*2,'LineWidth',1);
e1.Marker='o'
e1.MarkerFaceColor='k'
e1.MarkerEdgeColor='k'
e1.Color='k'
e1.LineStyle='none';

% sum r square
s2 = scatter(plot_names, sum_r, 90, '^', 'filled', 'MarkerFaceColor', '#7E2F8E')
e2 =errorbar(plot_names, bootr_sum, bootr_sum_std*-2, bootr_sum_std*2,'LineWidth',1);
e2.Marker='^'
e2.MarkerFaceColor='#7E2F8E'
e2.MarkerEdgeColor='#7E2F8E'
e2.Color='#7E2F8E'
e2.LineStyle='none';

%surge only r square 
s3 = scatter(plot_names, sl_rsq, 90, '^', 'filled', 'MarkerFaceColor', "#A2142F");
e3 =errorbar(plot_names, bootr_sl, bootr_sl_std*-2, bootr_sl_std*2,'LineWidth',1);
e3.Marker='^'
e3.MarkerFaceColor=	"#A2142F"
e3.MarkerEdgeColor=	"#A2142F"
e3.Color="#A2142F"
e3.LineStyle='none';

%q only r square
s4 = scatter(plot_names, q_rsq, 90,'^','filled','MarkerFaceColor',"#0072BD")
e4 =errorbar(plot_names, bootr_q, bootr_q_std*-2, bootr_q_std*2,'LineWidth',1);
e4.Marker='^'
e4.MarkerFaceColor=	"#0072BD"
e4.MarkerEdgeColor=	"#0072BD"
e4.Color= "#0072BD"
e4.LineStyle='none';

ax=gca;
ax.FontSize=16;
grid on
box on
ylabel('{\it R^2}','FontSize',16,'FontWeight','bold')
legend([s1,s2,s3,s4],'{\it R^2} multiple', '{\it R^2} sum', '{\it R^2 NTOOE}', '{\it R^2 Q}')
% title('r^2','FontSize',14)
hold off


%%% a coefficient
subplot(1,3,2);
scatter(plot_names, a, 90,'filled','MarkerFaceColor', "#A2142F")
hold on
e=errorbar(plot_names, boota, boota_std*-2, boota_std*2,'LineWidth',1);
hold on
e.Marker='o'
e.MarkerFaceColor="#A2142F"
e.MarkerEdgeColor="#A2142F"
e.Color="#A2142F"
e.LineStyle='none';
ylabel('{\it a}','FontSize',14,'FontWeight','bold')
ax=gca;
ax.FontSize=16;
grid on
box on
% title('Surge Regression Coefficient','FontSize',14)
hold off

%%%% b coefficient
subplot(1,3,3)
scatter(plot_names, b, 90, 'filled','MarkerFaceColor',"#0072BD")
hold on
e=errorbar(plot_names, bootb, bootb_std*-2, bootb_std*2, 'LineWidth',1);
e.Marker='o'
e.MarkerFaceColor="#0072BD"
e.MarkerEdgeColor="#0072BD"
e.Color="#0072BD"
e.LineStyle='none';
ax=gca;
ax.FontSize=16;
grid on
box on
ylabel('{\it b}','FontSize',16,'FontWeight','bold'); 
% title('Discharge Regression Coefficient','FontSize',14)


%% Residual Matrix (Figure S1)

%make matrix of rsq residuals
residual_matrix=zeros(8,8); %empty matrix for the rsq of each site pair

for i=2:numel(data);
    for ii=2:numel(data);      
        x=data(i).residuals - nanmean(data(i).residuals); % x is site 1
        y=data(ii).residuals - nanmean(data(ii).residuals); % y is site 2

        mdl=fitlm(x,y,'Linear');
        residual_matrix(i-1,ii-1) = mdl.Rsquared.Ordinary; %fill in matrix according to the site pair
    end
end

%plot the matrix
figure
imagesc(residual_matrix);
ax=gca;
ax.XTick=[1:8];
ax.YTick=[1:8];
set(ax,'XTickLabel',flip(plot_names));
set(ax,'YTickLabel',flip(plot_names));
ax.YDir='reverse';
ax.XDir='reverse';
ax.XTickLabelRotation=45;
ax.YTickLabelRotation=45;
ax.FontSize=14;
cmap=cbrewer2('seq',['Blues'],11);
colormap(cmap);
c=colorbar;
clim([0 1])
c.Label.String='\it R^2';
c.Label.FontSize=20;
c.Label.FontWeight='bold'
c.FontSize=20;
title('mean{\it R^2} of residuals')


%% Bar Graph of Flood Types (Figure 4c)

%make vectors for plotting
singles=[];
doubles=[];
triples=[];
four=[];

singles_std = [];
doubles_std = [];
triples_std = [];
four_std = [];

% use bootstrapped means of each flood type rather than observations
for nn=2:numel(bstats);
    %concatenate means
    singles = [singles; bstats(nn).singles_count_mean];
    doubles = [doubles; bstats(nn).doubles_count_mean];
    triples = [triples; bstats(nn).triples_count_mean];
    four = [four; bstats(nn).four_count_mean];

    %concatenate standard deviations
    singles_std = [singles_std; bstats(nn).singles_count_std];
    doubles_std = [doubles_std; bstats(nn).doubles_count_std];
    triples_std = [triples_std; bstats(nn).triples_count_std];
    four_std = [four_std; bstats(nn).four_count_std];

    %make matrixes for doing a stacked bar graph
    types(:,nn-1) = [bstats(nn).singles_count_mean; bstats(nn).doubles_count_mean; bstats(nn).triples_count_mean; bstats(nn).four_count_mean];
    err(:,nn-1) = [bstats(nn).singles_count_std*2; bstats(nn).doubles_count_std*2; bstats(nn).triples_count_std*2; bstats(nn).four_count_std*2];

end


%bar graph along river
figure
b=bar(plot_names,types','stacked');
hold on
e=errorbar(cumsum(types)',err','.k','LineWidth',1);


patchHand = findobj(b, 'Type', 'Patch'); 
b(1).FaceColor = "#FF0000";
b(2).FaceColor = "#0000FF";
b(3).FaceColor = "#0000CD";
b(4).FaceColor = "#00008B"; 
legend('extremes','2 components','3 components','4 components','FontSize',12);
ax=gca;
grid on
box on
ax.FontSize=18;
ylabel('Number of Floods','FontSize',18,'FontWeight','bold')

%% Bar Graph of Seasons (Figure S2)

%save seasonal data in vectors for plotting
win = [];
spr = [];
summ = [];
fall = [];
total=[];

for nn = 2:numel(data);
    win = [win, floods(nn).win];
    spr = [spr, floods(nn).spr];
    summ = [summ, floods(nn).summ];
    fall = [fall, floods(nn).fall];

    total = [total, numel(floods(nn).wl)];
end

%make matrix for stacked bar graph
seasons = [win;spr;summ;fall];
win_perc = win/total*100;
spr_perc = spr./total*100;
summ_perc = summ./total*100;
fall_perc = fall./total*100;

figure
b=bar(plot_names, seasons', 'stacked');
legend('winter','spring','summer','fall');

%% Number and Magnitude of Flood Types (Figure S3)

%concatenate water levels into one vector
wl1 = []; %single component floods
wl1_std =[]; %single component standard deviation 

wl2 = []; %2 components
wl2_std = [];

wl3 = []; %3 components
wl3_std = [];

wl4 = []; %4 components
wl4_std = [];

for nn = 2:numel(bstats)
    %concatenate means
    wl1 = [wl1; bstats(nn).singles_wl_mean];
    wl2 = [wl2; bstats(nn).doubles_wl_mean];
    wl3 = [wl3; bstats(nn).triples_wl_mean];
    wl4 = [wl4; bstats(nn).four_wl_mean];

    %concatenate standard deviations
    wl1_std = [wl1_std; bstats(nn).singles_wl_std];
    wl2_std = [wl2_std; bstats(nn).doubles_wl_std];
    wl3_std = [wl3_std; bstats(nn).triples_wl_std];
    wl4_std = [wl4_std; bstats(nn).four_wl_std];
end

%plot
figure
subplot(1,2,1) %first subplot is flood counts
hold on
s1=scatter(plot_names, singles, 150, 'o', 'filled');
s1.MarkerFaceColor="#D95319";
e=errorbar(plot_names, singles, singles_std*-2, singles_std*2,'LineWidth',1);
e.Color="#D95319";

s2=scatter(plot_names, doubles, 150, '^', 'filled');
s2.MarkerFaceColor="#0072BD";
e=errorbar(plot_names, doubles, doubles_std*-2, doubles_std*2,'LineWidth',1);
e.Color="#0072BD";

s3=scatter(plot_names, triples, 150,'square','filled');
s3.MarkerFaceColor="#77AC30";
e=errorbar(plot_names, triples, triples_std*-2, triples_std*2,'LineWidth',1);
e.Color="#77AC30";

s4=scatter(plot_names, four, 150, 'diamond', 'filled');
s4.MarkerFaceColor="#7E2F8E";
e=errorbar(plot_names,four, four_std*-2, four_std*2,'LineWidth',1);
e.Color="#7E2F8E";

legend([s1,s2,s3,s4],'Single','Double','Triple','Four')
% title('Flood Days by Flood Type')
ylabel('Flood Counts')
ax=gca
ax.FontSize=14;
grid on
box on

%second subplot is magnitudes
subplot(1,2,2)
hold on
s1=scatter(plot_names, wl1, 150, 'o', 'filled');
s1.MarkerFaceColor="#D95319";
e=errorbar(plot_names, wl1 , wl1_std*-2, wl1_std*2,'LineWidth',1);
e.Color="#D95319";

s2=scatter(plot_names, wl2, 150, '^', 'filled');
s2.MarkerFaceColor="#0072BD";
e=errorbar(plot_names, wl2 , wl2_std*-2, wl2_std*2,'LineWidth',1);
e.Color="#0072BD";

s3=scatter(plot_names, wl3, 150, 'square','filled');
s3.MarkerFaceColor="#77AC30";
e=errorbar(plot_names, wl3 , wl3_std*-2, wl3_std*2,'LineWidth',1);
e.Color="#77AC30";

s4=scatter(plot_names, wl4, 150, 'diamond','filled');
s4.MarkerFaceColor="#7E2F8E";
e=errorbar(plot_names, wl4 , wl4_std*-2, wl4_std*2,'LineWidth',1);
e.Color="#7E2F8E";

legend([s1,s2,s3,s4],'Single','Double','Triple','Four')

% title('Flood Magnitude by Flood Type')
ylabel('Sea Level (m, MHHW)')
ax=gca
ax.FontSize=14;
grid on
box on



%% Plot matrix of percents (Figure 4b)
%use bootstrapped means, not observations

%concatenate bootstrapped percent means into one matrix
percents = zeros(8,width(bstats(nn).percents_mean));
for nn=2:numel(bstats);
    percents(nn,:) = bstats(nn).percents_mean;
end

labels=categorical({'\it T',' \it Q','\zeta_0','\epsilon',' \it Q + T',...
    '\zeta_0 + \it T','\epsilon + \it T','\it Q + \epsilon','\zeta_0 + \epsilon',...
    '\it Q + \zeta_0','\it Q + T + \epsilon','\zeta_0 + \it T + \epsilon',...
    '\it Q + \zeta_0 + \it T','\it Q + \zeta_0 + \epsilon','\it Q + \zeta_0 + \it T + \epsilon'});

figure
pcolor(flipud(percents(1:9,1:16))); %only plot the first 14 columns because the 15th is the sum
ax=gca;
ax.YDir='reverse';
ax.XTickLabelRotation=45;
ax.YTickLabelRotation=45;
ax.FontSize=18;
ax.FontName = "Georgia"
ax.XTick=[1:18];
set(ax,'XTickLabel',(labels(1:15)));
set(ax,'YTickLabel',(flip(plot_names)));
% 
cmap=cbrewer2('seq',['Blues'],15);
colormap(cmap)
clim([0 25])
c=colorbar;
c.Location='southoutside';
c.Label.FontSize=16;
c.Label.FontWeight='bold';
c.FontSize=14;
c.Label.String='Percentage of Total Floods';
% clear ax

%% Percents at Philadelphia (Figure 4a)
%reorder labels because scatter wants to plot them categorical
labels = reordercats(labels,{'\it T',' \it Q','\zeta_0','\epsilon',' \it Q + T',...
    '\zeta_0 + \it T','\epsilon + \it T','\it Q + \epsilon','\zeta_0 + \epsilon',...
    '\it Q + \zeta_0','\it Q + T + \epsilon','\zeta_0 + \it T + \epsilon',...
    '\it Q + \zeta_0 + \it T','\it Q + \zeta_0 + \epsilon','\it Q + \zeta_0 + \it T + \epsilon'});

figure
e=errorbar(labels, bstats(7).percents_mean(1:15), (bstats(7).percents_lowCI(1:15)-bstats(7).percents_mean(1:15)),...
    (bstats(7).percents_upCI(1:15)-bstats(7).percents_mean(1:15)), 'LineStyle', 'none')
e.Marker = 'square'
e.Color = "k"
e.MarkerSize = 20;
e.MarkerFaceColor = "#0072BD"
e.MarkerEdgeColor = 'k'
e.LineWidth = 1.5;
ax = gca
ax.FontSize = 16
grid on
ylabel('Percent of Floods')
title('Philadelphia')



%% Q + NTOOE  at all sites (Figure 4d)
%this is variable column 10 in the percents matrix

QN_mean = [];
QN_lowCI = [];
QN_highCI = [];
for i = 2:numel(bstats);
    QN_mean = [QN_mean; bstats(i).percents_mean(10)];
    %error bar is measured in distance from the mean, so remove the mean
    %from the confidence intervals to get it to plot right

    QN_lowCI = [QN_lowCI; bstats(i).percents_lowCI(10) - bstats(i).percents_mean(10)];
    QN_highCI =[QN_highCI; bstats(i).percents_upCI(10) - bstats(i).percents_mean(10)];
end

figure
e=errorbar(plot_names, QN_mean, QN_lowCI, QN_highCI, 'LineStyle', 'none')
e.Marker = 'square'
e.Color = "k"
e.MarkerSize = 20;
e.MarkerFaceColor = "#0072BD"
e.MarkerEdgeColor = 'k'
e.LineWidth = 1.5;
ax = gca
ax.FontSize = 16
ax.YAxisLocation = "right"
xtickangle(-45)
grid on
box on
ylim([0 30])
ylabel('Percent of Floods')
title('Discharge and NTOOE Floods')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Repeat for Number of Floods Rather than Percents (Figure S4b)
%use bootstrapped means, not observations

%concatenate bootstrapped percent means into one matrix
numfl = zeros(8,width(bstats(nn).numfl_mean));
for nn=2:numel(bstats);
    numfl(nn,:) = bstats(nn).numfl_mean;
end

labels=categorical({'T','Q','NTOOE','R','Q + T','NTOOE + T','R + T','Q + R',...
    'NTOOE + R','Q + NTOOE','Q + T + R','NTOOE + T+ R','Q + NTOOE + T',...
    'Q + NTOOE + R','Q + NTOOE + T+ R'});

figure
pcolor(flipud(numfl(1:9,1:16))); %only plot the first 14 columns because the 15th is the sum
ax=gca;
ax.YDir='reverse';
ax.XTickLabelRotation=45;
ax.YTickLabelRotation=45;
ax.FontSize=18;
ax.FontName = "Georgia"
ax.XTick=[1:18];
set(ax,'XTickLabel',(labels(1:15)));
set(ax,'YTickLabel',(flip(plot_names)));

cmap=cbrewer2('seq',['Blues'],16);
colormap(cmap)
clim([0 80])
c=colorbar;
c.Location='southoutside';
c.Label.FontSize=16;
c.Label.FontWeight='bold';
c.FontSize=14;
c.Label.String='Number of Total Floods';

%% Numbers at Philadelphia (Figure S4a)
%reorder labels because scatter wants to plot them categorical
labels = reordercats(labels,{'T','Q','NTOOE','R','Q + T','NTOOE + T','R + T',...
    'Q + R','NTOOE + R','Q + NTOOE','Q + T + R','NTOOE + T+ R','Q + NTOOE + T',...
    'Q + NTOOE + R','Q + NTOOE + T+ R'});

figure
e=errorbar(labels, bstats(7).numfl_mean(1:15), (bstats(7).numfl_lowCI(1:15) - bstats(7).numfl_mean(1:15)),...
    (bstats(7).numfl_upCI(1:15)-bstats(7).numfl_mean(1:15)), 'LineStyle', 'none')
e.Marker = 'square'
e.Color = "k"
e.MarkerSize = 20;
e.MarkerFaceColor = "#0072BD"
e.MarkerEdgeColor = 'k'
e.LineWidth = 1.5;
ax = gca
ax.FontSize = 16
ax.YAxisLocation = "right"
xtickangle(-45)
grid on
box on
ylabel('Number of Floods')
title('Philadelphia')

%% Q + NTOOE Numbers at all sites (Figure S4d)
%this is variable column 10 in the  matrix


%do this for just Q + NTOOE which is variable 10
QN_mean = [];
QN_std = [];
QN_lowCI = [];
QN_upCI = [];
for i = 2:numel(bstats);
    QN_mean = [QN_mean; bstats(i).numfl_mean(10)];
    QN_std = [QN_std; bstats(i).numfl_std(10)];
    QN_lowCI = [QN_lowCI; bstats(i).numfl_lowCI(10) - bstats(i).numfl_mean(10)];
    QN_upCI = [QN_upCI; bstats(i).numfl_upCI(10) - bstats(i).numfl_mean(10)];
end

figure
e=errorbar(plot_names, QN_mean, QN_lowCI, QN_upCI, 'LineStyle', 'none')
e.Marker = 'square'
e.Color = "k"
e.MarkerSize = 20;
e.MarkerFaceColor = "#0072BD"
e.MarkerEdgeColor = 'k'
e.LineWidth = 1.5;
ax = gca
ax.FontSize = 16
ax.YAxisLocation = "right"
xtickangle(-45)
grid on
box on
ylabel('Number of Floods')
title('Discharge and NTOOE Floods')



