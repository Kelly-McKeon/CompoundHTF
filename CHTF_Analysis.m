clear all;
close all;

%timer
start = tic;
%% Set up Data

%NOAA stations
%first one is Lewes
noaaID = [8557380 8555889 8537121 8551910 8551762 8540433 8545240 8539094 8548989];
save(['NOAAfileID.mat'],'noaaID')
t0=2002;%start year 
tf=2023;%end year 

%USGS stations
usgsID = string("01463500");
gage = ["Trenton, NJ"];
save(['UGSSfileID.mat'],'usgsID')
cnt=1;
parameter{cnt}='00060';  %discharge parameters is 00060
varname{cnt}='Q';
vstr{cnt}='Discharge(cfs)';
dstart='2002-01-01';
dend='2023-12-31';

%% Download NOAA Data and run UTide- Comment out once it's saved
% All sl/td data is in GMT time, metric units, and MHHW datum
% UTide must be in GMT to work 
% flood thresholdsin STND datum


tstart = tic;
for id=1:numel(noaaID),disp(['Now on Site ', num2str(id)]);
    tic

        [datum]=noaaDatums(noaaID(id));
        [sl,dn]=noaaSealevel_GMT(noaaID(id),t0-1,tf+1); %get one year before and one year after for UTide
        [flood]=noaaFlood(noaaID(id));
        [noaa_td, dn_local] = noaaTide(noaaID(id),t0,tf); %NOAA tide to compare to UTide

        % detrend sl data first
        ind = ~isnan(sl);
        P = polyfit(dn(ind), sl(ind), 1); %fit linear trend
        m = P(1); %slope
        b = P(2); %intercept
        trend = m*dn + b; %make trend line
        sldt = sl-trend; %detrended
        clearvars ind


        % fit year by year analysis
        tdyr=nan*sl; %empty array for tide data to go in
        for i = t0:tf, disp(['Year ', num2str(i)]);
            [y,~,~] = datevec(dn);
            ind = find( (y >= i-1) & (y <= i+1)); %index of 3 year period to estimate tides

            if id ==2 %use entire timeseries at brandywine because has so many data gaps
                ind = find( (y >= t0) & (y <= tf));
            end

            if id == 6 && (i == 2015) %Marcus Hook missing year 2015
                ind = find((y >= i-3) & (y <= i+3));
            end

            if mean(~isnan(sldt(ind))) > 0.75 %only do full year analysis if >75% hourly data are available
                coef = ut_solv(dn(ind), sldt(ind), [], datum.lat, 'auto');
                indy = find(y==i); %only grab the middle year 
                [td_temp, ~] = ut_reconstr(dn(indy), coef);
                tdyr(indy) = td_temp;
            else
                disp(['not enough data', num2str(i)])
                ind = find( (y >= i-2) & (y <= i+2)); %if not enough data points expand window to 5 years
                coef = ut_solv(dn(ind), sldt(ind), [], datum.lat, 'auto');
                indy = find(y==i); %only grab the middle year 
                [td_temp, ~] = ut_reconstr(dn(indy), coef);
                tdyr(indy) = td_temp;
            end

        end

        ind2 = find( (y >= t0) & (y <= tf) ); %get only the study years
        % get seasonal component
        coef = ut_solv(dn(ind2), sldt(ind2)-tdyr(ind2), [], datum.lat, [{'SA'};{'SSA'}]);
        [td_seas, ~] = ut_reconstr(dn(ind2), coef);

        % add together
        td = tdyr(ind2) + td_seas;
        dn = dn(ind2);
        sl = sl(ind2);
        sldt = sldt(ind2);
        trend = trend(ind2);

        save(['noaa_tidegauge_yearly_t',num2str(datum.id),'.mat'],'datum','sl','sldt','trend','m','b','dn','td', 'flood', 'noaa_td', 'dn_local')
    toc

end 
tstop = toc(tstart);

disp(['time to run tidal analysis is ', num2str(tstop/60), ' minutes']) 


%% Download USGS data 

for id=1:numel(usgsID);
    tic
    datatable = getusgsdata((usgsID(id)), parameter, dstart, dend, varname);
    data_array = table2array(datatable(:,4));
    date_array = table2array(datatable(:,3));
    usgsdata(id).gage = gage(id);
    usgsdata(id).site_no = num2str(usgsID(id));
    usgsdata(id).date = date_array; 
    usgsdata(id).qaqc = datatable(1:end,5);

    %usgs data
    usgsdata(id).Q = data_array;
    usgsdata(id).Q = usgsdata(id).Q/35.3147; %convert from cfs to cms
     
    save(['usgs_gage_',num2str(usgsID(id)),'.mat'],'usgsdata')
    toc
end


%% Load in Downloaded Data and Create Data Structure

load NOAAfileID.mat 

for nn=1:numel(noaaID), disp(['Now on file number ', num2str(nn)])
    load(['noaa_tidegauge_yearly_t',num2str(noaaID(nn)),'.mat']) %load data from saved mat files


    %convert to daily 
    td = td+trend; %put the sea level rise in the tidal component
    sl_daily = max(reshape(sl,24,numel(sl)/24,[],1)); %reshape sea level to be daily max
    td_daily = max(reshape(td,24,numel(td)/24,[],1)); %reshape tide to be daily max 
    NTOOE = sl_daily-td_daily; %NTOOE are maximum daily sea level - maximum daily tidal prediction (remember SLR is in the tide here)


    %save site names in structure
    data(nn).na = datum.name'; 
    if nn==8;
        data(nn).na = 'Burlington';
    end

    %save data in structure
    data(nn).lo = datum.lon;
    data(nn).la = datum.lat;
    data(nn).date = dn_local; %hourly date nums in local timezone
    data(nn).dn = floor(mean(reshape(dn_local,24,numel(dn_local)/24,1))); %get daily datenums in same timezone as discharge
    data(nn).sl = sl; %hourly sea level (MHHW)
    data(nn).td = td; %hourly tide (MHHW) 
    data(nn).sl_daily = sl_daily'; %daily maximum sea level (MHHW)
    data(nn).td_daily = td_daily'; %daily maximum tide (MHHW)
    data(nn).NTOOE = NTOOE'; %daily NTOOE (MHHW)
    data(1).q = usgsdata(1).Q; %q

    % NWS flood thresholds
    data(nn).minor_thresh = flood.nws_minor - datum.MHHW ; %convert STND datum to MHHW 
    data(nn).moderate_thresh = flood.nws_moderate - datum.MHHW ;
    data(nn).major_thresh = flood.nws_major - datum.MHHW ;

    % Get WL percentiles
    data(nn).prct95 = prctile(data(nn).sl_daily,95); 
    data(nn).prct99 = prctile(data(nn).sl_daily,99);

    %MHHW above MSL
    data(nn).MHHW = datum.MHHW - datum.MSL;
    
    %set Brandywine and Ship John Shoal thresholds as Mahmoudi et al. 2024 values
    data(2).minor_thresh = 0.5166;
    data(3).minor_thresh = 0.4650;
end

%% Kendall's Tau

%compute for observations

for nn=2:numel(data);
    %remove means before computing
    q = data(1).q - nanmean(data(1).q);
    NTOOE = data(1).NTOOE - nanmean(data(1).NTOOE);
    NTOOEx = data(nn).NTOOE - nanmean(data(nn).NTOOE);
    
    % Q
    [r1,p1] = corr(NTOOEx, q,'type','Kendall','rows','complete');
    data(nn).kt_Q = r1; %tau statistic
    data(nn).kt_Qp = p1; %p value
    % NTOOE
    [r2,p2] = corr(NTOOEx, NTOOE, 'type', 'Kendall', 'rows', 'complete');
    data(nn).kt_NTOOE = r2; %tau statistic
    data(nn).kt_NTOOEp = p2; %p value

end


%% Multiple Regression 
for nn=2:numel(data);
    %remove nans and then means before regression
    y = data(nn).NTOOE;
    x1 = data(1).NTOOE;
    x2 = data(1).q;

    %only regress where there are data points in all three timeseries
    iii = find(~isnan(y) & ~isnan(x1) & ~isnan(x2));
    y = y(iii); x1=x1(iii); x2=x2(iii);
    y = y - mean(y); %remove means before regressing
    x1 = x1 - mean(x1);
    x2 = x2 - mean(x2);
    X = [x1,x2];

    %solve
    mdl = fitlm(X,y,'linear');
    a = mdl.Coefficients.Estimate(2);
    b = mdl.Coefficients.Estimate(3);
    r = mdl.Rsquared.Ordinary;
    residuals = y-a*x1-b*x2;


    % make new residual variable that's the same size as the timeseries
    data(nn).residuals = nan(size(data(nn).NTOOE));
    data(nn).residuals(iii) = residuals;
    data(nn).a = a;
    data(nn).b = b;
    data(nn).r = r;

    %just use this line to make sure the mean of everything is zero 
    data(nn).res_mean = mean(residuals);
    data(nn).NTOOE_mean = mean(a*x1);
    data(nn).Q_mean = mean(b*x2);

    %% Nash-Sutcliffe Equation
    Qo = y; %observed
    Qm = a*x1 + b*x2; %simulated

    n=[];
    d=[];
    for i = 1:numel(y)
        n_temp = (Qo(i) - Qm(i))^2 ;
        d_temp = (Qo(i) - mean(Qo))^2;

        n = [n, n_temp];
        d = [d, d_temp];

    end
    data(nn).NSE = 1 - (sum(n)/sum(d));

    %% Kling-Gupta Equation

    us = a*x1 + b*x2; %simulated
    uo = y; %observed

    pr = corrcoef(us,uo); % Pearson Correlation Coefficient
    pr = pr(1,2);
    beta = mean(us,'omitnan') / mean(uo,'omitnan'); %divide by zero so ignore this term
    alpha = var(us,'omitnan') / var(uo,'omitnan');

    data(nn).KGE = 1 - sqrt( ((pr-1)^2) + ((alpha-1)^2) );
    
%% Single Regressions
    %NTOOE
    sl_mdl = fitlm(x1,y,'Linear');
    data(nn).sl_r = sl_mdl.Rsquared.Ordinary;
        
    %discharge
    q_mdl = fitlm(x2,y,'Linear');
    data(nn).q_r = q_mdl.Rsquared.Ordinary;

    %sum 
    r_sum = sl_mdl.Rsquared.Ordinary + q_mdl.Rsquared.Ordinary;
    data(nn).r_sum = r_sum;

%% Flood Days

    %Find where daily maximum sea level over threshold
    ind = find(data(nn).sl_daily >= data(nn).minor_thresh);

    % Create Data Structure For Just Flood Data
    floods(nn).na = data(nn).na; %site names
    floods(nn).wl = data(nn).sl_daily(ind); %water level
    floods(nn).dn = data(nn).dn(ind); %date
    floods(nn).NTOOE = data(nn).NTOOE(ind); %NTOOE
    floods(nn).sl_Lewes = data(1).sl_daily(ind); %water level at Lewes during flood
    floods(nn).NTOOE_Lewes = data(1).NTOOE(ind); %NTOOE at Lewes during flood
    floods(nn).td = data(nn).td_daily(ind); %tide during flood
    floods(nn).q = data(1).q(ind); %discharge during flood
    floods(nn).mean_wl = mean(data(nn).sl_daily(ind),'omitnan'); %water level during flood
    floods(nn).residuals = data(nn).residuals(ind); %residuals on flood days

    %get how many floods in each season
    DJF=[12,1,2];
    MAM=[3,4,5];
    JJA=[6,7,8];
    SON=[9,10,11];

    [y,m]=datevec(floods(nn).dn);
    ind1=find(ismember(m,DJF));
    ind2=find(ismember(m,MAM));
    ind3=find(ismember(m,JJA));
    ind4=find(ismember(m,SON));

    floods(nn).win=numel(ind1);
    floods(nn).spr=numel(ind2);
    floods(nn).summ=numel(ind3);
    floods(nn).fall=numel(ind4);


end

%save data and clear variables to clean up workspace
save data.mat data floods
clearvars -except data floods start

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin Bootstrapping

%% bootstrapping loop (comment out once iterations are compiled)
load data.mat

y1 = data(1).NTOOE;
y2 = data(1).q;
i = 1000;

bstart = tic;
for nn=2:numel(data);
    tic
    x = data(nn).NTOOE;
    x2 = data(nn).sl_daily;
    x3 = data(nn).td_daily;
    thresh = data(nn).minor_thresh;

    [boot(nn).p95, boot(nn).a, boot(nn).b, boot(nn).r, boot(nn).nse, boot(nn).kge,...
        boot(nn).r_sl, boot(nn).r_q, boot(nn).r_sum, ...
        boot(nn).residuals, boot(nn).fd, boot(nn).numfl,...
        boot(nn).frac] = CHTF_bootstrapping(x,y1,y2,x2,x3,thresh,i);

    display(['now on boot ', num2str(nn)])
    toc
end
bstop = toc(bstart);
disp(['time to run entire bootstrap is ', num2str(bstop/60), ' minutes'])

%% save data and clear variables to clean up workspace
save bootstrap_data.mat boot 
clearvars -except start

%% load saved bootstrap data
load data.mat %this dataset has only the observations for (variables "floods" and "data")

%this dataset has the 1000 bootstrap iterations of the data
load bootstrap_data.mat 

%% Get summary statistics from the bootstrap data
% save these in bstats structure

for nn = 2:numel(boot);

    %%%%%%%%%%%%% coefficients
    
    %standard deviations (all)
    %multiple regression
    bstats(nn).stda = nanstd(boot(nn).a); % a coeff std
    bstats(nn).stdb = nanstd(boot(nn).b); % b coeff std
    bstats(nn).stdr = nanstd(boot(nn).r); % rsq std
    %single regressions
    bstats(nn).stdr_sl = nanstd(boot(nn).r_sl); % sl single regression
    bstats(nn).stdr_q = nanstd(boot(nn).r_q); % q single regression
    bstats(nn).stdr_sum = nanstd(boot(nn).r_sum); % single regression sum
    %kge and nse
    bstats(nn).std_kge = nanstd(boot(nn).kge); %kge std
    bstats(nn).std_nse = nanstd(boot(nn).nse); %nse std
   
    %means (all)
    %multiple regression
    bstats(nn).meana = nanmean(boot(nn).a); % a coeff mean
    bstats(nn).meanb = nanmean(boot(nn).b); % b coeff mean
    bstats(nn).meanr = nanmean(boot(nn).r); % rsq mean
    %single regressions
    bstats(nn).meanr_sl = nanmean(boot(nn).r_sl); % sl single regression
    bstats(nn).meanr_q = nanmean(boot(nn).r_q); % q single regression
    bstats(nn).meanr_sum = nanmean(boot(nn).r_sum); % single regression sum
    %kge and nse
    bstats(nn).mean_kge = nanmean(boot(nn).kge); %kge mean
    bstats(nn).mean_nse = nanmean(boot(nn).nse); %nse mean


    
    %%%%%%%%%%%%%%% residuals

    res_mean=[];
    res_std=[];
    res_acorr=[];
    
    fd_res_mean=[];
    fd_res_std=[];
    for ii = 1:numel(boot(nn).fd);  %get the 1000 iterations into one vector      
        
        %all (the mean of these should be zero, this is a sanity check)
        res_mean = [res_mean; boot(nn).residuals(ii).mean]; %compile means of the 1000 residuals timeseries
        res_std = [res_std; boot(nn).residuals(ii).std]; %compile standard deviations of the 1000 residuals timeseries
        %autocorrelation
        res_acorr = [res_acorr; boot(nn).residuals(ii).acorr]; %compile autocorrelations of the 1000 residual timeseries
        
        %flood days 
        fd_res_mean = [fd_res_mean; nanmean(boot(nn).fd(ii).res)];
        fd_res_std = [fd_res_std; nanstd(boot(nn).fd(ii).res)];
         
    end
    
    %compile mean/std of residuals
    %all
    bstats(nn).res_mean = mean(res_mean); %this is the mean of the means of the 1000 residual timeseries
    bstats(nn).res_std = mean(res_std); %this is the mean of the standard deviations of the 1000 residual timeseries
    bstats(nn).res_acorr = mean(res_acorr); %this is the mean of the 1000 autocorrelations from the 1000 residual timeseries

    %flood days
    bstats(nn).fd_res_mean = mean(fd_res_mean); %this is the mean of the mean of the residuals on flood days from the 1000 ts
    bstats(nn).fd_res_std = mean(fd_res_std); %this is the mean of the std of the residuals on flood days from 1000 ts

   
    
    
    
    %%%%%%%%%%%%%%%%%% water level contributions on flood days 

    %create empty vectors to populate in the loop
    fd_td_mean=[];
    fd_td_std=[];
    
    fd_NTOOE_mean=[];
    fd_NTOOE_std=[];

    fd_Q_mean=[];
    fd_Q_std=[];
    
    %get 1000 timeseries means and stds
    for ii = 1:numel(boot(nn).fd);                
        %tide
        fd_td_mean = [fd_td_mean; nanmean(boot(nn).fd(ii).td)]; %mean of flood day tides
        fd_td_std = [fd_td_std; nanstd(boot(nn).fd(ii).td)]; %std of the flood day tides

        %NTOOE
        fd_NTOOE_mean = [fd_NTOOE_mean; nanmean(boot(nn).fd(ii).cont_NTOOE)]; %mean of the flood day NTOEE contribution
        fd_NTOOE_std = [fd_NTOOE_std; nanstd(boot(nn).fd(ii).cont_NTOOE)]; %std of the flood day NTOOE contribution

        %Q
        fd_Q_mean = [fd_Q_mean; nanmean(boot(nn).fd(ii).cont_Q)]; %mean of the flood day Q contribution
        fd_Q_std = [fd_Q_std; nanstd(boot(nn).fd(ii).cont_Q)]; %std of the flood day Q contribution

    end    
    
    %get means and stds of the 1000 means
    bstats(nn).fd_mean_td = mean(fd_td_mean); %mean of the 1000 tide means
    bstats(nn).fd_std_td = std(fd_td_mean); %std of the means

    bstats(nn).fd_mean_NTOOE_cont = mean(fd_NTOOE_mean); %mean of the 1000 NTOOE mean contributions
    bstats(nn).fd_std_NTOOE_cont = std(fd_NTOOE_mean); %std of the 1000 NTOOE means

    bstats(nn).fd_mean_q_cont = mean(fd_Q_mean); %mean of the 1000 Q means
    bstats(nn).fd_std_q_cont = std(fd_Q_mean); %std of the 1000 Q means


    %%%%%%%%%%%%%%%% flood counts and magnitudes    
    fd_count = [];
    fd_mag = [];

    for ii=1:numel(boot(nn).fd);
        fd_count = [fd_count; boot(nn).fd(ii).count]; %number of floods in each of the 1000 timeseries
        fd_mag = [fd_mag; mean(boot(nn).fd(ii).sl)]; %mean magnitude of the floods in each of the 1000 timeseries
    end

    bstats(nn).fd_count_mean = mean(fd_count); %mean flood counts
    bstats(nn).fd_count_std = std(fd_count); %standard deviation of the flood counts

    bstats(nn).fd_mag_mean = mean(fd_mag); %means of flood day water levels
    bstats(nn).fd_mag_std = std(fd_mag); %standard deviation of means of flood day water levels 




    %%%%%%%%%%%%%%%%%% flood type prevalence and magnitude
    count1 = [];
    count2 = [];
    count3 = [];
    count4 = [];

    mag1 = [];
    mag2 = [];
    mag3 = [];
    mag4 = [];
    
    allcount = [];
    ccount = [];

    for ii = 1:numel(boot(nn).fd);
        count1 = [count1; boot(nn).fd(ii).singles]; %number of single component floods in each of the iterations
        mag1 = [mag1; boot(nn).fd(ii).singles_wl]; %mean magnitude of the single component floods in that iteration

        count2 = [count2; boot(nn).fd(ii).doubles]; %number of 2 component floods in each of the iterations
        mag2 = [mag2; boot(nn).fd(ii).doubles_wl]; %mean magnitude of the 2 component floods in that iteration

        count3 = [count3; boot(nn).fd(ii).triples]; %number of 3 component floods in each of the iterations
        mag3 = [mag3; boot(nn).fd(ii).triples_wl]; %mean magnitude of the 3 component floods in that iteration

        count4 = [count4; boot(nn).fd(ii).four]; %number of 4 component floods in each of the iterations
        mag4 = [mag4; boot(nn).fd(ii).four_wl]; %mean magnitude of the 4 component floods in that iteration

        allcount = [allcount; boot(nn).fd(ii).count];
        ccount = [ccount; count2+count3+count4];

    end
    

    %mean count for each flood type
    bstats(nn).count_mean = nanmean(allcount);
    bstats(nn).singles_count_mean = nanmean(count1);
    bstats(nn).doubles_count_mean = nanmean(count2);
    bstats(nn).triples_count_mean = nanmean(count3);
    bstats(nn).four_count_mean = nanmean(count4);
    bstats(nn).comp_count_mean = nanmean(ccount);

    %mean water level for each flood type (this is the mean of the flood
    %day water level means from each timeseries iteration)
    bstats(nn).singles_wl_mean = nanmean(mag1);
    bstats(nn).doubles_wl_mean = nanmean(mag2);
    bstats(nn).triples_wl_mean = nanmean(mag3);
    bstats(nn).four_wl_mean = nanmean(mag4);

    %std of counts for each flood type
    bstats(nn).count_std = std(allcount);
    bstats(nn).singles_count_std = nanstd(count1);
    bstats(nn).doubles_count_std = nanstd(count2);
    bstats(nn).triples_count_std = nanstd(count3);
    bstats(nn).four_count_std = nanstd(count4);
    bstats(nn).comp_count_std = nanstd(ccount);

    %std of water level for each flood type (this is the standard deviation
    %of the flood day water level magnitude from each timeseries iteration)
    bstats(nn).singles_wl_std = nanstd(mag1);
    bstats(nn).doubles_wl_std = nanstd(mag2);
    bstats(nn).triples_wl_std = nanstd(mag3);
    bstats(nn).four_wl_std = nanstd(mag4);




    %%%%%%%%%%%%%%%%% percentages of flood combinations

    %each column in the frac array represents a combination of flood
    %components set in the "bootstrap floods" function
    
    bstats(nn).frac_mean = zeros(1,width(boot(nn).frac)); 
    bstats(nn).frac_std = zeros(1, width(boot(nn).frac));
    
    
    for ii=1:width(boot(nn).frac);       
        bstats(nn).frac_mean(ii) = mean(boot(nn).frac(:,ii));
        bstats(nn).frac_std(ii) = std(boot(nn).frac(:,ii));
        bstats(nn).frac_lowCI(ii) = prctile(boot(nn).frac(:,ii),2.5);
        bstats(nn).frac_upCI(ii) = prctile(boot(nn).frac(:,ii), 97.5);
    end

end

%% Compile Flood Percentages
T_cols = [1,5,6,7,11,12,13,15]; %columns with T
NTOOE_cols = [3,6,9,10,12,13,14,15]; %columns with NTOOE
Q_cols = [2,5,8,10,11,13,14,15]; %columns with Q
R_cols = [4,7,8,9,11,12,14,15]; %columsn with Residual


for nn = 2:numel(boot);
    NTOOE_temp = []; %empty vector to save the percentages in
    Q_temp = [];
    T_temp =[];
    NTR_temp  = [];

    for ii = 1:numel(boot(nn).fd)
        NTOOE = boot(nn).frac(ii, NTOOE_cols); %number of NTOOE floods
        Q = boot(nn).frac(ii, Q_cols); %number of Q floods
        T = boot(nn).frac(ii,T_cols);
        NTR = boot(nn).frac(ii, R_cols);
        
        % get total floods for each component
        NTOOE_temp = [NTOOE_temp; sum(NTOOE)]; 
        Q_temp = [Q_temp; sum(Q)];
        T_temp = [T_temp; sum(T)];
        NTR_temp = [NTR_temp, sum(NTR)];
    end

    %preserve NTOOE percents
    bstats(nn).NTOOE_percent = mean(NTOOE_temp);
    bstats(nn).NTOOE_std = std(NTOOE_temp);
    bstats(nn).NTOOE_percent_lowCI = prctile(NTOOE_temp, 2.5);
    bstats(nn).NTOOE_percent_highCI = prctile(NTOOE_temp, 97.5);

    %preserve Q percents
    bstats(nn).Q_percent = mean(Q_temp);
    bstats(nn).Q_std = std(Q_temp);
    bstats(nn).Q_percent_lowCI = prctile(Q_temp, 2.5);
    bstats(nn).Q_percent_highCI = prctile(Q_temp, 97.5);

    %preserve T percents
    bstats(nn).T_percent = mean(T_temp);
    bstats(nn).T_std = std(T_temp);
    bstats(nn).T_percent_lowCI = prctile(T_temp, 2.5);
    bstats(nn).T_percent_highCI = prctile(T_temp, 97.5);

    %preserve NTR percents
    bstats(nn).NTR_percent = mean(NTR_temp);
    bstats(nn).NTR_std = std(NTR_temp);
    bstats(nn).NTR_percent_lowCI = prctile(NTR_temp, 2.5);
    bstats(nn).NTR_percent_highCI = prctile(NTR_temp, 97.5);
end



save bstats_data.mat bstats
clearvars -except data floods boot bstats start

%% timer
finish = toc(start);
disp(['time to run entire code is ', num2str(finish/60), ' minutes'])
