clear all;
close all;
tic

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

%% Download NOAA Data - Comment out once it's saved
% All sl/td data is in GMT time, metric units, and MHHW datum
% UTide must be in GMT to work 
% flood thresholdsin STND datum

for id=1:numel(noaaID);
    tic
        [datum]=noaaDatums(noaaID(id));
        [sl,dn]=noaaSealevel_GMT(noaaID(id),t0,tf);
        [flood]=noaaFlood(noaaID(id));
        [noaa_td, dn_local] = noaaTide(noaaID(id),t0,tf); %NOAA tide to compare to UTide

        %get UTide coefficients from sl
        coef = ut_solv(dn, sl, [], datum.lat, 'auto');
        %reconstruct
        [td, ~] = ut_reconstr(dn, coef);

        save(['noaa_tidegauge_',num2str(datum.id),'.mat'],'datum','sl','dn','td', 'flood', 'noaa_td', 'dn_local')
    toc

end 


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
    load(['noaa_tidegauge_',num2str(noaaID(nn)),'.mat']) %load data from saved mat files

    %convert data to daily 
    sl_daily = max(reshape(sl,24,numel(sl)/24,[],1)); %reshape sea level to be daily max
    td_daily = max(reshape(td,24,numel(td)/24,[],1)); %reshape tide to be daily max 
    NTOOE = sl_daily-td_daily; %NTOOE are maximum daily sea level - maximum daily tidal prediction


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
    data(1).q = usgsdata(1).Q; %q with seasonal cycle

    % NWS flood thresholds
    data(nn).thresh = flood.nws_minor - datum.MHHW ; %convert STND datum to MHHW 
    if nn==2 || nn==3 ;
        data(nn).thresh = prctile(data(nn).sl_daily,95); %no NWS thresholds for Brandywine and SJS
    end

    %Threshold by percentile
    data(nn).prct95 = prctile(data(nn).sl_daily,95); 

%% Multiple Regression 
%get observational results to compare to bootstrapping
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
    %a sanity check to remove the tide mean and add it back in the
    %classification step
    data(nn).res_mean = mean(residuals);
    data(nn).NTOOE_mean = mean(a*x1);
    data(nn).Q_mean = mean(b*x2);

    
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
    ind = find(data(nn).sl_daily >= data(nn).thresh);

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
clearvars -except data floods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin Bootstrapping

%% bootstrapping loop (comment out once iterations are compiled)
load data.mat

y1 = data(1).NTOOE;
y2 = data(1).q;
i = 1000;

tic
for nn=2:numel(data);
    x = data(nn).NTOOE;
    x2 = data(nn).sl_daily;
    x3 = data(nn).td_daily;
    thresh = data(nn).thresh;

    [boot(nn).a, boot(nn).b, boot(nn).r, boot(nn).r_sl, boot(nn).r_q, ...
        boot(nn).r_sum, boot(nn).residuals, boot(nn).fd, boot(nn).numfl,...
        boot(nn).percents] = bootstrap_floods(x,y1,y2,x2,x3,thresh,i);
end
toc

%save data and clear variables to clean up workspace
save bootstrap_data.mat boot 
clear all;

%% load saved bootstrap data
load data.mat %this dataset has only the observations for (variables "floods" and "data")

%this dataset has the 1000 bootstrap iterations of the data
%see "bootstrap floods" function 
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
   
    %means (all)
    %multiple regression
    bstats(nn).meana = nanmean(boot(nn).a); % a coeff mean
    bstats(nn).meanb = nanmean(boot(nn).b); % b coeff mean
    bstats(nn).meanr = nanmean(boot(nn).r); % rsq mean
    %single regressions
    bstats(nn).meanr_sl = nanmean(boot(nn).r_sl); % sl single regression
    bstats(nn).meanr_q = nanmean(boot(nn).r_q); % q single regression
    bstats(nn).meanr_sum = nanmean(boot(nn).r_sum); % single regression sum

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
        
        %flood days (this is the residual value on flood days, not the
        %residuals from the flood-only regression)
        fd_res_mean = [fd_res_mean; nanmean(boot(nn).fd(ii).res)];
        fd_res_std = [fd_res_std; nanstd(boot(nn).fd(ii).res)];
        %no autocorrelation in flood day residuals because unevenly spaced 
    end
    
    %compile mean/std of residuals
    %all
    bstats(nn).res_mean = mean(res_mean); %this is the mean of the means of the 1000 residual timeseries
    bstats(nn).res_std = mean(res_std); %this is the mean of the standard deviations of the 1000 residual timeseries
    %%%%% or is a better statistic the standard deviation of the 1000 residual timeseries means
    bstats(nn).res_acorr = mean(res_acorr); %this is the mean of the 1000 autocorrelations from the 1000 residual timeseries

    %flood days
    bstats(nn).fd_res_mean = mean(fd_res_mean); %this is the mean of the mean of the residuals on flood days from the 1000 ts
    bstats(nn).fd_res_std = mean(fd_res_std); %this is the mean of the std of the residuals on flood days from 1000 ts

   
    
    
    
    %%%%%%%%%%%%%%%%%% water level contributions on flood days 

    % Get the mean contribution for each of the 1000 iterations and
    % take the mean of those 1000 means to get the entire site mean
    % contribution. Take the std of those 1000 means to get the overall std

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

    for ii = 1:numel(boot(nn).fd);
        count1 = [count1; boot(nn).fd(ii).singles]; %number of single component floods in each of the iterations
        mag1 = [mag1; boot(nn).fd(ii).singles_wl]; %mean magnitude of the single component floods in that iteration

        count2 = [count2; boot(nn).fd(ii).doubles]; %number of 2 component floods in each of the iterations
        mag2 = [mag2; boot(nn).fd(ii).doubles_wl]; %mean magnitude of the 2 component floods in that iteration

        count3 = [count3; boot(nn).fd(ii).triples]; %number of 3 component floods in each of the iterations
        mag3 = [mag3; boot(nn).fd(ii).triples_wl]; %mean magnitude of the 3 component floods in that iteration

        count4 = [count4; boot(nn).fd(ii).four]; %number of 4 component floods in each of the iterations
        mag4 = [mag4; boot(nn).fd(ii).four_wl]; %mean magnitude of the 4 component floods in that iteration
    end

    %mean count for each flood type
    bstats(nn).singles_count_mean = nanmean(count1);
    bstats(nn).doubles_count_mean = nanmean(count2);
    bstats(nn).triples_count_mean = nanmean(count3);
    bstats(nn).four_count_mean = nanmean(count4);

    %mean water level for each flood type (this is the mean of the flood
    %day water level means from each timeseries iteration)
    bstats(nn).singles_wl_mean = nanmean(mag1);
    bstats(nn).doubles_wl_mean = nanmean(mag2);
    bstats(nn).triples_wl_mean = nanmean(mag3);
    bstats(nn).four_wl_mean = nanmean(mag4);

    %std of counts for each flood type
    bstats(nn).singles_count_std = nanstd(count1);
    bstats(nn).doubles_count_std = nanstd(count2);
    bstats(nn).triples_count_std = nanstd(count3);
    bstats(nn).four_count_std = nanstd(count4);

    %std of water level for each flood type (this is the standard deviation
    %of the flood day water level magnitude from each timeseries iteration)
    bstats(nn).singles_wl_std = nanstd(mag1);
    bstats(nn).doubles_wl_std = nanstd(mag2);
    bstats(nn).triples_wl_std = nanstd(mag3);
    bstats(nn).four_wl_std = nanstd(mag4);





    %%%%%%%%%%%%%%%%% combinations of flood components 

    %make empty arrays to fill in the next loop
    %each column in the percents/numfl arrays represents a combination of flood
    %components set in the "bootstrap floods" function

    %the final column of the percents matrix is the total of all the
    %percentages... so this will tell you how many floods could have been
    %caused by several different flood factor combinations
    
    bstats(nn).numfl_mean = zeros(1,width(boot(nn).numfl)); %raw number of floods in each timeseries iteration
    bstats(nn).numfl_std = zeros(1, width(boot(nn).numfl));
    
    bstats(nn).percents_mean = zeros(1,width(boot(nn).percents)); %percentage of floods observed in that iteration
    bstats(nn).percents_std = zeros(1,width(boot(nn).percents));
    
    for ii=1:width(boot(nn).percents);
        bstats(nn).percents_mean(ii) = mean(boot(nn).percents(:,ii));
        bstats(nn).percents_std(ii) = std(boot(nn).percents(:,ii));
        bstats(nn).percents_lowCI(ii) = prctile(boot(nn).percents(:,ii),2.5);
        bstats(nn).percents_upCI(ii) = prctile(boot(nn).percents(:,ii), 97.5);
        
        bstats(nn).numfl_mean(ii) = mean(boot(nn).numfl(:,ii));
        bstats(nn).numfl_std(ii) = std(boot(nn).numfl(:,ii));
        bstats(nn).numfl_lowCI(ii) = prctile(boot(nn).numfl(:,ii),2.5);
        bstats(nn).numfl_upCI(ii) = prctile(boot(nn).numfl(:,ii), 97.5);
    end

end

save bstats_data bstats
clearvars -except data floods boot bstats


toc