%% READ FIRST
% x is the thing you want to estimate (NTOOE at point x upriver)
% y1 is the downstream boundary condition (NTOOE at Lewes)
% y2 is the upstream boundary condition (Q at Trenton)
% x2 is MHHW sea level at point x upriver
% x3 is MHHW predicted tide (daily max MHHW) at point x upriver
% thresh is the flood threshold 
% i is the number of times to bootstrap

%function outputs are 
% a coeff from multiple linear regression
% b coeff from multiple linear regression
% rsq from multiple linear regression
% nash sutcliffe efficiency test
% kling-gupta efficiency test
% rsq from NTOOE single linear regression
% rsq from Q single linear regression
% sum of rsqs from single linear regressions
% residuals from multiple linear regression
% flood days data structure
% number of floods of each type matrix
% percent of floods of each type matrix 


%% function

function [b_95, b_a, b_b, b_r, b_nse, b_kge, r_sl, r_q, r_sum, ... 
    residuals, fd, numfl, allfrac] = CHTF_bootstrapping(x, y1, y2, x2, x3, thresh, i); 
b_a=[];
b_b=[];
b_r=[];
b_nse =[];
b_kge = [];

b_95=[];

r_q = [];
r_sl = [];
r_sum = [];

tic
    for ii=1:i;
        %% make new timeseries
        rmax = length(x); %set the boundaries of the random index vector which will be 1: the length of the timeseries
        rmin = 1;
        ind = randi([rmin rmax],size(x)); %generate a random set of indexes that are the same length as timeseries 
        b_X = x(ind); %sample with replacement NTOOE at x
        b_L = y1(ind); %NTOOE lewes
        b_Q = y2(ind); %Q trenton
        b_SL = x2(ind); %sea level at x, not NTOOE
        b_td = x3(ind); %predicted tide at x 

        %standardize by removing nans and subtracting the mean 
        yy = b_X; 
        xx1 = b_L;
        xx2 = b_Q;

        %use only timesteps that have data in all three timeseries
        iii = find(~isnan(yy) & ~isnan(xx1) & ~isnan(xx2));  %get indexes where all three ts have values
        yy = yy(iii); xx1 = xx1(iii); xx2 = xx2(iii); %remake the timeseries (will be shorter)
        yy = yy - mean(yy); %remove means before regressing
        xx1 = xx1 - mean(xx1);
        xx2 = xx2 - mean(xx2);
        X = [xx1,xx2];

        %make full size vectors 
        XX1 = nan(size(b_X)); XX1(iii)=xx1;
        XX2 = nan(size(b_X)); XX2(iii)=xx2;
        
        %% solve
        mdl=fitlm(X,yy,'linear');
        a= mdl.Coefficients.Estimate(2);
        b= mdl.Coefficients.Estimate(3);
        r= mdl.Rsquared.Ordinary;
        
        %save coefficients
        b_a=[b_a;a];
        b_b=[b_b;b];
        b_r=[b_r;r];

        %% nash-sutcliffe
        Qo = yy;
        Qm = a*xx1 + b*xx2;

        n=[];
        d=[];
        for j = 1:numel(yy);
            n_temp = (Qo(j) - Qm(j))^2 ;
            d_temp = (Qo(j))^2;
    
            n = [n, n_temp];
            d = [d, d_temp];
        end
        nse = 1 - (sum(n)/sum(d));

        b_nse = [b_nse; nse];

        %% Kling-Gupta equation
        us = a*xx1 + b*xx2; %simulated
        uo = yy; %observed
    
        pr = corrcoef(us,uo);
        pr = pr(1,2);
        beta = mean(us,'omitnan') / mean(uo,'omitnan');
        alpha = var(us,'omitnan') / var(uo,'omitnan');

        kge = 1 - sqrt( ((pr-1)^2) + ((alpha-1)^2) );
        b_kge = [b_kge; kge];

        
        %% get residuals
        res = yy - a*xx1 - b*xx2; %get residuals in this step while they all have matches
        resnan = nan(size(b_X)); %make residual vector the same size as the data
        resnan(iii) = res; %add back in residuals where there are values, so now residuals are nans where one ts was missing values
        residuals(ii).res = res; %get individual timeseries of residuals to calculate later
        
        residuals(ii).mean = mean(res,'omitnan'); %this is a sanity check, means should be zero
        residuals(ii).std = std(res,'omitnan'); %standard deviation of residuals
        residuals(ii).acorr = r1auto(res); %get autocorrelation of residual timeseries
        
        
        %% single regressions
        
        %NTOOE
        sl_mdl = fitlm(xx1,yy,'Linear');
        r_sl = [r_sl; sl_mdl.Rsquared.Ordinary];
        
        %discharge
        q_mdl = fitlm(xx2,yy,'Linear');
        r_q = [r_q; q_mdl.Rsquared.Ordinary];

        %sum 
        tempsum = sl_mdl.Rsquared.Ordinary + q_mdl.Rsquared.Ordinary;
        r_sum = [r_sum; tempsum];
        
        %% get number and magnitudes of flood days in bootstrapped series 
        fd_ind = find(b_SL >= thresh); %find where total sea level is over thresh
        fd(ii).count = numel(fd_ind); %number of floods
        fd(ii).sl = b_SL(fd_ind); %magnitude of flood (MHHW sea level)
        fd(ii).q = b_Q(fd_ind); %discharge during flood days
        fd(ii).NTOOE_L = b_L(fd_ind); %Lewes NTOOE during flood days
        fd(ii).NTOOE_X = b_X(fd_ind); %Site x NTOOE during flood days
        fd(ii).td = b_td(fd_ind); %max daily predicted tide at point x upriver
        fd(ii).res = resnan(fd_ind); %residuals during floods
        fd(ii).cont_NTOOE = XX1(fd_ind) * a; %get NTOOE contribution on flood days using this iterations a coeff
        fd(ii).cont_Q = XX2(fd_ind) * b; %get discharge contribution on flood days using this iterations b coeff
        fd(ii).cont_td = b_td(fd_ind)-mean(b_td(iii)); %get tide contribution without mean tide included 
        % (remember iii is the index for only where data exists)
        % remove the tide mean because all the other variables have mean
        % excluded
  

    %% classify floods 

    %save flood day water level and tide
    %separate tide and MSL here so MSL isn't counted twice 
    floods = fd(ii).sl; %sl MHHW during flood
    tide = fd(ii).cont_td; %td w/o MSL during flood
    MSL = nanmean(b_SL); %mean sea level determined by bootstrapped sl timeseries

    %%%%%%% convert the boundary conditions to water contributions at the
    %%%%%%% site using the coefficients from the current iteration
   
    %single variables
    Q_days = fd(ii).q * b ; %FT1 
    SL_days = fd(ii).NTOOE_L * a; %FT2
    res_days = fd(ii).res; %FT3
    
    %two components
    %single variables with tides
    Q_tide = Q_days + tide;
    SL_tide = SL_days + tide;
    res_tide = res_days + tide;
    %single variables with residuals
    Q_res = Q_days + res_days;
    SL_res = SL_days + res_days;
    %other
    Q_SL = Q_days + SL_days;
    
    %three conditions
    Q_tide_res = Q_days + res_days + tide;
    SL_tide_res = SL_days + res_days + tide;
    Q_SL_tide = Q_days + SL_days + tide;
    Q_SL_res = Q_days + SL_days + res_days;
    
    %all four
    Q_SL_res_tide = Q_days + SL_days + res_days + tide;
        
    %%%%%%%%%%%%%%%%%%%% READ ABOUT THIS METHOD
    % we classify floods into 1,2,3, and 4 component floods
    % this method allows floods to be double counted within the component
    % class, but not across classes

    % frac matrix shows the percentage of all the flood types at that site
    % numfl matrix counts up the raw number of floods in each category
    
    % add in the MSL to all of the counts
    %% Tide Only - FT1
    tind = find((tide + MSL >= thresh)); 
    numfl(ii,1) = numel(tind);

    %% Discharge Only - FT2
    qind = find((Q_days + MSL >= thresh)); 
    numfl(ii,2) = numel(qind);

    %% NTOOE - FT3
    sind = find((SL_days + MSL >= thresh));
    numfl(ii,3) = numel(sind);

    %% Residual - FT4
    resind = find((res_days + MSL >= thresh));
    numfl(ii,4) = numel(resind);

     %% concatenate single component floods
    singles = [tind; qind; sind; resind];
    fd(ii).singles = numel(unique(singles)); %get number of single component floods
    fd(ii).singles_wl = mean(fd(ii).sl(singles)); %get mean water level during single component floods
    

    %% Singles Fraction
    %make empty matrix of size #components x #floods in this category
    %column 1 is tide
    %column 2 is Q
    %column 3 is NTOOE
    %column 4 is residual
    %column 5 is Q + T
    %column 6 is NTOOE + T
    %column 7 is R + T
    %column 8 is Q + R
    %column 9 is NTOOE + R
    %column 10 is Q + NTOOE
    %column 11 is Q + T + R
    %column 12 is NTOOE + T + R
    %column 13 is Q + NTOOE + T
    %column 14 is Q + NTOOE + R
    %column 15 is Q + NTOOE + T + R

    frac = zeros(numel(floods), 15);

    singlesu = unique(singles);

    for dd = 1:numel(singlesu);
        dn_ind = singlesu(dd);
        
        % tide only
        for tt = 1:numel(tind);
            if tind(tt) == dn_ind;
                frac(dn_ind,1) = 1/(sum(dn_ind==singles))/numel(floods);
            end
        end

        % q only
        for qq = 1:numel(qind);
            if qind(qq) == dn_ind;
                frac(dn_ind,2) = 1/(sum(dn_ind==singles))/numel(floods);
            end
        end


        % NTOOE only
        for ss = 1:numel(sind);
            if sind(ss) == dn_ind;
                frac(dn_ind,3) = 1/(sum(dn_ind==singles))/numel(floods);
            end
        end        

        % res only
        for rr = 1:numel(resind);
            if resind(rr) == dn_ind;
                frac(dn_ind,4) = 1/(sum(dn_ind==singles))/numel(floods);
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Two components 


    %% Discharge + Tide -- FT5
    qtind = find((Q_tide + MSL >= thresh));
    qtind = setdiff(qtind, singles); %exclude any floods that have already been classified under single component 
    numfl(ii,5) = numel(qtind);

    %% NTOOE + Tide -- FT6
    stind = find((SL_tide + MSL >= thresh));
    stind = setdiff(stind,singles);
    numfl(ii,6) = numel(stind);

    %% Residual + Tide -- FT7
    restind = find((res_tide + MSL >= thresh));
    restind = setdiff(restind,singles);
    numfl(ii,7) = numel(restind);

     %% Discharge + Residual -- FT8
    qrind = find((Q_res + MSL >= thresh));
    qrind = setdiff(qrind,singles);
    numfl(ii,8) = numel(qrind);

    %% NTOOE + Residual -- FT9
    srind = find((SL_res + MSL >= thresh));
    srind = setdiff(srind,singles);
    numfl(ii,9) = numel(srind);

    %% Q + NTOOE -- FT10
    qsind = find((Q_SL + MSL >= thresh));
    qsind = setdiff(qsind,singles);
    numfl(ii,10) = numel(qsind);

    %% concatenate 2 component floods

    doubles = [qtind;stind;restind;qrind;srind;qsind];
    fd(ii).doubles = numel(unique(doubles)); %get number of days of 2 component floods, without double counting them
    fd(ii).doubles_wl = mean(fd(ii).sl(doubles)); %get mean water level during 2 component floods


    %% Doubles fraction
    %column 1 is tide
    %column 3 is NTOOE
    %column 2 is q
    %column 4 is residual
    %column 5 is Q + T
    %column 6 is NTOOE + T
    %column 7 is R + T
    %column 8 is Q + R
    %column 9 is NTOOE + R
    %column 10 is Q + NTOOE
    %column 11 is Q + T + R
    %column 12 is NTOOE + T + R
    %column 13 is Q + NTOOE + T
    %column 14 is Q + NTOOE + R
    %column 15 is Q + NTOOE + T + R

    doublesu = unique(doubles);

    for dd = 1:numel(doublesu);
        dn_ind = doublesu(dd);
        
        % Q + Tide floods
        for qt = 1:numel(qtind);
            if qtind(qt) == dn_ind
                frac(dn_ind,5) = 1/(sum(dn_ind==doubles))/numel(floods);
            end
        end

        % NTOOE + Tide floods
        for st = 1:numel(stind);
            if stind(st) == dn_ind;
                frac(dn_ind,6) = 1/(sum(dn_ind==doubles))/numel(floods);
            end
        end

        % Res + Tide floods
        for rt = 1:numel(restind);
            if restind(rt) == dn_ind;
                frac(dn_ind,7) = 1/(sum(dn_ind==doubles))/numel(floods);
            end
        end

        % Q + Res floods
        for qr = 1:numel(qrind);
            if qrind(qr) == dn_ind;
                frac(dn_ind,8) = 1/(sum(dn_ind==doubles))/numel(floods);
            end
        end
        
        % NTOOE + Res floods
        for sr = 1:numel(srind);
            if srind(sr) == dn_ind;
                frac(dn_ind,9) = 1/(sum(dn_ind==doubles))/numel(floods);
            end
        end

        % NTOOE + Q floods
        for qs = 1:numel(qsind);
            if qsind(qs) == dn_ind;
                frac(dn_ind,10) = 1/(sum(dn_ind==doubles))/numel(floods);
            end
        end


    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Three components
    
    %% Q + tide + residual
    qtrind = find((Q_tide_res + MSL >= thresh));
    qtrind = setdiff(qtrind,doubles); %exclude floods already classified as doubles
    qtrind = setdiff(qtrind,singles); %exclude floods already classified as singles
    numfl(ii,11) = numel(qtrind);

    %% Sugre + Tide + Residual
    strind = find((SL_tide_res + MSL >= thresh));
    strind = setdiff(strind,doubles);
    strind = setdiff(strind,singles);
    numfl(ii,12) = numel(strind);
    
    %% Q + NTOOE + Tide
    qstind = find((Q_SL_tide + MSL >= thresh));
    qstind = setdiff(qstind, doubles);
    qstind = setdiff(qstind,singles);
    numfl(ii,13) = numel(qstind);
    
    %% Q + NTOOE + Residual
    qsrind = find((Q_SL_res + MSL >= thresh));
    qsrind = setdiff(qsrind,doubles);
    qsrind = setdiff(qsrind,singles);
    numfl(ii,14) = numel(qsrind);

    %% concatenate triples 
    triples = [qtrind;strind;qstind;qsrind];
    fd(ii).triples = numel(unique(triples)); %get number 
    fd(ii).triples_wl= mean(fd(ii).sl(triples)); %get water level

    %% Triples fraction
    %column 1 is tide
    %column 2 is NTOOE
    %column 3 is q
    %column 4 is residual
    %column 5 is Q + T
    %column 6 is NTOOE + T
    %column 7 is R + T
    %column 8 is Q + R
    %column 9 is NTOOE + R
    %column 10 is Q + NTOOE
    %column 11 is Q + T + R
    %column 12 is NTOOE + T + R
    %column 13 is Q + NTOOE + T
    %column 14 is Q + NTOOE + R
    %column 15 is Q + NTOOE + T + R

    triplesu = unique(triples);

    for dd = 1:numel(triplesu);
        dn_ind = triplesu(dd);
        
        % Q + Tide + Residual floods
        for qtr = 1:numel(qtrind);
            if qtrind(qtr) == dn_ind;
                frac(dn_ind,11) = 1/(sum(dn_ind==triples))/numel(floods);
            end
        end

        % NTOOE + Tide +Residual floods
        for str = 1:numel(strind);
            if strind(str) == dn_ind;
                frac(dn_ind,12) = 1/(sum(dn_ind==triples))/numel(floods);
            end
        end

        % Q + NTOOE + Tide floods
        for qst = 1:numel(qstind);
            if qstind(qst) == dn_ind;
                frac(dn_ind,13) = 1/(sum(dn_ind==triples))/numel(floods);
            end
        end

         % Q + NTOOE + REsidual floods
        for qsr = 1:numel(qsrind);
            if qsrind(qsr) == dn_ind;
                frac(dn_ind,14) = 1/(sum(dn_ind==triples))/numel(floods);
            end
        end       


    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Four Fraction

    qstrind = find((Q_SL_res_tide + MSL >= thresh)); 
    qstrind = setdiff(qstrind,singles);
    qstrind = setdiff(qstrind,doubles);
    qstrind = setdiff(qstrind,triples); %exclude all the relevant preclassifications

    numfl(ii,15) = numel(qstrind);
    fd(ii).four = numel(qstrind);
    fd(ii).four_wl = mean(fd(ii).sl(qstrind));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Assign Fraction
    %column 1 is tide
    %column 2 is NTOOE
    %column 3 is q
    %column 4 is residual
    %column 5 is Q + T
    %column 6 is NTOOE + T
    %column 7 is R + T
    %column 8 is Q + R
    %column 9 is NTOOE + R
    %column 10 is Q + NTOOE
    %column 11 is Q + T + R
    %column 12 is NTOOE + T + R
    %column 13 is Q + NTOOE + T
    %column 14 is Q + NTOOE + R
    %column 15 is Q + NTOOE + T + R
    

    for dd = 1:numel(qstrind);
        dn_ind = qstrind(dd);
        frac(dn_ind,15) = 1/1/numel(floods); %no repeat combinations here
    end
    
    allfrac(ii,:) = sum(frac)*100;

    end
toc
end
