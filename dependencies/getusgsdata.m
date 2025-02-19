%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% McKeon and Piecuch 2025
% Compound minor floods and the role of discharge in the Delaware River
% Estuary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to download USGS discharge data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data]=getusgsdata(site,parameter,dstart,dend,varname)
warning('off')

front='https://waterdata.usgs.gov/nwis/dv?';%start of usgs url data portal
partxt=[];
for i=1:length(parameter); %loop for parameter string
    partxt=[partxt, 'cb_', parameter{i}, '=on&'];
end

url=strcat(front,partxt,'format=rdb&site_no=',site, ...
    '&period=&begin_date=', dstart,'&end_date=',dend)

options=weboptions('ContentReader',@readtable);
options.Timeout = 10 
data=webread(url,options);%downloading data using url

%%parcing variables from downloaded table
dn=datenum(table2array(data(2:height(data),'datetime')));% date number string
dname=data.Properties.VariableNames;%getting variable names in table

%loop to find column for each parameter number specified
for i=1:length(parameter);
    ind=find(contains(dname,parameter{i}));
    %internal loop for cases when 2 columns of data for parameter
    %(e.g. top and bottom observations)
    for k=1:length(ind)/2;
        eval([varname{i} '(:,k)=str2double(table2array(data(2:height(data),dname{ind(k*2-1)})));']);
    end
end


end


