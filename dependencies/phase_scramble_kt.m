%% READ FIRST
% x is NTOOE
% y is Q
% i is the number of times to bootstrap


%make sure all inputs are column vectors


function [kt, ktp]=phase_scramble_kt(x,y,i); 
kt=[];
ktp=[];

tic
    for ii=1:i;

        %only phase scramble where all datasets have data
        iii = find(~isnan(x) & ~isnan(y));
        x=x(iii); % NTOOE at tide gauge
        y=y(iii); %Q at Trenton

        pmat = [x,y];

        pmat2=phase_scramble_ts(pmat);

        [r, p] = corr(pmat2(:,1), pmat2(:,2), 'Type', 'Kendall');
        kt = [kt; r];
        ktp = [ktp; p];

    end
toc