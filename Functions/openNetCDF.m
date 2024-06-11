%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% McKeon and Piecuch 2024
% Compound high-tide floods and the role of discharge in the Delaware River
% Estuary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scripts to open NetCDF files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataFields = openNetCDF(fileName)

ncid=netcdf.open(fileName,'NC_NOWRITE');
[ndim, nvar, natt, unlim] = netcdf.inq(ncid);
for n=0:(nvar-1) 
 [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,n);
 eval(['dataFields.',varname,'=netcdf.getVar(ncid,n);'])
 for m=0:(natts-1)
  eval(['dataFields.',varname,'_',netcdf.inqAttName(ncid,n,m),...
      '=netcdf.getAtt(ncid,n,netcdf.inqAttName(ncid,n,m));'])
 end
end
netcdf.close(ncid);

return