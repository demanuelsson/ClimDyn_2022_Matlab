%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  x=cosweight(sst, nya);
%   
%   This function weights data in matrix 'sst' by the cosine of 
%   latitudes in matrix 'nya'
%     University of Washington online archive Atmospheric Science
%     their archive is no longer online
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sst=era_z500_c5;
% nya=era_lat_c;


function x=cosweight(sst, nya)

ndim = sum(size(size(sst)))-1;
if ndim == 3
  [ntim, nlat, nlon] = size(sst);

% Assume the second dimension is y

  sst = reshape(shiftdim(sst, 1), nlat, nlon*ntim);
  [m,n] = size(nya);
  if m == 1; nya = nya'; end
  nya = nya .* (3.1415927)/180;
  sst = sst .* sqrt(abs(cos(nya * ones(1, nlon*ntim))));
  sst = shiftdim(reshape(sst, nlat, nlon, ntim), 2);

else

%  [tmax,ngrid]=size(sst);
%  ny=length(nya);
%  nx=ngrid/ny;
%  for i=1:ny
%    sst(:,(nx*(i-1)+1):(nx*i))=sst(:,(nx*(i-1)+1):(nx*i))*cos(nya(i)*pi/180);
%  end
  error('2D matrices not done correctly, use 3D');
end

x=sst;