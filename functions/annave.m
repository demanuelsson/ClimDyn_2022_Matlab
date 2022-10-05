function [x, clim] = annave(sst, clim, plusorminus)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [x, clim] = annave(y);
% 
%    x = matrix with annual average removed
%    y = input matrix
%    clim = climatology
%    University of Washington online archive Atmospheric Science
%    their archive is no longer online
%%%%%%%%%%%%%%%%%


n = ndims(sst);

szsst = size(sst);
sst = reshape(sst, szsst(1), prod(szsst(2:n)));

[tmax,ngrid]=size(sst);

%  Get climatology
if nargin < 2; 
  plusorminus=-1;
  for m=1:12
    clim(m,:)=mean2(sst(m:12:(tmax),:));
    
    %clim(m,:)=nanmean(sst(m:12:(tmax),:));
    
    
  end
end;

%  Remove climatology from each time step
x=zeros(tmax,ngrid);
for m=1:tmax;
  l=rem(m-1,12)+1;
  x(m,:)=sst(m,:)+plusorminus*clim(l,:);
end

%  Put data back into output format
x = reshape(x, szsst);
clim = reshape(clim, [12 szsst(2:n)]);

