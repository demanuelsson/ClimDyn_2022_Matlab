function [xkeep, ykeep] = keep_var(lim, x, y);
%%%%%%%%%%%%%%%%%%%%%%%
%    [xkeep, ykeep] = keep_var(lim, x, y);
% 
%    where lim = [minx maxx miny maxy] to be kept.
% 
%    If lon and lat are not input, then they are assumed
%    to be global variables under the names XAX and YAX.
%    University of Washington online archive Atmospheric Science
%    their archive is no longer online
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin == 1;
  global XAX YAX
  x = XAX;
  y = YAX;
end

if (lim(2) <= lim(1) | lim(4) <= lim(3))
  error(['lim must be input as [minx maxx miny maxy];  '...
         'try keep_var2'])
end

xkeep = find(x >= lim(1)-eps & x <= lim(2)+eps);
ykeep = find(y >= lim(3)-eps & y <= lim(4)+eps);