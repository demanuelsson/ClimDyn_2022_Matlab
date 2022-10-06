function out = varimax(E,Lamda,L,A,norm)

% out = varimax(E,Lamda,L,A,norm)
%
% Varimax rotation of factor loadings.
%
% Input:
%      E = eigenvectors of unit length.
%  lamda = eigenvalues of the corresponding eigenvectors.
%      L  = number of EOFs to use in rotation (the program uses the L
%           largest PCs assuming the input loadings (eigenvectors)
%           are arranged in descending order).
%      A  = area weights (fractional weight of each "gridbox").
%   norm  = true (1) if loadings are to be normalized. In case of 
%           rotating loading vectors based upon correlation matrix.
%
% Output structure with elements:
%   rload	= rotated factor loadings.  
%   rexpv	= fractional variance explained by each loading vector.
%	 cscore	= score-coefficient matrix
%	 h			= communalities
%	 rownorm	= normalising factors for vectors - empty if norm false

%------- Originally from Yochanan Kushnir (May 1991) --------------

%disp(sprintf('Performing Factor Rotation \n'))
[M,N]=size(E);
var_total=sum(Lamda(1:L)/sum(Lamda));

% the loading vectors of unrotated EOFs, M by L
for k=1:L
   E(:,k)=sqrt(Lamda(k))*E(:,k);
end
Frot=E(:,1:L);

% Calculate communalities and normalize loadings:
for m=1:M, Var = Frot(m,:).^2; h(m,1) = sum(Var');, end;
fvar = sum(h .* A);
Lvec=[];
if norm
   %disp('Varimax: normalization enabled')
   for l=1:L; Frot(:,l) = Frot(:,l)./sqrt(h);, end;
   Lvec = sqrt(1./h);			% Bugfix - Sam Dean - 24-Mar-2011
end;

%----------------------------------------------------------------------
% Iterate to achieve the Varimax criterion until conversion:
nit = -1; ncm = 0; tvi = 0.;
while ncm < 6 
   % Calculate the varimax criterion:
   for l=1:L;
      sf1(l)=sum((Frot(:,l).^2).*A);
      sf2(l)=sum((Frot(:,l).^4).*A);
   end  
   tvf = sum(sf2.*M - sf1.^2)/(M*M);
   if nit >= 0;
      if abs(tvf-tvi) <= 0.000001; ncm = ncm +1; end;
   end
   nit = nit + 1; 	tvi = tvf;
   cstep=sprintf('Iteration No. %g Varimax criterion = %g',nit+1,tvf);
   
   % Rotate columns i and j:
   for i=1:L-1
      for j=i+1:L
         p = Frot(:,i);		q = Frot(:,j);       
         u = (p+q) .* (p-q);		v = 2*(p .* q);
         s = (u+v) .* (u-v);   	t = 2*(u .* v);
         a = sum(u .* A);		b = sum(v .* A);
         c = sum(s .* A); 		d = sum(t .* A);
         xn = d - 2*(a*b)/M;		xo = c-(a*a-b*b)/M;
         xr = sqrt(xn*xn + xo*xo);
         if xr > .001;
            cos4t = xo/xr;			% Define rotation angle
            cos2t = sqrt((1. + cos4t)/2.);
            cos1t = sqrt((1. + cos2t)/2.);
            sin1t = sqrt( 1. - cos1t^2  );
            if sin1t > .001;
               if xn < 0.; sin1t = -sin1t; end;
               Frot(:,i) = cos1t*p + sin1t*q;	% Rotate columns i & j
               Frot(:,j) = cos1t*q - sin1t*p;
            end;
         end;
      end;   % Go to next j
   end;     % Go to next i
end;       % Go for the next iteration until convergence  
%----------------------------------------------------------------------

% Un-normalize rotated loadings:
if norm
   for l=1:L; Frot(:,l) = Frot(:,l).*sqrt(h); end;
end;

% Recalculate communalities:
for m=1:M; Var = Frot(m,:).^2; h(m,1) = sum(Var'); end;
fvar = sum(h .* A);

% Calculate partial variance explained by each PC
for l=1:L; p = Frot(:,l) .^2; Var(l) = sum( p .* A ); end;
[pvar(L:-1:1),k] = sort(Var); Frot(:,L:-1:1) = Frot(:,k); 
for l=1:L; p = Frot(:,l) .^2; Var(l) = sum( p .* A ); end;

% the percentage of partial variance 
Var=Var/sum(Var)*100; 

% the percentage of explained variance in terms of original field
Var=Var*var_total;

% Calculating score-coefficient matrix by least-square method:
for l=1:L; Cscor(:,l) = Frot(:,l) .* sqrt(A); end;
Cscor = Frot*inv(Cscor'*Cscor);
for l=1:L; Cscor(:,l) = Cscor(:,l) .* sqrt(A); end;

disp(cstep)		% The last one
out = makestruct('rload',Frot, 'rexpv',Var, 'cscore',Cscor,'comm',h, 'rownorm',Lvec);
