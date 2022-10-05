function out = eof_j(ldata,nsave,scaling,atype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = eof(ldata,nsave,scaling,atype)
%
% Get eigenvectors/values by eig or SVD method
% Input:
%	ldata	-	n by p data array
%	nsave	-	How many vectors (etc) to save/return - 0 => all
%	scaling -	if true (1), scale columns of ldata by std(ldata) first
%	atype	-	analysis type: 'eig' or 'svd' (default depends on n/p)
%
% Output:
%  Structure with components
%	vect		-	Eigenvectors (Most important "nsave")
%	val		-	Eigenvalues (  ..     ..       ..   )
%	series	-	Principal components (loadings, time series)
%	regmap	-	Regression map version of vectors
%	totvar	-	Total variance of input data
%	expv		-	Percent explained variance
%	rank		-	Rank info: [sum(val>1e-7) spatial d.o.f]
%	std		-	standard deviations of columns of ldata (=1 if scaling)
%     [from James Renwick's library]
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,p]=size(ldata);

if nargin < 4 | isempty(atype)
	if n./p <= 2./3, atype='svd'; else atype='eig'; end
end
if nargin < 3 | isempty(scaling), scaling=0; end
if nargin < 2 | isempty(nsave), nsave=0; end

[n,p]=size(ldata);
msave=nsave;
if nsave<=0, msave=p; end
atype=lower(atype);
disp(sprintf('Performing %s EOF calculations, saving %g',atype,msave))
if scaling, disp('EOFs calculated for normalised data'), end

%   Get anomalies

lmean=mean(ldata);
% ldata=make_anoms(ldata,scaling);
lstd=std(ldata);

%   Get eigen-form

flipped=false;
if strcmp(atype,'svd')
	[u,s,v]=svd(ldata,0);			% Go for "economy size"
	%disp('Have singular decomposition')
	
	val=(diag(s).^2)'/n;
	totvar=sum(val);
	expv=100*val(1:msave)/totvar;
	vect=v(:,1:msave);
	series=u(:,1:msave)*s(1:msave,1:msave);
	clear u s
else
	if p./n >= 2
		disp(sprintf('Transposing %ix%i data matrix to save space on covariance matrix',n,p))
		ldata=ldata';			% save memory
		flipped=true;
	end
	llcov=cov(ldata,1);
	%disp('Have covariance matrix')
	
	[v,d]=eig(llcov);		%   Get eigen-form
	%disp('Have eigen-decomposition')
	
	[v,eval]=sorteig(v,d);
 	clear d   
	z=sum(eval);
	expv=100*eval(1:msave)/z;
	vect=v(:,1:msave);
	val=eval;					% Return all of them
	totvar=sum(val);
	series = ldata * vect;
	if flipped
		zs=series;
		series=vect;
		vect=zs;
		clear zs
		series=scale(series);
		vect=normalise(vect);
	end
end

rstats = [sum(val > 1.0e-7) (sum(val).^2)/sum(val.^2)];

z=expv(1:min(6,msave));
z=round([z sum(z)]);
nz=length(z);
disp(sprintf(['E.V: ',repmat('%i+',1,nz-2),'%i=%i  (Rank|s.DoF: %i,%i)'],z,round(rstats)))

%   Finally, covariance patterns, and give it all back

%disp(sprintf('Making covar patterns'))
regmaps=vect * diag(sqrt(val(1:msave)));

out=makestruct('vect',vect,'val',val,'series',series,'regmap',regmaps, ...
	'totvar',totvar,'expv',expv,'rank',rstats,'mean',lmean,'std',lstd);
