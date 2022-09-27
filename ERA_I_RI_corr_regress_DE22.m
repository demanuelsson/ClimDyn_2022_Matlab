%%%%%%%%%%%%%%%%%%%%%%%%%
% ERA_I_RI_corr_regress_DE21.m
% Daniel Emanuelsson, 2021
% MATLAB 2020b
%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
% * corrcoef_df.m   UoW Steig
% * cosweight.m     UoW online archive Atmospheric Science
% * annave.m        UoW online archive Atmospheric Science
% * star_coord_WA.m 
% * fig.m fileexchange
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  original file era_wais_regression_test.m
%  Wilks test
tic
clear
close all

site='RICE';
cor_nr=14;  % (14 dD d18O RICE)
%%%%%%%%%%%%%%%%%
      
if cor_nr==14
   %iso='dD';
   iso='d18O';
end


%%%%%%%%%%%%%%%%%%%%%%%
crop_nr=0; % (1/0) crop?
lock_scalebar=1; % (1/0)
%rcontour_p=0;% (1/0) SAT reg contour

% corr_time=1; % (1) annual, (2) monthly
% show_max=0;   % (1/0) show max label
 max45=1;


%%%%%%%%%%Figure format
% figure_format=2; % (1) - EPS, (2) - PNG 
%%%%%%%%%%

corr_label=1; %(1/0) %%%%%%%%%%% Show colorbar and label
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_value=1; % sign (1)  positive (2) negative r-values in the ABS/Ross Sea region

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
era_name_nr=2; % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Z500 *
% 2. 2mT *

if era_name_nr==1
    rcontour_pc1=1; %%%%%%%%%% (1/0) on/off Show regression contours and stippling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    rcontour_pc1=0;
end

if era_name_nr==1
name='z500';
name_str='Z500';
elseif era_name_nr==2
name='2mT'; 
name_str='SAT';
end

yr_s=1979;
%yr_e=2009;
yr_e=2011;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sea_nr=3;  %%%%%%%%%%          Season
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sea_nr==1
    season='annual';
elseif sea_nr==2 
    season='MAMJJASON';
elseif sea_nr==3    
    season='AMJJASON';
end
proj='stereo';
%proj='mercator';
% lat1=-90;
% lat2=-40; 



%        Size:       480x241x421
%        Dimensions: longitude,latitude,time
 addpath C:\Users\benman\matlab_lib\Data\ERA_interim
% ncdisp('ERA_int_monthly_z500.nc') 
% ncdisp('ERA_int_monthly_2m_T.nc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load ERA-interim data

name_c='ERA_int_monthly_z500_2.nc';

era_time=ncread(name_c,'time');
era_long=ncread(name_c,'longitude');
era_lat=ncread(name_c,'latitude');

% time
dayssince111=era_time/24;
datevalue=dayssince111+datenum(1900,1,1);
 date_vec=datevec(double(datevalue)); 
 yyyy=date_vec(:,1); 
 mm=date_vec(:,2);
 era_year_num=yyyy+(mm/12-1/12); 
%

if strcmp(name,'2mT')==1
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of daily means
    % surface 2mT 
    name_c='ERA_int_monthly_2m_T_c.nc';
    era_T=ncread(name_c,'t2m'); % 2m temp
    era_T=era_T- 273.15;
    %letter='c';
    letter='d';
    lat1=-90; 
    lat2=-30; 
    if yr_e==2009
    letter='b';    
    end
 
 
 elseif strcmp(name,'z500')==1 %&& strcmp(season,'annual')==1
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of dail means
    % geopotentail 500 hPa  
     
     
    name_c='ERA_int_monthly_z500_2.nc'; 
    era_z500=ncread(name_c,'z'); 
    era_z500=era_z500/9.80665;
       letter='c';
    %letter='b'; 
    if yr_e==2009
    letter='a';    
    end
  
        if strcmp( proj,'stereo')==1
        lat1=-90;
        lat2=-10;
        end
 end

% era_count=396; % that overlaps with RICE record 1979-2011
era_count=find(era_year_num==yr_e)+11;  % data in months
%era_start=find(era_year_num==1980-1);

% annual 
era_start=find(era_year_num==yr_s);
div=ceil((era_count-era_start)/12);

    if strcmp(name,'2mT')==1
       era_z500=era_T;

    end
toc
%%%%
%% resize grid
 tic
 reanalysis_nr=1;
 if reanalysis_nr==1 || reanalysis_nr==4 % 
 
 [a b c]=size(era_z500);
 
 %scale=0.5;
 scale=0.75;
 %scale=1;
 
    for i=1:c
        
    era_z500_v(:,:,i)=resizem(era_z500(:,:,i),scale);
    
    end 
 
    
  [av bv cv]=size(era_z500_v);
%     era_lat_v=resizem(era_lat,scale);
%     era_long_v=resizem(era_long,scale);  


% xq = (1:bv)*scale + 0.5 * (1 - scale);

xq = (1:bv)*1/scale; % Lat

yq = (1:av)*1/scale; % Lon

% resize lat long

x=[1:length(era_lat)]';  % Lat
[p,s]=polyfit(x,era_lat,1);

era_lat_v=polyval(p,xq');
clear x p s

x=[1:length(era_long)]';
[p,s]=polyfit(x,era_long,1);


era_lon_v=polyval(p,yq');    
%     clear a b
    
    era_z500=era_z500_v;
    era_lat=era_lat_v;
    era_long=era_lon_v;
 
 downscale_str=['downscale_',num2str(scale)];   
    
 end
 toc
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%
 tic

 %   M=era_z500;
 M=era_z500(:,(91:end),:); % SH 
 era_lat_c=era_lat;
[nlon, nlat, ntim]= size(M);
era_z500_c = permute(M, [3 2 1]);


%%% -=(1)=-
% Remove monthly climatology for all cells 

 era_z500_c2 = double(reshape(era_z500_c, ntim, nlat*nlon));  % One time series (column) for each grid point
 
 %  lat 1, 2, 3, 4,..........33...lat 1, 2, 3, 4.....33       until 33X240=7920
 %  
 %1
 %time
 
 [era_z500_c3,clim_z500] = annave(era_z500_c2);   % Removes seasonal cycle

 era_z500_c4= reshape(era_z500_c3,  ntim, nlat, nlon );

 
%%% -=(2)=- 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Since the grid size decreases as you move towards the pole, 
% weight each grid box (i.e., multiply the time series at each grid box) 
% by the square root of the cosine of latitude (the weights are based on 
% the square root of the cosine so that the covariance matrix is weighted by the cosine of latitude).

%  z = cosweight(z, lat);
% weighted to account for the decrease in area towards the pole.
%era_z500_c2=cosweight(era_z500_c,era_lat_c); % Original UoW function    (time x lat x lon) 
era_z500_c5=cosweight(era_z500_c4,era_lat_c(91:end));
%%%%%%%%%%%

% Back to old format again
era_z500_c6= reshape(era_z500_c5,  ntim, nlat, nlon );
M = permute(era_z500_c6, [3 2 1]);  %  

% S_lon=266;S_lat=40;%check
% era_z500_c4 ( S_lon, S_lat-20,2) 
%  era_z500_c4=cosweight(era_z500_c4,era_lat(91:end)); % UoW function
%  era_z500_c4 ( S_lon, S_lat-20,2) % check should changes value before and after

%%% Annual and seasonal means

[m n t]=size(M(:,:,1:era_count)); %
for i= 1:m
    for j= 1:n

        dummy1=reshape(squeeze(M(i,j,era_start:era_count)), 12, div);          
%           
     % Annual    
    if  strcmp(season,'annual')==1
         
%          if corr_time==1
          ERA_M_annual=mean(dummy1,1); % annual ERA-values
          ERA_M_annual_reg(i,j,:)=ERA_M_annual(1:end); % save annual values for regression
          ERA_M_annual_dummy=detrend(ERA_M_annual(1:end)); % change to 1 to get 1979 annual
          ERA_M_annual_detrend(i,j,:)=ERA_M_annual_dummy;            
%          end

    elseif strcmp(season,'MAMJJASON')==1
        
          ERA_M_season_mamjjason=nanmean(dummy1(3:11,:),1);  %3. seasonal (JJA)
          ERA_M_season_mamjjason=ERA_M_season_mamjjason(1:end);
          ERA_M_season_reg_mamjjason(i,j,:)= ERA_M_season_mamjjason; % for regression
          ERA_M_season_dummy_mamjjason=detrend(ERA_M_season_mamjjason); 
          ERA_M_season_detrend_mamjjason(i,j,:)=ERA_M_season_dummy_mamjjason;  
          
    elseif strcmp(season,'AMJJASON')==1
        
          ERA_M_season_amjjason=nanmean(dummy1(4:11,:),1);  %3. seasonal (JJA)
          ERA_M_season_amjjason=ERA_M_season_amjjason(1:end);
          ERA_M_season_reg_amjjason(i,j,:)= ERA_M_season_amjjason; % for regression
          ERA_M_season_dummy_amjjason=detrend(ERA_M_season_amjjason); 
          ERA_M_season_detrend_amjjason(i,j,:)=ERA_M_season_dummy_amjjason;         
        
          
    end  
    
    end
end



% %%
% % save anomalies
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save_anomalies=1;
% 
% if save_anomalies==1
%     scale=0.5;
% 
%     filedir ='C:\PHD\matlab_storage_of_output_files\';
% 
%     name_str=['ERAinterim_downscaled_',name,'_scale_',num2str(scale),'_2020_anomalies_',num2str(yr_s),'-',num2str(yr_e)];
%  
%     savefilename =[filedir, name_str '.mat'];
% 
%    
%     save(savefilename,'M','era_lat','era_long','era_year_num');
%     
%         
% end

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% isotopes record
tic
if strcmp(iso,'dD')==1 || strcmp(iso,'d18O')==1 
% if corr_time==1    % Annual 
%load('C:\PHD\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c19.mat'); % May 2017
%load('C:\PHD\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c23.mat');

    iso_alt_nr=3;
        
        if iso_alt_nr==1
             load('C:\Users\benman\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c23.mat');
             iso_season_str='annual';
        elseif iso_alt_nr==2
            load('C:\Users\benman\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_c24.mat'); % MAMJJASON 
%             load('C:\Users\benman\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_clim_c24.mat'); % clim removed, no difference
            iso_season_str='nosummer';
        elseif iso_alt_nr==3    
            load('C:\Users\benman\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_AMJJASON_c24.mat'); % AMJJASON 
%             load('C:\Users\benman\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_AMJJASON_clim_c24.mat'); % clim removed, no difference
            iso_season_str='nosummer';           
        end

stacked_record_annual_Ma=MA_save;
date_annual=stacked_record_annual_Ma(:,1);
    
% if  strcmp(season,'annual')==1 && corr_time==1
   start_t=find(date_annual==yr_s);
   end_t=find(date_annual==yr_e); 
   

% if corr_time==1    % Annual 
  %  X=detrend(stacked_record_annual_Ma((start_t:end_t),2)); 
    X=detrend(stacked_record_annual_Ma((start_t:end_t),3));% d18O 
end
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERA-data choose annual or season 

if strcmp(season,'annual')==1
%    season='Annual';
   Y=ERA_M_annual_detrend;
   
elseif strcmp(season,'MAMJJASON')==1
    
   Y=ERA_M_season_detrend_mamjjason;
   
elseif strcmp(season,'AMJJASON')==1
    
   Y=ERA_M_season_detrend_amjjason;   
   
end
  
[A B C]=size(Y);
R=zeros(2,2,A,B);P=zeros(2,2,A,B);
for i=1:A
    for j=1:B
        Q=squeeze(Y(i,j,:));
      [R(:,:,i,j) P(:,:,i,j)]=corrcoef_df(X,Q, 'rows','pairwise');
     % moving average
     %   [R(:,:,i,j) P(:,:,i,j)]=corrcoef_df(moving_average(X,5),moving_average(Q,5), 'rows','pairwise');
    end
end

% for RI
% R(:,:,find(floor(era_long)==199),find(floor(era_lat)==-79))
% P(:,:,find(floor(era_long)==199),find(floor(era_lat)==-79))

toc 

%%
 tic
 
% Wilks test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contour_alt_nr=2; % (1) p regular local, (2) with Wilks 2006, 2016 test

if contour_alt_nr==2
 

p_in= squeeze(P(1,2,:,:));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alt 2
% from Github 
%find the locations where the null hypothesis 
%could be rejected locally
% alp_c=.05;
alp_c=.1;
rej_null = find(p_in<alp_c); % changes from grid to vector

no_rej_null = length(rej_null); 
%sort from smallest to largest the p values that
%would suggest the null hypothesis could be rejected
%locally
[rej_null_sort, rej_null_ind] = sort(p_in(rej_null)); 
[rn_x,rn_y] = ind2sub([360,90],rej_null(rej_null_ind));
p_local_rej = ones(360,90);
p_fdr_rej = ones(360,90);
no_grid_points = 360*90;
n_fdr = 0; 

%loop over all the cases where it the null hypothesis
%is assumed could be rejected locally
for i = 1:no_rej_null
    pi_c = p_in(rej_null(rej_null_ind(i)));
    p_local_rej(rn_x(i),rn_y(i))= pi_c;
    % compute fdr threshold
    
    
   % fdr_thres = 2*alp_c*i/no_grid_points; % correct?
   % fdr_thres = 2*alp_c*((no_grid_points+1)-rej_null(rej_null_ind(i)))/no_grid_points; % modified otherwise the sort/rank isn't taken into acount
   % fdr_thres = 2*alp_c*((no_rej_null+1)-rej_null_ind(i))/no_grid_points; % modified otherwise the sort/rank isn't taken into acount
    fdr_thres = alp_c*((no_rej_null+1)-rej_null_ind(i))/no_grid_points; % remove 2, eq. 3 Wilks 2016
    
    
    
    %find local p values that are less than fdr threshold
    if(pi_c<fdr_thres)
        p_fdr_rej(rn_x(i),rn_y(i)) = pi_c;
        n_fdr = n_fdr+1;
    end 
end

end
toc

%% Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
% axesm define projection
tic
f1=fig('units','inches','width',10.5,'height',10.5,'font','Helvetica','fontsize',18,'border','on'); 

proj='stereo';
%proj='mercator';

axesm( proj,'MapLatLimit',[lat1 lat2],'Grid','on','ParallelLabel','on','Frame','on',... %
       'MeridianLabel','on','FontWeight','bold','FontSize',18,...
       'mlabellocation',[0:30:179,0:-30:-180]); 
   %Frame - around figure
   
set(gca,'box','off','XColor',[1,1,1],'YColor',[1,1,1]);  % turns off border box and axes  
  
gridm('GLineStyle','--','Gcolor',[.6 .6 .6],'GLineWidth',1.5,...
    'MLineLimit',[lat1 lat2],'Galtitude', .02)

hold on

%%%%%%%%%%%%%%%%%%%%
c_min=squeeze(R(1,:,:));
c_min=c_min(2,:);
c_min=min(c_min);
c_max=squeeze(R(1,:,:)); 
c_max=c_max(2,:);
c_max=max(c_max);

c_max_c=c_max;

% Masking out areas that are not statistically significant
Rc=R;

% Rc2=R;

% find lat long max 
% c_max=max(c_max,[],2);
% find(c_max==c_max_c)   (54-lat,
%c_max=X_r(2,:);


X_r=squeeze(Rc(1,2,:,:))';
% X_r=X_r((181:241),:); % from -45S to -90 max

if max45==1
%         X_r((1:181),:)=NaN;  %%% MAX below -45 degrees, and east ant E 0-150
%         X_r(:,(1:201))=NaN; 
        X_r((1:46),:)=NaN;  %%% MAX below -45 degrees, and east ant E 0-150
%         X_r(:,(1:201))=NaN; 
        
else
        X_r=X_r;  %%% MAX below -45 degrees, and east ant E 0-150
%         X_r(:,(1:201))=NaN;   
end

era_lat_c2=era_lat(91:end);

if r_value==1
X_r_max=nanmax(nanmax(X_r));
X_r_long=nanmax(X_r); % 1 long
X_r_lat=nanmax(X_r,[],2); %2

elseif r_value==2
X_r_max=nanmin(nanmin(X_r));
X_r_long=nanmin(X_r); % 1 long
X_r_lat=nanmin(X_r,[],2); %2
end
    
xr=find(X_r_long==X_r_max); 
xr_lon=era_long(xr);

xr=find(X_r_lat==X_r_max);
    if max45==1
    xr_lat=era_lat_c2(xr);
    else
    xr_lat=era_lat_c2(xr);
    end
c_limit=max(abs(c_min),abs(c_max));
c_limit_c=c_limit;

% significant level
% pp_wr_c=p_fdr_rej'
% pp_wr_cx=pp_wr_c(61,227);
% p_level(pp_wr_cx)

% if lock_scalebar==1 && strcmp(name,'2mT') && (strcmp(season,'annual')) && (strcmp(iso,'dD'))
%     c_limit=0.7484;
% %    c_limit=1;
% elseif lock_scalebar==1 && strcmp(season,'annual') && strcmp(iso,'dD')    
%     c_limit=0.7351;   
% end
if lock_scalebar==1
  c_limit=0.701;  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rc(P>=0.1)=NaN; % maksing comment out if you dont want it
Rc_c=squeeze(Rc(1,2,:,:))';
 
% %%%%%%%%%%%%% fill in seam 
% era_long_c=era_long;
% era_lat_c=era_lat;
% era_long_c(481,:)=359.99;
% Rc_c(:,481)=Rc_c(:,1);
% %%%%%%%%%%%%%%%%%%



hSurf=surfm(double(era_lat_c2),double(era_long),squeeze(Rc_c));
colormap(b2r(-c_limit,c_limit));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(name,'2mT') && (strcmp(season,'annual'))
    color_alt=1;
elseif strcmp(name,'z500')
    color_alt=5;
else
    color_alt=1;    
end


if color_alt==1
    colormap(brewermap(256,'*RdBu')); % mainly used
elseif color_alt==2
    colormap(flipud(cbrewer2('div','Spectral',10))); 
elseif color_alt==3
    colormap(flipud(cbrewer2('div','PiYG',20)));
elseif color_alt==4
    colormap(flipud(cbrewer2('div','RdYlGn',20))); 
elseif color_alt==5    
    colormap(flipud(cbrewer2('div','PuOr',20)));  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hSurfAx=(gca);
cRange= caxis;

if corr_label==1
    h=colorbar('EastOutside');     
 
    set(h, 'Position', [.76 .145 .015 .75], 'FontSize',18, 'FontWeight', 'bold');  
 

 %%%%%%%%%%%% Move Colorbar  %%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%% move colobar and corr title
 
        if strcmp(name,'2mT')
            x_move_colorbar=0.060;  % 0.16; 0.11 0.085
        else   
            x_move_colorbar=0.092;  % 0.16; 0.11
        end
 
 %         if  strcmp(iso,'Accumulation')==1 || strcmp(iso,'SST_1')==1 || strcmp(iso,'SST_2')==1
            pos_c = get(h,'position'); 
            pos_c(1,1) = pos_c(1,1)+x_move_colorbar;    
            pos_c(2)=pos_c(2)+0.14;  
%             pos_c(3)=pos_c(3)- 0.0052;  % widthcolorbar

            pos_c(3)=pos_c(3)+0.014; 

            pos_c(4)=pos_c(4)-0.30;  % height colorbar
            set(h,'pos',pos_c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if era_name_nr==2
        txt_h=axestext_c(1.001,0.580,'Correlation','rotation',-90,'FontSize',18, 'FontWeight', 'bold'); % 2mT
   elseif era_name_nr==1     
        txt_h=axestext_c(1.04,0.580,'Correlation','rotation',-90,'FontSize',20, 'FontWeight', 'bold'); % z500
    else
        txt_h=axestext_c(1.03,0.58,'Correlation','rotation',-90,'FontSize',18, 'FontWeight', 'bold'); 
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%contour_alt_nr=1; % (1) p, (2) with wilks 2006, 2016 test
if contour_alt_nr==1
% contourm function keeps projection defined above
    h1= contourm( double(era_lat_c2),double(era_long),squeeze(P(1,2,:,:))', [0.05],'-k','ShowText','off',...
    'Linecolor',[.1 .1 .1],'LineWidth', 2);
elseif contour_alt_nr==2
    h1= contourm( double(era_lat_c2),double(era_long), p_fdr_rej',  [0.05],'-k','ShowText','off',...
    'Linecolor',[.1 .1 .1],'LineWidth', 2);
end

% plot3m(double(era_lat_c2(74)),double(era_long(151)),'.r','MarkerSize', 35)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
caxis(cRange);  
alpha(1)

%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(iso,'dD')==1
    iso_str='{\delta}D';
elseif strcmp(iso,'d18O')==1
    iso_str='{\delta}^1^8O';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 coast_nr=1; % (1/0) On/Off

if coast_nr==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load coast
%%%%%%%%%%%%%
% to be able to use coastline for continents in combination with ant
% coastline and grounding line from bedmap
in_c=find(lat<-60);
lat_cr=lat;
lon_cr=long;
lat_cr(in_c)=NaN;
lon_cr(in_c)=NaN;

%%%%%%%%%%%%%
color_code=[.4 .4 .4];

 plot3m(lat_cr,lon_cr,'-','LineWidth', 2,'color',color_code);
%plot3m(lat(in_c),long(in_c),'.k','MarkerSize', 1);

addpath C:\Users\benman\matlab_mapping
%addpath C:\PHD\matlab_mapping
bedmap2('gl','LineWidth', 2, 'color',color_code );
bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
show_title=3; %(1/0/2/3/4) 2- just iso string (3) use(4) season

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if show_title==1

    hp1=title([site,' ',iso_str,' ',name,' ',season,' ERA-Interim Correlation  ', num2str(round(date_annual(start_t))),' - ',num2str(round(date_annual(end_t)))],...
    'FontSize',20, 'FontWeight', 'bold');    

elseif show_title==2
  
    hp1=title([iso_str],...
    'FontSize',22, 'FontWeight', 'bold'); 

elseif show_title==3

    hp1=title(['r(',iso_str,', ',name_str,')'],...
    'FontSize',30, 'FontWeight', 'bold');   

elseif show_title==4

    hp1=title([iso_str,' ',name,' ',season],...
    'FontSize',24, 'FontWeight', 'bold');
end


    if strcmp(proj,'stereo')==1

    pos=get(hp1,'Position');
    
        if  era_name_nr==2 
            pos(2)=pos(2)-0.35;
            
        elseif era_name_nr==1     % z500
            pos(2)=pos(2)-0.35;%0.08
        else
            pos(2)=pos(2)-0.30;%0.08
        end
        
    set(hp1,'Position',pos)

    elseif strcmp(proj,'mec')==1
    pos(2)=pos(2)-0.10;
    set(hp1,'Position',pos)  
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Max corr point  %%%%%%%%%%%%%%%

if xr_lat<=-45 %&&  (strcmp(iso,'dD_1213B')==1 ||  strcmp(iso,'dD')==1)    

%  if strcmp(iso,'Accumulation')==1 &&  strcmp(name,'tp')==1 
%  max_position_str= [330 175 210 50];
%  else
 max_position_str= [160 115 190 50];
%  end
 
 
 show_max=1;
 
 if show_max==1 && strcmp(name,'z500')
         plotm(double(xr_lat),double(xr_lon),'.','MarkerSize',30,'color',[.4,.6,.1]);  % max corr point

  
%          TextBox = uicontrol('style','text');
%          set(TextBox,'String',[' r = ',num2str(round(X_r_max*1000)/1000)],...
%              'position',max_position_str,'FontWeight','bold','FontSize',22 ); % x position y position xsize ysize
%           set(TextBox,'foregroundcolor', [0 0 0], ...
%          'backgroundcolor', [1 1 1]);
 end
     
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   site_coor=star_coord_WA(site); % Site Marker


%         plotm(site_coor(1),site_coor(2),'o','MarkerSize',12,'MarkerFaceColor','none','MarkerEdgeColor',[.4,.6,.1],'LineWidth',3.3); 
%         % plotm(site_coor(1),site_coor(2),'o','MarkerSize',16,'MarkerFaceColor','none','MarkerEdgeColor',[.9,.3,.9],'LineWidth',3.2);    

        %     plotm(site_coor(1),site_coor(2),'o','MarkerSize',22,'MarkerFaceColor','none','MarkerEdgeColor',[.4,.6,.1],'LineWidth',2.5);
        
         if  era_name_nr==2 
            font_size_star=30;
         elseif era_name_nr==1
            font_size_star=36; 
         end
        
        plotm(site_coor(1),site_coor(2),'.','MarkerSize',font_size_star,'MarkerFaceColor','none','MarkerEdgeColor',[.9,.1,.9]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Regression contours %%%%%%%%%%
% Surfaces saved from ERA-I regression file
if rcontour_pc1==1
   p_level_rcontour_psa=0.05; % 0.05 or 0.01
   
    if era_name_nr==2
        filename_con=['pc_2mT_ind_psa_sum_',num2str(p_level_rcontour_psa),'_c5.mat'];
        nr_lat=121;
    elseif era_name_nr==1
        filename_con=['pc_z500_ind_psa_sum_',num2str(p_level_rcontour_psa),'_','AMJJASON','_c5.mat']; 
%         filename_con=['pc_z500_ind_psa_sum_',num2str(p_level_rcontour_psa),'_','MAMJJASON','_c5.mat']; 
        %filename_con=['pc_z500_ind_psa_sum_',num2str(p_level_rcontour_psa),'_c5.mat']; 
        % _c3 test with SAM and PSA2, doesn't explain corr pattern as well
        nr_lat=91;
    
    end
   
   folder_c='C:\Users\benman\matlab_storage_of_output_files\';
   load([folder_c,filename_con]);
  
% 1.prints  all contours 
    color_code=[0 0.8 .2]; 

    h12= contourm( double(era_lat(nr_lat:end)),double(era_long),ind_psa_sum_c',[2,2],'-','ShowText','off',...
    'Linecolor',color_code,'LineWidth', 2.5); 

% 2. green prints over 0 and -2
  % color_code= [1 .5 .1];    % overlap out-of-phase   
  % color_code= [1 0.8 .05];
%    color_code= [0.1 0.1 .1];
%     color_code= [0.3 0.3 .3];
   color_code= [1 0 0.8];
   h22= contourm( double(era_lat(nr_lat:end)),double(era_long),ind_psa_sum_c',[0,0],'-','ShowText','off',...
   'Linecolor',color_code,'LineWidth', 2.5);

% 3. prints over -2
%%%%%%%%
    color_code=[0 0.9 .2]; %
    h12= contourm(double(era_lat(nr_lat:end)),double(era_long),ind_psa_sum_c',[-2,-2],'-','ShowText','off',...
    'Linecolor',color_code,'LineWidth', 2.5);

%%%%%%%%%%%
ind_psa_sum_nonactive_c=ind_psa_sum_nonactive;
%ind_psa_sum_nonactive_c=ind_psa_sum_c;
% ind_psa_sum_nonactive_c(ind_psa_sum_nonactive_c==9999)=NaN;
ind_psa_sum_nonactive_c(ind_psa_sum_nonactive_c==0)=NaN;
% ind_psa_sum_nonactive_c(ind_psa_sum_nonactive_c==-2)=1;
ind_psa_sum_nonactive_c(ind_psa_sum_nonactive_c==2)=1;

 stipple_nr=0; % (1/0) on/off
 
 if  (stipple_nr==1 && era_name_nr==5) || (stipple_nr==1 && era_name_nr==1)
 
 [c1,c2,c3]=size(ind_psa_sum_nonactive_c);
 
 wr_c=double(era_lat(nr_lat:end));
 le_c2=length(wr_c);
 for i=1:2:c1 %360 lon
     for  k=1:2:le_c2% c2; % lat
         
         
         ind_cp=isnan(ind_psa_sum_nonactive_c(i,k));
         
%          land_ind1 = landmask(double(HadISST_lat(k)), double(HadISST_lon(i)),'Antarctica'); % check if lat long falls on land
%          land_ind2 = landmask(double(HadISST_lat(k)), double(HadISST_lon(i)),'North and South America'); % check if lat long falls on land
         
         if ind_cp==1 %&& land_ind1==0 &&  land_ind2==0
             
      %   s1=plotm(double(era_lat(k)), double(era_long(i)),'k','Marker','.','MarkerSize',4);  % marker size change comes through clearer in saved fig
         s1=plotm(wr_c(k), double(era_long(i)),'k','Marker','.','MarkerSize',4);
         end
     
     
     end
 end
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% font size for colorbar 
% placed here because it gets over writen if place earlier 

if corr_label==1
    set(h, 'FontSize',18, 'FontWeight', 'bold');
end


if strcmp(name,'z500')==1
    plabel on;
    mlabel('MLabelParallel',-18,'PLabelLocation',[-75 -45 -15 0],'PLabelMeridian',100) 

elseif strcmp(name,'2mT')==1 
    plabel on;
    mlabel('MLabelParallel',-37 ,'PLabelLocation',[-75 -60 -45 -30],'PLabelMeridian',100)
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off

if strcmp(name,'2mT')==1
    letter_pos=[220 760 60 60];  
else
    letter_pos=[200 760 60 60];  
end

         TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',40 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);


label_size=25;
            if strcmp(iso,'SST_1')==1 
            label_string='CTP SST/CP El Nino';
     t1=text(-0.82, 2.1, label_string,'FontWeight','bold','FontSize',label_size);                      
            elseif strcmp(iso,'SST_2')==1 
            label_string='CSTP SST/CP La Nina';

     t1=text(-0.82, 2.1, label_string,'FontWeight','bold','FontSize',label_size);           
            end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save figure
if corr_label==0
    filename=[site,'_',iso,'_',name,'_',season,'_ERA_interim_correlation _',...
        num2str(round(date_annual(start_t))),'_',num2str(round(date_annual(end_t))),proj];
elseif corr_label==1 % w. colorbar
    filename=[site,'_',iso,'_',name,'_',season,'_ERA_interim_correlation _',...
        num2str(round(date_annual(start_t))),'_',num2str(round(date_annual(end_t))),proj,'_crop_c2'];
end

filedir ='C:\Users\benman\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);
orient landscape

export_fig('-pdf','-painters', '-depsc','-opengl', '-r190',savefilename_c); % EPS works  func needs ghostscript  
  
%%%%%%%%%%%%%%%%%%%%

