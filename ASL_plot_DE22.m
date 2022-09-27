%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASL_plot_DE22.m
% and index
% Daniel Emanuelsson, 2021
% MATLAB 2021a
%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
% * corrcoef_df.m   UoW Steig
% * fig.m fileexchange
%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Calculate mean z500 AS Box contour area
% the start of the code from era_rice_corr_regress_annual_seasonal  to line 355
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
tic

dataset=1; % (1) ERA-interim (2) NCEP-NCAR  %%%%%%%%%%%%%
sea_nr=6;
max_ASL_point=1; % (0/1)

if dataset==1
    dataset_str='ERA-Interim';
    yr_s=1980;
  %  addpath C:\PHD\ERA_interim
    addpath C:\Users\benman\matlab_lib\Data\ERA_interim
    era_time=ncread('ERA_int_monthly_SST.nc','time');
    era_long=ncread('ERA_int_monthly_SST.nc','longitude');
    era_lat=ncread('ERA_int_monthly_SST.nc','latitude');
    era_sst=ncread('ERA_int_monthly_SST.nc','sst');
% era_t850=ncread('ERA_int_monthly_t850.nc','t');

% 'ERA_int_monthly_meridional_wind.nc');
% era_date=((double(era_time)./24)./365.2400)+1900;

    dayssince111=era_time/24;
    datevalue=dayssince111+datenum(1900,1,1);
    date_vec=datevec(double(datevalue)); 
    yyyy=date_vec(:,1); 
    mm=date_vec(:,2);
     era_year_num=yyyy+(mm/12-1/12); 

    era_z500=ncread('ERA_int_monthly_z500_2.nc','z'); % GPH z500 _2 goes to 2015
    era_z500=era_z500/9.80665;
  

elseif dataset==2
    dataset_str='NCEP-NCAR';
    yr_s=1949;
    %%%%%%%%%%%%
    % NCEP-NCAR

         nameoffile='ncep_ncar_hgt.mon.mean.nc';
        [ncep_ncar_lon, ncep_ncar_lat, ncep_ncar_time, ncep_ncar_year_num, ncep_ncar_hgt, mm]=ncep_ncar_load_lat_lon(nameoffile);
    
        % use excisting names
        era_long=ncep_ncar_lon;
        era_lat=ncep_ncar_lat;
        era_time=ncep_ncar_time; 
        era_year_num=ncep_ncar_year_num; 
        
        ncep_ncar_z500=ncread( nameoffile,'hgt'); 
        ncep_ncar_z500=squeeze(ncep_ncar_z500(:,:,6,:)); % lon,lat,level,time 
        era_z500=ncep_ncar_z500;
    
end


    era_count=find(era_year_num==2014)+11;  % data in months
%era_start=find(era_year_num==1980-1);
    era_start=find(era_year_num== yr_s-1); % annual 1980
    
div=ceil((era_count-era_start)/12);

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
 M=era_z500;
 name='z500';
 

[m n t]=size(M(:,:,1:era_count)); %
for i= 1:m
    for j= 1:n

        dummy1=reshape(squeeze(M(i,j,era_start:era_count)), 12, div);

% % %   ONDJFM      summer/fall
%         season_jfma=dummy1(1:3,:);
%         season_nd=dummy1(10:12,:);
%         season_ndjfma= [season_nd(1:div-1) ;season_jfma(:,2:div) ]; % 1980 first djf season, w. only D from 1979
%         ERA_M_season_ndjfma=nanmean(season_ndjfma(1:3,:),1); % 1980-2011
%         
%         ERA_M_season_reg_ndjfma(i,j,:)= ERA_M_season_ndjfma; % for regression
%         ERA_M_season_dummy_ndjfma=detrend(ERA_M_season_ndjfma); 
%         ERA_M_season_detrend_ondjfm(i,j,:)=ERA_M_season_dummy_ndjfma; % for corr      


% %     
% 
%    % AMJJASON Extended-Winter
          ERA_M_season_amjjason=nanmean(dummy1(4:11,:),1);  % 2. seasonal
          ERA_M_season_amjjason=ERA_M_season_amjjason(1:end); % 
          ERA_M_season_reg_amjjason(i,j,:)= ERA_M_season_amjjason; % for regression
          ERA_M_season_dummy_amjjason=detrend(ERA_M_season_amjjason); 
          ERA_M_season_detrend_amjjason(i,j,:)=ERA_M_season_dummy_amjjason;
% 
% 
%    % MAMJJA Winter Raphael etal. 2014 seasons
% %           ERA_M_season_amjjas=nanmean(dummy1(3:8,:),1);  % 2. seasonal (MAM)
% %           ERA_M_season_amjjas=ERA_M_season_amjjas(2:end); % use data starting from 1980
% %           ERA_M_season_reg_amjjas(i,j,:)= ERA_M_season_amjjas; % for regression
% %           ERA_M_season_dummy_amjjas=detrend(ERA_M_season_amjjas); 
% %           ERA_M_season_detrend_amjjas(i,j,:)=ERA_M_season_dummy_amjjas;
%           
%     % SONDJF Winter Raphael etal. 2014 seasons       
%           
% %         season_jf=dummy1(1:2,:);
% %         season_sond=dummy1(9:12,:);
% %         season_sondjf= [season_sond(1:div-1) ;season_jf(:,2:div) ]; % 1980 first djf season, w. only D from 1979
% %         ERA_M_season_sondjf=nanmean(season_sondjf(1:3,:),1); % DJF 1980-2011
% %         
% %         ERA_M_season_reg_sondjf(i,j,:)= ERA_M_season_sondjf; % for regression
% %         ERA_M_season_dummy_sondjf=detrend(ERA_M_season_sondjf); 
% %         ERA_M_season_detrend_sondjf(i,j,:)=ERA_M_season_dummy_sondjf; % for corr        
%           
%     
% %  % 1. DJF    
%         season_jf=dummy1(1:2,:);
%         season_d=dummy1(12:12,:);
%         season_djf= [season_d(1:div-1) ;season_jf(:,2:div) ]; % 1980 first djf season, w. only D from 1979
%         ERA_M_season_djf=nanmean(season_djf(1:3,:),1); % DJF 1980-2011
%         
%         ERA_M_season_reg_djf(i,j,:)= ERA_M_season_djf; % for regression
%         ERA_M_season_dummy_djf=detrend(ERA_M_season_djf); 
%         ERA_M_season_detrend_djf(i,j,:)=ERA_M_season_dummy_djf; % for corr
% %              
% %         ERA_M_annual_JF=nanmean(dummy1(1:2,:),1);  % seasonal (JF)
% %         ERA_M_annual_D=nanmean(dummy1(12,:),1);  % seasonal (D)
% % 
% %               
% %    % other seasons that aren't cut off at the beginning of the year           
% % 
% % %    % MAM
%           ERA_M_season_mam=nanmean(dummy1(3:5,:),1);  % 2. seasonal (MAM)
%           ERA_M_season_mam=ERA_M_season_mam(2:end); % use data starting from 1980
%           ERA_M_season_reg_mam(i,j,:)= ERA_M_season_mam; % for regression
%           ERA_M_season_dummy_mam=detrend(ERA_M_season_mam); 
%           ERA_M_season_detrend_mam(i,j,:)=ERA_M_season_dummy_mam;
% % % %           
% %     % JJA      
%           ERA_M_season_jja=nanmean(dummy1(6:8,:),1);  %3. seasonal (JJA)
%           ERA_M_season_jja=ERA_M_season_jja(2:end);
%              ERA_M_season_reg_jja(i,j,:)= ERA_M_season_jja; % for regression
%           ERA_M_season_dummy_jja=detrend(ERA_M_season_jja); 
%           ERA_M_season_detrend_jja(i,j,:)=ERA_M_season_dummy_jja;          
% % %           
% % %     %SON
%           ERA_M_season_son=nanmean(dummy1(9:11,:),1);  % 4.seasonal (SON)
%           ERA_M_season_son=ERA_M_season_son(2:end);
%           ERA_M_season_reg_son(i,j,:)= ERA_M_season_son; % for regression
%           ERA_M_season_dummy_son=detrend(ERA_M_season_son); 
% %           ERA_M_season_dummy_son=ERA_M_season_son; 
%           ERA_M_season_detrend_son(i,j,:)=ERA_M_season_dummy_son;    
% % % %           

% %    % AMJ fall-winter
%           ERA_M_season_fma=nanmean(dummy1(2:4,:),1);  % 2. seasonal (amj)
%           ERA_M_season_fma=ERA_M_season_fma(2:end); % use data starting from 1980
%           ERA_M_season_reg_fma(i,j,:)= ERA_M_season_fma; % for regression
%           ERA_M_season_dummy_fma=detrend(ERA_M_season_fma); 
%           ERA_M_season_detrend_fma(i,j,:)=ERA_M_season_dummy_fma;

            
%           
%  
%         %individual months, mean not necassary as we just have one month
%           ERA_M_season_march=nanmean(dummy1(3:3,:),1);  % 
%           ERA_M_season_march=ERA_M_season_march(2:end);
%              ERA_M_season_reg_march(i,j,:)= ERA_M_season_march; % for regression
%           ERA_M_season_dummy_march=detrend(ERA_M_season_march); 
%           ERA_M_season_detrend_march(i,j,:)=ERA_M_season_dummy_march;          
         
           
     % Annual   will have one extra year          
          ERA_M_annual=mean(dummy1,1); % annual ERA-values
          ERA_M_annual_reg(i,j,:)=ERA_M_annual(1:end); % save annual values for regression
          ERA_M_annual_dummy=detrend(ERA_M_annual(1:end)); % change to 1 to get 1979 annual
          ERA_M_annual_detrend(i,j,:)=ERA_M_annual_dummy;
          
          
          % Commment/uncomment code to capture variability in z500 in the AS 
          ERA_M_annual_std=nanstd(dummy1,1); 
          ERA_M_annual_reg_std(i,j,:)=ERA_M_annual_std(1:end); % save annual values for regression
        
 
    %     clear ERA_M_season_dummy_djf ERA_M_season_dummy_mam ERA_M_season_dummy_jja ERA_M_season_dummy_son    
    end
end
 
 toc
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% 
contour_use=0; % (0/1) % changed to grid
ASL_box=1; % (0/1)  Hoskins et al. 2013 (2) AS Amundsen Sea 
%%%%%%%%%%%%%%%%%%%%%%%%%%


  site= 'RI';
  
if sea_nr==1
    season='annual';
    %season='DJF';
    %season='MAM';
    %season='JJA';
    %season='SON';
elseif sea_nr==6 
    season='AMJJASON';
end

seas=1; % 1- annual
area='area_1';

% if strcmp(season,'annual')
% coord=cord_WA_poly(site);
% else
%  coord=cord_WA_poly_seasonal(site,season,area);   
% 
% % end
% 
% lon1=coord(1,1); % lon
% lon2=coord(1,2);
% lon3=coord(1,3);
% lon4=coord(1,4);
% lon5=coord(1,5);
% lon6=coord(1,6);
% 
% lat1=coord(2,1); %lat replace with
% lat2=coord(2,2);
% lat3=coord(2,3);
% lat4=coord(2,4);
% lat5=coord(2,5);
% lat6=coord(2,6);
% 
% lon_b=[lon1 lon2 lon3 lon4  lon5 lon6 lon1]; % polygone
% lat_b=[ lat1  lat2  lat3  lat4  lat5  lat6  lat1];
% 
% 


era_long=double(era_long);
era_lat=double(era_lat);

[long, latg] = meshgrid(era_long, era_lat);


if contour_use==1 && ASL_box==0
% 
% 
%  lon_b_east=lon_b+360;

contour_lim=98;
limit_c=2;
anom_tresh=100;
PPA_days=3;
% fac_3=0;
% prec_dur=2;
apro='apro_1';
 
     load(strcat('blocking\treshhold_count\contours\ind_grid_contour_',num2str(contour_lim),'_',...
site,'_',num2str(limit_c),'_',num2str( anom_tresh),...
     '_',num2str(PPA_days),'_',season,'.mat')); 
%  lat_b=contour_lat_lon(2,:);
%  lon_b_east= contour_lat_lon(1,:);

in=index_grid;

 
elseif ASL_box==1 && contour_use==0
% ASL sector z500


%%ASL%%%%%%%%%%

lon_b_east=[170 170 290 290 170];
lat_b=[ -75 -60 -60 -75 -75];
box_str='ASL';


    
 in = inpolygon(latg, long, lat_b, lon_b_east);
 in=in';    
   
 elseif ASL_box==2 && contour_use==0
    lon_b_east=[210 210 270 270 210];
    lat_b=[ -75 -60 -60 -75 -75];
    box_str='AS';
        in = inpolygon(latg, long, lat_b, lon_b_east);
        in=in';    

end
 



%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_ASL_point_str='';

if strcmp(season,'annual')==1
    
ERA_M_detrended=ERA_M_annual_detrend; % detrended time series
ERA_M_annual_reg_c=ERA_M_annual_reg;

    
elseif strcmp(season,'DJF')==1
    
ERA_M_detrended=ERA_M_season_detrend_djf; % detrended time series
ERA_M_annual_reg_c=ERA_M_season_reg_djf;  
        
elseif strcmp(season,'MAM')==1
    
ERA_M_detrended=ERA_M_season_detrend_mam; % detrended time series
ERA_M_annual_reg_c=ERA_M_season_reg_mam; 

elseif strcmp(season,'JJA')==1
    
ERA_M_detrended=ERA_M_season_detrend_jja; % detrended time series
ERA_M_annual_reg_c=ERA_M_season_reg_jja; 

elseif strcmp(season,'SON')==1
    
ERA_M_detrended=ERA_M_season_detrend_son; % detrended time series
ERA_M_annual_reg_c=ERA_M_season_reg_son; 

elseif strcmp(season,'AMJJASON')==1 && max_ASL_point==0

ERA_M_detrended=ERA_M_season_detrend_amjjason; % detrended time series
ERA_M_annual_reg_c=ERA_M_season_reg_amjjason;   

elseif strcmp(season,'AMJJASON')==1 && max_ASL_point==1 % just keeping one point the max
 
 ERA_M_season_detrend_amjjason([(1:226),(228:360)],:,:)=NaN;   
 ERA_M_season_detrend_amjjason(:,[(1:151),(153:180)],:)=NaN; 

 ERA_M_season_reg_amjjason([(1:226),(228:360)],:,:)=NaN;   
 ERA_M_season_reg_amjjason(:,[(1:151),(153:180)],:)=NaN;  
 
ERA_M_detrended=ERA_M_season_detrend_amjjason; 
ERA_M_annual_reg_c=ERA_M_season_reg_amjjason;     
    
% ERA_M_detrended=ERA_M_season_detrend_amjjason(227,151,:); 
% ERA_M_annual_reg_c=ERA_M_season_reg_amjjason(227,151,:); 

max_ASL_point_str='_max';
end    


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cut out SH part to get on the same size grid

if ASL_box==0 && contour_use==1

        ERA_M_detrended_c=ERA_M_detrended(:,161:238,:);
        ERA_M_annual_reg_c2=ERA_M_annual_reg_c(:,161:238,:);

elseif  (ASL_box==1 || ASL_box==2 ) && contour_use==0 
    
        ERA_M_detrended_c=ERA_M_detrended(:,:,:);
        ERA_M_annual_reg_c2=ERA_M_annual_reg_c(:,:,:);
    
end

% ERA_M_annual_reg_c

%%%
% lon_b=[lon_box1 lon_box1 lon_box1_2 lon_box1_2 lon_box1];
% lat_b=[ lat_box1_2 lat_box1 lat_box1 lat_box1_2 lat_box1_2];
% 
% 
% lon_box1c=180+lon_box1+180;
% lon_box1_2c=180+lon_box1_2+180;
% 
% ind_lon = find(era_long<=lon_box1_2c & era_long>=lon_box1c);
% ind_lat =find(era_lat<=lat_box1 & era_lat>=lat_box1_2);

le=length(ERA_M_annual);

era_x_polygon_mean_detrended=zeros(le,1);
era_x_polygon_mean_reg=zeros(le,1);

[a,b,c]=size(ERA_M_detrended_c);

for i=1:c
    
   era_x_event_mean_detrended=ERA_M_detrended_c(:,:,i); 
   era_x_event_mean_detrended_c=era_x_event_mean_detrended;
   era_x_event_mean_detrended_c(~in) = NaN; 
   
     era_x_event_mean_reg=ERA_M_annual_reg_c2(:,:,i); 
     era_x_event_mean_reg_c=era_x_event_mean_reg;
     era_x_event_mean_reg_c(~in) = NaN; 
   
     era_x_polygon_mean_detrended(i)=nanmean(nanmean(era_x_event_mean_detrended_c,2),1);
     era_x_polygon_mean_reg(i)=nanmean(nanmean(era_x_event_mean_reg_c,2),1);
   
% 
% era_x_event_mean_reg=nanmean(nanmean(ERA_M_annual_reg(ind_lon,ind_lat,i),2),1);
% 
% era_z500_AS_time_series_detrended(i,1)=era_x_event_mean_detrended;
% era_z500_AS_time_series_reg(i,1)=era_x_event_mean_reg;

end



if strcmp(season,'annual')==1

era_MA_z500_annual_time_series=zeros(le,3);
    if le==36
    era_MA_z500_annual_time_series(:,1)=1979:2014;
    elseif  le==67
    era_MA_z500_annual_time_series(:,1)=1948:2014;
    end
era_MA_z500_annual_time_series(:,2)=era_x_polygon_mean_detrended;
era_MA_z500_annual_time_series(:,3)=era_x_polygon_mean_reg;


elseif strcmp(season,'AMJJASON')==1
era_MA_z500_annual_time_series=zeros(36,3);
era_MA_z500_annual_time_series(:,1)=1979:2014;
era_MA_z500_annual_time_series(:,2)=era_x_polygon_mean_detrended(1:36);
era_MA_z500_annual_time_series(:,3)=era_x_polygon_mean_reg(1:36);

else  %1980-2013

era_MA_z500_annual_time_series=zeros(34,3);
era_MA_z500_annual_time_series(:,1)=1980:2013;
era_MA_z500_annual_time_series(:,2)=era_x_polygon_mean_detrended(1:34);
era_MA_z500_annual_time_series(:,3)=era_x_polygon_mean_reg(1:34);    
    
end  


save_nr=1;

if contour_use==1 && ASL_box==0 && save_nr==1
    savefilename = strcat('blocking\era_MA_z500_',site,'_contour_',num2str(contour_lim),...
    '_',num2str(limit_c),'_',num2str( anom_tresh),...
     '_days1_',num2str(PPA_days),season,'_',apro,'.mat'); 
    % %matlab workspace
    save(savefilename,'era_MA_z500_annual_time_series');
 
elseif (ASL_box==1 || ASL_box==2) && contour_use==0 && save_nr==1
     
     
     %savefilename = strcat('blocking\',name,'ASL_sector_P_Hoskings.mat'); 
     
    savefilename = strcat('C:\Users\benman\matlab_storage_of_output_files\',name,'_', dataset_str,'_',season,'_', box_str,'_sector',max_ASL_point_str,'.mat'); 
     
    % %matlab workspace
    save(savefilename,'era_MA_z500_annual_time_series');
     
end
     

% Trends linear least-square
% date_annual=[1979:2013];

 if strcmp(season,'annual')==1
     date_annual=[1979:2014];
start_t=find(date_annual==1979);
end_t=find(date_annual==2014);

 else
     date_annual=[1980:2013];
start_t=find(date_annual==1980);
end_t=find(date_annual==2010);
 end

  y=era_MA_z500_annual_time_series((start_t:end_t),3);

date_annual(start_t);
date_annual(end_t);

x=date_annual(start_t:end_t);
%y=BPA_counts_box(start_t:end_t);
%y=anomaly(D_del_stacked_annual(start_t:end_t));% standardized run for
%table

% [s b]=trend(y,[],1);
% s10=s*10 % trend per dec


s=regstats(y,x,'linear','all');
s.tstat;
s10=s.tstat.beta(2)*10; % per dec
p_iso_rec=s.tstat.pval;
p_iso_rec(2)
s10;
% The second element is the p-value for a significantly non-zero slope.
pc=p_iso_rec(2) 
gh=[s10 pc];
gh

figure
plot(x,y,'*-k')

 clear filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig ASL
%load('z500_ASL_c')
% load('C:\PHD\matlab_storage_of_output_files\z500_ERA-Interim_ASL_sector_P_Hoskings.mat')

sea_nr=6;
max_nr=1;
if sea_nr==1 && max_nr==0
    load('C:\Users\benman\matlab_storage_of_output_files\z500_ERA-Interim_ASL_sector_P_Hoskings.mat')
    season='annual';
    max_str='';
elseif sea_nr==6 && max_nr==0
    load('C:\Users\benman\matlab_storage_of_output_files\z500_ERA-Interim_AMJJASON_ASL_sector.mat')
    season='AMJJASON';
    max_str='';
elseif sea_nr==6 && max_nr==1    
    load('C:\Users\benman\matlab_storage_of_output_files\z500_ERA-Interim_AMJJASON_ASL_sector_max.mat')    
    season='AMJJASON';
    max_str='_max';
end

h=fig('units','inches','width',18,'height',5,'font','Helvetica','fontsize',18,'border','on');
   box on
   hold on   
%%%%%%%%%%%%%% 
color_code=[0.0,0.0,0.4];
years=[1979:2011]';
yyaxis left

% h7=plot(years, z500_ASL((1:33),1),'--','LineWidth',4,'MarkerSize',14,'Color',color_code); 
% hold on

h7=plot(era_MA_z500_annual_time_series((1:33),1), era_MA_z500_annual_time_series((1:33),3),'-','LineWidth',2.5,'MarkerSize',14,'Color',color_code); 


     % second y-axis color
       ax = gca;
       ax.YColor = color_code;    



         ylabel('Geopotential height (m)','FontWeight','bold','FontSize',20 );
                 
%%%%%%%%%%
        iso_alt_nr=2;
        
        if iso_alt_nr==1
             load('C:\Users\benman\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c23.mat');
        elseif iso_alt_nr==2
%             load('C:\Users\benman\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_c24.mat'); % MAMJJASON
            load('C:\Users\benman\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_AMJJASON_c24.mat'); % AMJJASON 
%             load('C:\Users\benman\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_clim_c24.mat'); % clim removed, no difference
%             
        end
% 
% 
% 
 yyaxis right
% 
color_code2=[0.5,0.5,0.5];
plot(MA_save((1980:end),1),MA_save((1980:end),3),':','LineWidth',4,'Color',color_code2)

     ylabel(['{\delta}^1^8O (',char(8240),')'],'FontWeight','bold','FontSize',20 );

% set(gca,'YLim',[-1 1])
 set(gca,'XLim',[1978 2012])
%  set(gca,'YLim',[5020 5260])
set(gca,...
          'linewidth',2,...
          'FontWeight','bold' );

      
%       %%%%%%%%%%%%%%%%%%%%%%%%
%      % second y-axis color
       ax = gca;
       ax.YColor = color_code2;     


 xlabel('Year','FontWeight','bold','FontSize',20 );
         
 title_str='Z500 ASL';   
 axestext_c(0.9999999,1.0,title_str,'FontWeight','bold','FontSize',20 );          
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x y width height] for the figure, increase margins
    pos1 = get(gca, 'Position');
    pos1(1) = 0.2;
    pos1(3) = 0.6;
    set(gca, 'Position', pos1)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %letter='a';
  letter='b';
 %letter_position=   [240 380 60 60];
 letter_position=   [200 380 50 60];
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_position,'FontSize',32); % x position y position xsize ysize ,'FontWeight','bold'
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
% save fig.
end_str='_c2';
filename=['ASL_plot','_',season,end_str,max_str]; 
filedir ='C:\Users\benman\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);
% save as png 
orient landscape
%export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
export_fig('-pdf','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c); % PNG-nocrop'
%

%%

% Y=detrend(ASL_SectorP_x(1:33,1));
% Y=detrend(ASL_RelCenPres_x(1:33,1));
Y=detrend(era_MA_z500_annual_time_series(1:33,2));

%
X=detrend(MA_save(1980:2012,3));

[r,p, nv]=corrcoef_df_c2(X,Y);
rp=[r(2),p(2)]
p_level(p(2))


