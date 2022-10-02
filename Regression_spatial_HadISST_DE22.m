%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regression_spatial_HadISST_DE22.m
% Matlab 2021a
% Daniel Emanuelsson 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP environment ----
close all
fprintf('[\bSetting up environment ... ]\b')
tic
% CLEAR environment
clearvars -except config ; close all

% ESTABLISH configuation
% If running from master script ELSE user input the config file
if exist('config','var')
    eval(config)
else 
    addpath('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab\configfiles\')

  eval('Config_reg_SST')

 
end
toc




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input

type_nr=2;    % (1) annual (2) monthly

param_nr=2;% (1-SST/ 2 SIC)

% if param_nr==1
%     name='SST';
% elseif param_nr==2
    name='SIC';
% end
 
 site='RICE';
 
 set_yr=1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if set_yr==1    
        yr_s=1979;
        yr_e=2014;     
   elseif set_yr==8       
        yr_s=2000;
        yr_e=2013;
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
 
 show_colorbar=1;
  
 figure_format=2; %(1)  EPS (2) PNG
  
 sea_nr=1;

 season='annual';

 
 %%%%%%%%%%%%%%%%%%%%
 iso_nr=15;  
 %%%%%%%%%%%%%%%%%%%
 
 
if iso_nr==7 %***************
iso='PCs';
elseif iso_nr==9 %***************
iso=[]; % just trend
elseif iso_nr==14 %***************
iso='CP EP index';
elseif iso_nr==15 %***************
iso='SIC PCs';
end
 
show_maximum_point=0;
 
show_title=1; % (1/0)
area_2_box=0 ; % 0/1 for RICE SST corr Box and text for Area 2
  
if  strcmp(name,'SIC')==1
  
  lat1=-90;
  lat2=-50;
  box_use=0;  % (0/1) 

  lon1=-300;
   lon2= 60;
  
  
% elseif  strcmp(name,'SST')==1  
%   proj='merc';
%   %proj='stereo';
%   lat1=-90;
%   lat2=20;
%   
%   lon1=-300;
%   %lon2= -30;
%    lon2= 60;
%   
%    box_use=0;  % (0/1)
end
  
%   letter=[];
 
  lock_scalebar=1;
  
  
%   if  strcmp(season,'DJF')==1
%   letter='a';
%   elseif strcmp(season,'MAM')==1
%   letter='b';
%   elseif strcmp(season,'JJA')==1
%   letter='c';
%   elseif strcmp(season,'SON')==1
%   letter='d';
%   end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SST_dataset=1; % (1) HadISST 

if SST_dataset==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HadISST data
% version 1.1 Rayner et al. 2003


addpath C:\PHD\HadISST_data\ % Ienovo
%    ncdisp('HadISST_ice_c.nc');
%    ncdisp('HadISST_sst.nc');

 HadISST_time=ncread('HadISST_ice_c.nc','time'); %units         = 'days since 1870-1-1 0:0:0'
 % missing values  missing_value = -1e+30 ( looks more like -1000)\
 %monthly values
 
 HadISST_time_c=HadISST_time+datenum(1870,1,1);  % to 2015 month 7
 
 date_vec=datevec(double(HadISST_time_c)); 
 yyyy=date_vec(:,1); 
 mm=date_vec(:,2);
 HadISST_year_num=yyyy+(mm/12-1/12); 

 HadISST_time_bnds =ncread('HadISST_ice_c.nc','time_bnds');
 HadISST_lat=ncread('HadISST_ice_c.nc','latitude');
 HadISST_lon=ncread('HadISST_ice_c.nc','longitude');
 
 
end 
 
    


 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  if strcmp(name,'SST')==1
%  %%%%%%%%%% file    _c contains data until 2015
%   %name='SST';
%   if SST_dataset==1
%   HadISST_sst=ncread('HadISST_sst_c.nc','sst'); % 360x180x1728 lon x lat x time
%   M=HadISST_sst;
%   
%   elseif SST_dataset==2
%   ERSST_sst=ncread('sst.mnmean.v4.nc','sst');
%   HadISST_sst=ERSST_sst;
%   M=ERSST_sst;
%   end
  
 if strcmp(name,'SIC')==1
  
%   name='ice'; 
   HadISST_ice=ncread('HadISST_ice_c.nc','sic'); % 360x180x1728 lon x lat x time
   M=HadISST_ice;
 
 end
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % 4. 0 Seasonal and Annual
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %
 % HadISST_count=1705; % that overlaps with RICE record 1979-2011
HadISST_count=find(HadISST_year_num==yr_e)+11; % values in months
% HadISST_start=289; %1894
HadISST_start=find(HadISST_year_num==yr_s); % 1950 hardly any data before this
% HadISST_start=1309; %1979 satelite
% (HadISST_count-HadISST_start)/12

% era_count=1693; 
% div=118;
div=ceil((HadISST_count-HadISST_start)/12); %1979

%%%%%%%%%%%%%%%%%%%%%%%
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%


mm_in=mm(HadISST_start:HadISST_count);  % seasonal index

 if sea_nr==1
    mm_in_c=[HadISST_start:HadISST_count]';
    
%  elseif sea_nr==2  % DJF
%      mm_in_c=[find(mm_in==12);  find(mm_in==1); find(mm_in==2)];
%      mm_in_c=sort( mm_in_c);
%      mm_in_c=mm_in_c(3:end-1); % dont use JF for first year as there is no D
%      
% elseif sea_nr==3  % MAM
%      mm_in_c=[find(mm_in==3);  find(mm_in==4); find(mm_in==5)];
%      mm_in_c=sort( mm_in_c);
%      mm_in_c=mm_in_c(4:end);
%      
% elseif sea_nr==4  % JJA
%      mm_in_c=[find(mm_in==6);  find(mm_in==7); find(mm_in==8)];
%      mm_in_c=sort( mm_in_c);
%      mm_in_c=mm_in_c(4:end);
%      
%  elseif sea_nr==5  % SON
%      mm_in_c=[find(mm_in==9);  find(mm_in==10); find(mm_in==11)];
%      mm_in_c=sort( mm_in_c);
%      mm_in_c=mm_in_c(4:end);
        
 end
%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

% % cosweight
% S_lon=17;S_lat=170;%check
% M ( S_lon, S_lat-5,2) 
%  M=cosweight_c(M,HadISST_lat); % UoW function
%  M ( S_lon, S_lat-5,2)
 
 
%%%

HadISST_year_num(HadISST_start);
HadISST_year_num(HadISST_count);

% HadISST_year_num()-HadISST_year_num()
% HadISST_start:HadISST_count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run from start if year is changed above
% (360x180x1728) (Long,Lat,month)
%%
[m n t]=size(M(:,:,HadISST_start:HadISST_count));% time(months) % 1894-2011
for i= 1:m
    for j= 1:n

        dummy1=reshape(squeeze(M(i,j,HadISST_start:HadISST_count)), 12, div);

%         xr_c=find(squeeze(M(i,j,HadISST_start:HadISST_count))<=-800);
%         fill_index(i,j,:)= xr_c;

%           Summer
%      % 1. ONDJFM     
%         season_jfm=dummy1(1:3,:);
%         season_ond=dummy1(10:12,:);
%         season_ondjfm= [season_ond(1:div-1) ;season_jfm(:,2:div) ]; % 1980 first djf season, w. only D from 1979
%         HadISST_M_season_ondjfm=nanmean(season_ondjfm(1:3,:),1); % DJF 1980-2011
%         HadISST_M_season_reg_ondjfm(i,j,:)=HadISST_M_season_ondjfm; % means for regression
%         HadISST_M_season_dummy_ondjfm=detrend(HadISST_M_season_ondjfm); 
%         HadISST_M_season_detrend_ondjfm(i,j,:)=HadISST_M_season_dummy_ondjfm; % d


%           Summer
%      % 1. SONDJF   Raphael 2014  Retreat
%         season_jfm=dummy1(1:2,:);
%         season_ond=dummy1(9:12,:);
%         season_ondjfm= [season_ond(1:div-1) ;season_jfm(:,2:div) ]; % 1980 first djf season, w. only D from 1979
%         HadISST_M_season_ondjfm=nanmean(season_ondjfm(1:3,:),1); % DJF 1980-2011
%         HadISST_M_season_reg_ondjfm(i,j,:)=HadISST_M_season_ondjfm; % means for regression
%         HadISST_M_season_dummy_ondjfm=detrend(HadISST_M_season_ondjfm); 
%         HadISST_M_season_detrend_ondjfm(i,j,:)=HadISST_M_season_dummy_ondjfm; % d
% 
% 
% 
% 
% 
% %  % Winter season Raphael 2014
%  %  Advance MAMJJA

if strcmp(season,'MAMJJA')==1

          HadISST_M_season_amjjas=nanmean(dummy1(3:8,:),1);  %3. seasonal (JJA)
          HadISST_M_season_amjjas=HadISST_M_season_amjjas(1:end);
          HadISST_M_season_reg_amjjas(i,j,:)=HadISST_M_season_amjjas; % means for regression         
          HadISST_M_season_dummy_amjjas=detrend(HadISST_M_season_amjjas); 
          HadISST_M_season_detrend_amjjas(i,j,:)=HadISST_M_season_dummy_amjjas;
%         
%         
%  % Winter
       % JJA      
%           HadISST_M_season_amjjas=nanmean(dummy1(6:11,:),1);  %3. seasonal (JJA)
%           HadISST_M_season_amjjas=HadISST_M_season_amjjas(2:end);
%           HadISST_M_season_reg_amjjas(i,j,:)=HadISST_M_season_amjjas; % means for regression         
%           HadISST_M_season_dummy_amjjas=detrend(HadISST_M_season_amjjas); 
%           HadISST_M_season_detrend_amjjas(i,j,:)=HadISST_M_season_dummy_amjjas;        
        
    
        
% 
% %          % 1. DJF     
    
elseif strcmp(season,'DJF')==1

        season_jf=dummy1(1:2,:);
        season_d=dummy1(12:12,:);
        season_djf= [season_d(1:div-1) ;season_jf(:,2:div) ]; % 1980 first djf season, w. only D from 1979
        HadISST_M_season_djf=nanmean(season_djf(1:3,:),1); % DJF 1980-2011
        HadISST_M_season_reg_djf(i,j,:)=HadISST_M_season_djf; % means for regression
        HadISST_M_season_dummy_djf=detrend(HadISST_M_season_djf); 
        HadISST_M_season_detrend_djf(i,j,:)=HadISST_M_season_dummy_djf; % detrended data for corr
% %         
% %         
% %            % MAM

elseif strcmp(season,'MAM')==1
          HadISST_M_season_mam=nanmean(dummy1(3:5,:),1);  % 2. seasonal (MAM)
          HadISST_M_season_mam=HadISST_M_season_mam(2:end); % use data starting from 1980
          HadISST_M_season_reg_mam(i,j,:)=HadISST_M_season_mam; % means for regression          
          HadISST_M_season_dummy_mam=detrend(HadISST_M_season_mam); 
          HadISST_M_season_detrend_mam(i,j,:)=HadISST_M_season_dummy_mam;
% %           
% %     % JJA      
elseif strcmp(season,'JJA')==1
          HadISST_M_season_jja=nanmean(dummy1(6:8,:),1);  %3. seasonal (JJA)
          HadISST_M_season_jja=HadISST_M_season_jja(2:end);
         HadISST_M_season_reg_jja(i,j,:)=HadISST_M_season_jja; % means for regression         
          HadISST_M_season_dummy_jja=detrend(HadISST_M_season_jja); 
          HadISST_M_season_detrend_jja(i,j,:)=HadISST_M_season_dummy_jja;          
% %           
% %     %SON
elseif strcmp(season,'SON')==1
          HadISST_M_season_son=nanmean(dummy1(9:11,:),1);  % 4.seasonal (SON)
          HadISST_M_season_son=HadISST_M_season_son(2:end);
          HadISST_M_season_reg_son(i,j,:)=HadISST_M_season_son; % means for regression          
          HadISST_M_season_dummy_son=detrend(HadISST_M_season_son); 
          HadISST_M_season_detrend_son(i,j,:)=HadISST_M_season_dummy_son;    
%         
%     % S month with maximum SIE
%           
%           HadISST_M_season_s=nanmean(dummy1(11:11,:),1);  % 5. individual month
%           HadISST_M_season_s=HadISST_M_season_s(2:end);
%           HadISST_M_season_reg_s(i,j,:)=HadISST_M_season_s; % means for regression
%           HadISST_M_season_dummy_s=detrend(HadISST_M_season_s); 
%           HadISST_M_season_detrend_s(i,j,:)=HadISST_M_season_dummy_s; 
% %     
%     
elseif strcmp(season,'annual')==1
         HadISST_M_annual_temp=nanmean(dummy1,1); % annual HadISST-values

         
%          HadISST_M_annual_temp=anomaly(HadISST_M_annual_temp); % try
%          normalizing
         HadISST_M_annual_dummy=detrend(HadISST_M_annual_temp(1:end));
         
         HadISST_M_annual(i,j,:)=HadISST_M_annual_temp(1:end);  %reg
         HadISST_M_annual_detrend(i,j,:)= HadISST_M_annual_dummy; % detrended
%          
end

%   clear       dummy1 HadISST_M_annual_dummy HadISST_M_annual_temp

    end
end
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(iso,'dD')==1 || strcmp(iso,'d18O')==1 || strcmp(iso,'d-excess')==1
    
%     iso=' dD ';

% isotopes stacked record
% load('ISO_rice_annual_c3.mat') % only load to get date vector
% date_annual=[1980:2012]';
date_annual=[1900:2011]';

if  strcmp(season,'annual')==1 ||  strcmp(season,'MAMJJA')==1
   start_t=find(date_annual==yr_s);
elseif strcmp(season,'DJF')==1 || strcmp(season,'MAM')==1 || strcmp(season,'JJA')==1 || strcmp(season,'SON')==1   
    start_t=find(date_annual==yr_s+1);
end
% start_t=86; % # 86-1979
end_t=find(date_annual==yr_e); 

load('stacked_record\stacked_record_annual_Ma_c14.mat') % c12 stacked, c14 just 1213B and 2013 deep
% updated age-scale

    if   strcmp(iso,'dD')==1 
    X=detrend(stacked_record_annual_Ma((start_t:end_t),2));

    elseif strcmp(iso,'d18O')==1
    X=detrend(stacked_record_annual_Ma((start_t:end_t),3));  
    
    elseif strcmp(iso,'d-excess')==1
    X=detrend(stacked_record_annual_Ma((start_t:end_t),4));         
    end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(iso,'Accumulation')==1
% iso=' Accumulation ';
% Accumulation
% accumulation record (:,3) is detrended
% load('rice_accumulation\Accumulation_we_RICE_top40m_1979_c4.mat')
% iso=' Accumulation ';
% 
% start_t=find(Accumulation_we_RICE_top40m_1979(:,1)==1980);
% % start_t=86; % # 86-1979
% end_t=find(Accumulation_we_RICE_top40m_1979(:,1)==2010); 
% 
% date_annual=Accumulation_we_RICE_top40m_1979(:,1);
% X=Accumulation_we_RICE_top40m_1979(:,3);
% X=X(start_t:end_t);

load('rice_accumulation\Accumulation_we_RICE_top40m_c7.mat')

date_annual=Accumulation_we_RICE_top40m(:,1);
if  strcmp(season,'annual')==1 ||  strcmp(season,'MAMJJA')==1
   start_t=find(date_annual==yr_s);
elseif strcmp(season,'DJF')==1 || strcmp(season,'MAM')==1 || strcmp(season,'JJA')==1 || strcmp(season,'SON')==1   
    start_t=find(date_annual==yr_s+1);
end
end_t=find(date_annual==yr_e); 


X=Accumulation_we_RICE_top40m(:,3);   %  (:,3) is detrended
% X=X(find(Accumulation_we_RICE_top40m_1979(:,1)==1979):find(Accumulation_we_RICE_top40m_1979(:,1)==2011));
X=X(start_t:end_t);


elseif strcmp(iso,'d_excess_1213B') % just using 1213B core
   
load('stacked_record\shallow_core_1213B_annual_ISO_RICE_Ma_c.mat')
date_annual=stacked_record_annual_Ma(:,1);
if  strcmp(season,'annual')==1 ||  strcmp(season,'MAMJJA')==1
   start_t=find(date_annual==yr_s);
elseif strcmp(season,'DJF')==1 || strcmp(season,'MAM')==1 || strcmp(season,'JJA')==1 || strcmp(season,'SON')==1   
    start_t=find(date_annual==yr_s+1);
end
end_t=find(date_annual==yr_e); 
X=detrend(stacked_record_annual_Ma((start_t:end_t),4)); % 4 d-excess     
 

elseif strcmp(iso,'dD_1213B') % just using 1213B core
load('stacked_record\shallow_core_1213B_annual_ISO_RICE_Ma_c.mat')
date_annual=stacked_record_annual_Ma(:,1);
if  strcmp(season,'annual')==1 ||  strcmp(season,'MAMJJA')==1
   start_t=find(date_annual==yr_s);
elseif strcmp(season,'DJF')==1 || strcmp(season,'MAM')==1 || strcmp(season,'JJA')==1 || strcmp(season,'SON')==1   
    start_t=find(date_annual==yr_s+1);
end
end_t=find(date_annual==yr_e); 
X=detrend(stacked_record_annual_Ma((start_t:end_t),2)); % 4 dD 1213B 

elseif strcmp(iso,'d18O_1213B') % just using 1213B core
   
load('stacked_record\shallow_core_1213B_annual_ISO_RICE_Ma_c.mat')
date_annual=stacked_record_annual_Ma(:,1);
if  strcmp(season,'annual')==1 ||  strcmp(season,'MAMJJA')==1
   start_t=find(date_annual==yr_s);
elseif strcmp(season,'DJF')==1 || strcmp(season,'MAM')==1 || strcmp(season,'JJA')==1 || strcmp(season,'SON')==1   
    start_t=find(date_annual==yr_s+1);
end
end_t=find(date_annual==yr_e); 
X=detrend(stacked_record_annual_Ma((start_t:end_t),3)); % 3 d18O 1213B 


elseif strcmp(iso,'PPAp')==1

%      load('PPAp_RI_3days_2.mat')  % 95 100 
     load('PPA_RI_mean_3days.mat') % PPA AS

    date_annual=[1979:2014]';
    start_t=find(date_annual==yr_s); 
    end_t=find(date_annual==yr_e); 
    
     y=detrend(PPA_box_mean_store(start_t:end_t));
     
elseif iso_nr==7

%     load('C:\PHD\matlab_storage_of_output_files\PCs_z500.mat')  
%load('C:\PHD\matlab_storage_of_output_files\PCs_z500_lim0-360E_-30_-90_1979-2015.mat');
% load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2015_annual');
load('C:\PHD\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2014_annual'); % 2014
        if type_nr==1 
            time_c=MA_PCs_save(:,1);
            yr_2=find( time_c==yr_e);
           
            
            
            
        elseif type_nr==2
            time_c=MA_PCs_save_monthly(:,1);
            yr_2=find( time_c==yr_e)+11;
        end
            
     yr_1=find( time_c==yr_s); 
  
                 
     
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            PC_nr=4; % (2) SAM (3) PSA1 (4) PSA2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
             label_2='';
             iso='z500';
   
            if PC_nr==2 %SAM
                
                    letter='c';
                    
                      
                if sea_nr==1
                    if param_nr==1
                    letter='a';
                        fact=-1;  %  SAM+ 2015
                    elseif param_nr==2  %SIC
                    letter='d';
%                     fact=1;%%%%%%%%%%%%%%%%%change polarity SAM-%%%%%%%%%%%%%%%%%%%%%%%
                        % fact=-1; % to corespond to SIC PC regression fig.
                        % 2015
                         fact=1; % to corespond to SIC PC regression fig. 2014
                    end
                 
                elseif sea_nr==2
                letter='b';
                elseif sea_nr==3
                 letter='b';       
               elseif sea_nr==4
                letter='b';        
               elseif sea_nr==5
               letter='b';  
                end
                
            
                
          
                
                
            elseif PC_nr==3
                
                        
                
                if sea_nr==1
                   if param_nr==1 
                    letter='d';
                    fact=-1;   %PSA1+
                   elseif param_nr==2
                     letter='e';
                    % fact=1;   %PSA1   % to corespond to SIC PC regression fig. 2015
                     fact=1;   %PSA1   % to corespond to SIC PC regression fig. 2014
                   end
                 
                elseif sea_nr==2
                letter='c';
                elseif sea_nr==3
                 letter='c';       
               elseif sea_nr==4
                letter='c';        
               elseif sea_nr==5
               letter='c';  
                end
                 
%                             if sea_nr==2
%                             fact=1;   %PSA1+ (as it doesnt resemble PSA1 in DJF)
%                             else
%                             fact=1;   %PSA1+   
%                             end
                 
            elseif PC_nr==4
                 
                                      % Changed sign convetion            
                
                % use -1 because wave train eminates from Aus positive SSTA?
                % fact=1; %PSA2-  
                
                 if sea_nr==1
                     if param_nr==1
                        letter='e';
                        fact=1; %PSA2+ 
                     elseif param_nr==2
                        letter='f';
                        % fact=-1; % to corespond to SIC PC regression fig.2015
                         fact=-1; % to corespond to SIC PC regression fig.2014
                         
                     end
                elseif sea_nr==2
                letter='d';
                elseif sea_nr==3
                 letter='d';       
               elseif sea_nr==4
                letter='d';        
               elseif sea_nr==5
               letter='d';  
                 end
                
     
                
                % use as anticyclone in the Ross Sea likely to be associated with enriched isotopes. 
            end
            
            
            
         if type_nr==1   && sea_nr==1  % annual
            y=MA_PCs_save((yr_1:yr_2),PC_nr)*fact;
            
              
         elseif  type_nr==1  && sea_nr>=2  % seasonal one value per year
             
          
            y_c2=MA_PCs_save_monthly((mm_in_c),PC_nr)*fact;
            
            
             y_c3=reshape(  y_c2,3,length( y_c2)/3);
            
              y_c4=nanmean( y_c3)';
              
              y= y_c4; % seasonal start 1980 allready taken out
             
            %%%%%%%%%%%%%%%%%%%%%%% 
            
         elseif type_nr==2 && sea_nr==1 %  monthly
            y=MA_PCs_save_monthly((yr_1:yr_2),PC_nr)*fact;
            
           elseif type_nr==2 && sea_nr>=2 %  monthly seasonal 
               
             y=MA_PCs_save_monthly((mm_in_c),PC_nr)*fact;   
            
         end
         
elseif iso_nr==8
     
    nino_nr=4; % (34) nino-3.4 (4) nino-4 (12) nino 1+2
     
     PC_nr=999; % just to get a number for enso
     label_2='';
%      iso='z500';
     
     letter='a'; 
     
     
    if nino_nr==4    % nino-4
    load('nino_4_M_c.mat'); % 2015
    nino_xc=nino_4_M;
     iso='ENSO_nino-4';
    elseif nino_nr==34
    load('nino_34_M.mat'); % nino 3.4 
    nino_xc=nino_34_M;  
    iso='ENSO_nino-34';
    
    elseif nino_nr==12 % nino 1+2
     load('C:\PHD\matlab_lib\Data\Nino_1p2_noaa.mat');
     nino_xc=Nino_1p2_noaa(:,(2:13));
     iso='ENSO_nino-1+2';
    
    end
     

% sea_nr=1; %% (1), annual (2) DJF
        date_annual=[1870:2014]';
        start_t=find(date_annual==yr_s); 
        end_t=find(date_annual==yr_e); 

% yr_1
% yr_2        
        
         if sea_nr==1
                s1=1;s2=12;
         elseif sea_nr==3
                s1=3;s2=5;
        elseif sea_nr==4
                s1=6;s2=8;
        elseif sea_nr==5
            s1=9;s2=11;    
         end

    if sea_nr==2
    %%%%%%%%DJF
    s1=1; s2=2; s3=12;
    nino_sea_jf=nino_xc((2:end),(s1:s2));
    nino_sea_d=nino_xc((1:end-1),(s3)); 
    nino_sea_djf=horzcat(nino_sea_d, nino_sea_jf);
    nino_x=nanmean(nino_sea_djf,2);
    %%%%%%%%%
    else
    nino_sea=(nino_xc((1:end),(s1:s2))); % ok
    nino_x=nanmean(nino_sea,2);
    end
    %%%%%%%%%%%%%%%%%% restrict years first
        nino_sea_c=nino_sea ((start_t:end_t),:);
    
        nino_sea_c2=permute(nino_sea_c,[2 1]);
    
        [w_a, w_b]=size(nino_sea_c2);
        nino_sea_c3=reshape(nino_sea_c2, w_b*w_a,1);
         nino_sea_c4=zscore( nino_sea_c3);
    

%     time_series=nino_x(start_t:end_t); % 1979-2009
        y=nino_sea_c4;  
    
  elseif iso_nr==9
      
      
            wr_yr=[1979:2014]';
            yr_1=find(wr_yr==yr_s);       
            yr_2=find(wr_yr==yr_e);
            
                if sea_nr==1 && type_nr==1% annual
 
                y=wr_yr(yr_1:yr_2);
                elseif sea_nr>=2  && type_nr==1        % Seasonal (1980-
                y=wr_yr(yr_1+1:yr_2);
                elseif  type_nr==2% annual and seasonal monthly
                yc=[1:HadISST_count]';
                y=yc(mm_in_c);
                
                
%                elseif sea_nr>=2 && type_nr==2% annual monthly
%                 y=[1:HadISST_count]';
                
        
                end
                
                
                 PC_nr=999; % just to get a number for enso    
                
                
            label_2='yr^-^1';
            
                 if sea_nr==1 && param_nr==1
                 letter='a';
                 elseif sea_nr==1 && param_nr==2
                     if set_yr==8
                         letter='a';
                     elseif set_yr==9
                         letter='d'; 
                     elseif set_yr==1
                         letter='a'; 
                     end
                 
%                  elseif sea_nr==2
%                  letter='a';
%                 elseif sea_nr==3
%                  letter='b'; 
%                  elseif sea_nr==4
%                  letter='c'; 
%                  elseif sea_nr==5
%                  letter='d'; 
                 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
   elseif iso_nr==10      %PDO
          
            if sea_nr==1 && type_nr==1
                
            
            date_annual=[1900:2015]';

            load('PDO_annual_c.mat')
            start_t=find(date_annual==yr_s); 
            end_t=find(date_annual==yr_e); 
               
            yc=PDO_x(start_t:end_t);
            y=zscore(yc);
            
            elseif type_nr==2  %%%% Monthly
                
                
                
                if  index_nr==1
                %%%% PDO Monthly 
                load('PDO_index_c'); % monthly
                
%                 wr_c=length(PDO_index_c(:,1));
%                 wr_c1=wr_c*12;
%                 wr_c2=[1:wr_c1]';
                
             %   80*12   % 1979
              %  115*12 -1 % 2014
              %  116*12 -1 % 2015
                
                PDO_index_c2=permute(PDO_index_c,[2, 1]);
                PDO_monthly_x=reshape(PDO_index_c2,[12*116 1]);
                Pwr=[1900:1/12:2016]';
                % y=zscore(PDO_monthly_x(960:1391));
                 y=zscore(PDO_monthly_x(949:1380));
                 
                %%%%%%%%%% 
                elseif index_nr==2 % IPO
                    
                 % IPO TPI filtered
                  %  Henley et al. 2015
                    load('C:\PHD\matlab_lib\Data\IPO_TPI_unfiltered.mat')


                    IPO_TPI_unfiltered_c=IPO_TPI_unfiltered(:,(2:13));
                    IPO_TPI_time=IPO_TPI_unfiltered(:,1);
%                     IPO_TPI_annual=nanmean(IPO_TPI_unfiltered_c,2);

                    [ac, bc, cc]=size(IPO_TPI_unfiltered_c);
                    IPO_TPI_unfiltered_c2=permute(IPO_TPI_unfiltered_c((1:146),:),[2,1]);
                    IPO_TPI_unfiltered_c3=reshape(IPO_TPI_unfiltered_c2,12*146,1);
                    Pwr=[1870:1/12:2016]';
                    y=zscore(IPO_TPI_unfiltered_c3(1309:1740));
                 
                 
                %%%%%%%%%
                end
            
            end
       
          PC_nr=999; % just to get a number for enso    
          letter='b';  
                
            label_2='°C yr^-^1';
            
            
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
   elseif iso_nr==11      %SAM Fogt        
     
       if sea_nr==1 && type_nr==1
       
            load('SAM_Fogt.mat')
            SAM_Fogt_c=SAM_Fogt(41:140,:);% SAM Fogt 1905-2004   9999 values at the start 
            SAM_annual_Fogt=mean(SAM_Fogt_c,2); 
            date_annual=[1905:2004]';     

            start_t=find(date_annual==yr_s); 
            end_t=find(date_annual==yr_e); 
               
            yc=SAM_annual_Fogt(start_t:end_t);
            y=zscore(yc);

      

       end

       
       
                 PC_nr=999; % just to get a number for enso    
                 letter='b';  
                
                label_2='°C s.d.^-^1';
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
%    elseif iso_nr==12      %SST Aus                 
%                 if sea_nr==1 && type_nr==1
%                 load('C:\PHD\matlab_storage_of_output_files\HadISST_SST_Aus_box.mat')  % annual SST south of Aus
%                 
%                 
%                         date_annual=[1870:2014]';
%                         start_t=find(date_annual==yr_s); 
%                         end_t=find(date_annual==yr_e); 
%                         
%                         yc=HadISST_M_sst_box4(start_t:end_t);
%                         y=zscore(yc);
%                 
%                 
%                 
%                 end
%                 
%                 
%                  PC_nr=999; % just to get a number for enso    
%                  letter='b';  
%                 
%                 label_2='°C s.d.^-^1';
                
    elseif iso_nr==14      %CP EP index             
       
        el_nino_type=1; %  (1) CP (2) EP (12) CP alt2 nino-3 removed
                
       if el_nino_type==1
           
%           load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_CP_lim-70_120_20_-20_1950-2014_annual_c.mat') % 1950
          load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_CP_lim-70_120_20_-20_1900-2014_annual_c.mat') % 1900  fact (-1) used
           
           
       elseif el_nino_type==2
           
           %load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_EP_lim-70_120_20_-20_1950-2014_annual_c.mat')
           load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_EP_lim-70_120_20_-20_1900-2014_annual_c.mat')
           
           
           
      elseif el_nino_type==12 % not in use
          
           load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_CP_alt2_lim-70_120_20_-20_1950-2014_annual.mat')
           
       end
       
                   time_c=MA_PCs_save_monthly(:,1);
                    yr_2=find( time_c==yr_e)+11;
                    yr_1=find( time_c==yr_s); 
       
                PC_nr=2;% leading PC
                fact=-1; % 1950 +1, 1900 -1
                y=MA_PCs_save_monthly((yr_1:yr_2),PC_nr)*fact;
       
                 PC_nr=999; % just to get a number for enso   
                 
                 if param_nr==1
                 letter='b';  
                 elseif param_nr==2
                 letter='b';      
                 end
                label_2='°C s.d.^-^1';
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    elseif iso_nr==15   % SIC PCs
        
        SIC_PCs_alt=5;%%%%%%%%%%%%%%%
        
    
          load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SIC_lim-150--30_-64_-75_1979-2014_annual.mat'); % use 
          
  
                    time_c=MA_PCs_save_monthly(:,1);
                    yr_2=find( time_c==yr_e)+11;
                    yr_1=find( time_c==yr_s); 
                %%%%%%%%%%%%%%%%%%%%%%    
                PC_nr=2;%  PC# (2) PC1 Leading SAM-related (3) PC2 PSA1-related (4) PC3
                %%%%%%%%%%%%%%%%%%%%%%%%
                
 

                 %%%%%%%%%%%%%%%%
                 if PC_nr==2
                 letter='a';  
%                  fact=-1; %
                  fact=1; %%%%%%%%%%%%% leave SIC PCs unchanged
                 
                 elseif PC_nr==3
                 letter='b'; 
%                  fact=-1; % 
                  fact=-1; %%%%%%%%%%%%% 
                
                 elseif PC_nr==4
                 letter='c'; 
                 fact=-1; %                 
                 
                 
                 end
                 
                y=MA_PCs_save_monthly((yr_1:yr_2),PC_nr)*fact;
                 
                label_2='SIC s.d.^-^1';
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    elseif iso_nr==16   % SST PCs tropical Pacific       
          

                    eof_sst_nr=2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if eof_sst_nr==1
    load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_lim-240--80_30_-30_1871-2014_annual_trend.mat'); %                      
    elseif eof_sst_nr==2
        
%     load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_lim-240--80_30_-30_1950-2014_annual_trend.mat'); %
    load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_lim-240--80_20_-20_1950-2014_annual_trend.mat'); %
    
    elseif eof_sst_nr==3  % 1979 
%     load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_lim-240--80_30_-30_1979-2014_annual_trend.mat'); %    
      load('C:\PHD\matlab_storage_of_output_files\HadISST_PCs_SST_lim-240--80_20_-20_1979-2014_annual_trend.mat'); % narrower
    end

           
             
                    time_c=MA_PCs_save_monthly(:,1);
                    yr_2=find( time_c==yr_e)+11;
                    yr_1=find( time_c==yr_s); 
                %%%%%%%%%%%%%%%%%%%%%%    
                PC_nr=4;%  PC# (2) PC1, (3) PC2, (4) PC3, (5) PC4
                %%%%%%%%%%%%%%%%%%%%%%%%
                
                
                lat1=-30;
                lat2=30;
                lon1=-240;
                lon2=-60;

                 %%%%%%%%%%%%%%%%
                 if PC_nr==2
                 letter='a';  
%                  fact=-1; %
                  fact=1; %%%%%%%%%%%%% leave SIC PCs unchanged
                 
                 elseif PC_nr==3
                 letter='c'; 
%                  fact=-1; % 
                  if yr_s==1950 
                    fact=-1; %%%%%%%%%%%%%
                  elseif yr_s==1979
                    fact=1;
                    
                  else
                    fact=1;  
                  end
                  
                 
                 elseif PC_nr==4
                 letter='e'; 
                 
                    if yr_s==1950 
                    fact=1; %
                    elseif yr_s==1979 
                    fact=-1; %    
                    end
                    
                 elseif PC_nr==5
                    letter='g';
                    fact=1; %
                 
                 end
                 
                y=MA_PCs_save_monthly((yr_1:yr_2),PC_nr)*fact;
                 
                label_2='°C s.d.^-^1';
                
                
                
   
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%     elseif iso_nr==17   % SST ERSST PCs tropical Pacific                
%     load('C:\PHD\matlab_storage_of_output_files\ERSST_PCs_SST_lim120-280_20_-20_1950-2014_annual_trend.mat'); %          
%           time_c=MA_PCs_save_monthly(:,1);
%                     yr_2=find( time_c==yr_e)+11;
%                     yr_1=find( time_c==yr_s); 
%                 %%%%%%%%%%%%%%%%%%%%%%    
%                 PC_nr=2;%  PC# (2)PC1 (3)PC2 (4)PC3
%                 %%%%%%%%%%%%%%%%%%%%%%%%
%                 
%                 
%                 lat1=-30;
%                 lat2=30;
%                 lon1=-240;
%                 lon2=-60;
% 
%                  %%%%%%%%%%%%%%%%
%                  if PC_nr==2
%                  letter='a';  
% %                  fact=-1; %
%                   fact=1; %%%%%%%%%%%%% leave SIC PCs unchanged
%                  
%                  elseif PC_nr==3
%                  letter='c'; 
% %                  fact=-1; % 
%                   if yr_s==1950 || yr_s==1979
%                     fact=1; %%%%%%%%%%%%%
%                   else
%                     fact=1;  
%                   end
%                   
%                  
%                  elseif PC_nr==4
%                  letter='e'; 
%                  fact=1; %                 
%                  
%                  
%                  end
%                  
%                 y=MA_PCs_save_monthly((yr_1:yr_2),PC_nr)*fact;
%                  
%                 label_2='°C s.d.^-^1';    
           
end


%% regression of deseasonal 

method_nr=2;  % Use 2 regstats it works fine

if type_nr==1 % interannual
    if  sea_nr==1 && param_nr==1 % SST
    X=HadISST_M_annual(:,:,1:end);
    
    elseif  sea_nr==1 && param_nr==2 % SIC extended winter season
    X=HadISST_M_season_reg_amjjas(:,:,1:end);   
    
    elseif sea_nr==2
    X=HadISST_M_season_reg_djf(:,:,1:end);  
    elseif sea_nr==3
    X=HadISST_M_season_reg_mam(:,:,1:end);  
    elseif sea_nr==4
    X=HadISST_M_season_reg_jja(:,:,1:end);  
    elseif sea_nr==5
    X=HadISST_M_season_reg_son(:,:,1:end); 
    end     

elseif type_nr==2
    
    if SST_dataset==1 % HadISST
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if  set_yr==1 &&  param_nr==1
        load('HadISST_sst_c4.mat') % from HadISST_cosweight_clim.m 1979-
        ntim=432;
    
    
        elseif (set_yr==6 &&  param_nr==1)  || (iso_nr==16 && param_nr==1) % 1950-2014
            
        yr_s_nr=1;
        
            if yr_s_nr==1
            load('C:\PHD\matlab_lib\Data\HadISST_sst_1950_2014_clim.mat')
            ntim=780;
    
            elseif yr_s_nr==2
            % 1871-2014
            load('C:\PHD\matlab_storage_of_output_files\clim_SST_1871_2014.mat');
            ntim=1728;
            HadISST_sst_c4=era_z500_c4;
        
            end
            
        
        elseif (set_yr==1 || set_yr==8  || set_yr==9) &&  param_nr==2
        load('C:\PHD\matlab_storage_of_output_files\HadISST_clim_sic_1979_2014.mat');
        ntim=432;    
            

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif SST_dataset==2 % ERSST
      
            load('C:\PHD\matlab_storage_of_output_files\clim_SST_ERSST_1950_2014.mat');
            ntim=780;
            HadISST_sst_c4=era_z500_c4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    
    

    
    end
    
    
    if SST_dataset==1
        nlat=180;
        nlon=360;
        
    elseif SST_dataset==2    
        nlat=89;
        nlon=180;
        
    end
    
    HadISST_sst_c5 = reshape(HadISST_sst_c4,ntim, nlat, nlon); 
    HadISST_sst_c6=permute(HadISST_sst_c5,[3 2 1]);
  %  X=HadISST_sst_c6(:,:,(yr_1:yr_2));
  
        if  sea_nr==1
            
            if set_yr==1 || set_yr==10
                
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%              if yr_e==2014   
%              X=HadISST_sst_c6(:,:,(1:end));
%              elseif yr_e==2013
%              X=HadISST_sst_c6(:,:,(1:end-12));    
%              end
             
             
             
             if yr_s==1871 || yr_s==1950
             X=HadISST_sst_c6(:,:,(1:end));
             elseif yr_s==1979 && param_nr==1
             X=HadISST_sst_c6(:,:,(349:end));
             elseif yr_s==1979 && param_nr==2
             X=HadISST_sst_c6(:,:,(1:end)); 
             
             
             end
             
             
             %%%%%%%%%%%%%%%
             
            elseif set_yr==8
                if yr_e==2014
                X=HadISST_sst_c6(:,:,(253:end));   % 2000-2014 IPO-
                elseif yr_e==2013
                X=HadISST_sst_c6(:,:,(253:end-12));   % 2000-2013 IPO-
                end
            elseif set_yr==9
             X=HadISST_sst_c6(:,:,(1:252));   % 2000-2014 IPO-
            
            end
        elseif  sea_nr>=2
%             999
            X=HadISST_sst_c6(:,:,(mm_in_c));
        end
        
end
    
    
[A B C]=size(X);
% b=zeros(2, A,B);bINT=zeros(2,2,A,B);r=zeros(324,A,B);rINT=zeros(324,2,A,B);
b=zeros(2, A,B);
bINT=zeros(2,2,A,B);
r=zeros(C,A,B); % residual
rINT=zeros(C,2,A,B);
STATS=zeros(1,4,A,B);  % 3 p-value

 p_all=zeros(A,B,1);  s_all=zeros(A,B,1);  % just one layer needed?
 
%  s_all_c_in=s_all;


for i=1:A
    for j=1:B
%         [b(:,i,j), bINT, r(:,i,j), rINT, STATS(:,:,i,j)]=regress(X, [squeeze(Y(i,j,:)) ones(size(X))]);

%%%%%%%%%%%%%%%%  % might need a matrix with ones to get rid of warning
%             if method_nr==1
%             [b(:,i,j), bINT, r(:,i,j), rINT, STATS(:,:,i,j)]=regress( [squeeze(X(i,j,:)) ], y);
%         
            
%%%%%%%%%%%%%%%%%%
%%% Replace NaNs
            if  isnan(X(i,j,:))==1
               X(i,j,:)=1000;     
%                 s_all_c_in(i,j)=1;
            end
  %%%%%%%%%%%%%%%%%%%%

        % Alt method
            if method_nr==2
            
                
            
            s=regstats(squeeze(X(i,j,:)),y,'linear','all');
            s.tstat;
            s_all(i,j,:)=s.tstat.beta(2); % beta- regresion coeficient 
            p_iso_rec=s.tstat.pval;
            p_all(i,j,:)= p_iso_rec(2); % The second element is the p-value for a significantly non-zero slope.
            
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%% Put (z500, u850) field first
            %%%%% check slop in at the end of file



    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Figure
%%%%%%%%%%%%%%%%%%%%%%%
%%% Monthly values

    

        proj_nr=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        proj= 'stereo';

   

    if  type_nr==2 &&  iso_nr==9 % TREND  per decade  % monthly
        
                if param_nr==1 % SST
                s_all_c1=s_all*10*12;
                elseif param_nr==2 % SIC in % /decade
                s_all_c1=s_all*10*12*100;
                end
    elseif type_nr==2 && (iso_nr==7  || iso_nr==8  || iso_nr==10) % Index regression 
        
    %1. Too high values using monthly PDO??
    %  UNIT NOT PER YEAR dont multiply by 12!!!!!!!!!!!!!!!
        s_all_c1=s_all;  % C per standard dev.
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Annual values
    elseif iso_nr==9 &&  type_nr==1% Trend per decade % annual
        s_all_c1=s_all*10;   
     
    else    % C per standard dev.
        s_all_c1=s_all;
    end

[c_max, c_min, c_limit]=Ma_max_min(s_all_c1);

c_max_c=c_max;

c_limit_c=c_limit;

    if lock_scalebar==1
         
         
        if  iso_nr==7 || iso_nr==15 % PCs
         % c_limit=0.5;
                if param_nr==1
                c_limit=1.0;
                elseif param_nr==2
                c_limit=0.1;
                end
                
                
           adj_nr=25;
         
        elseif  iso_nr==8 % nino
%         c_limit=1.5;
         c_limit=1.0;
          adj_nr=25;
          
        elseif iso_nr==16 || iso_nr==17
           
         %c_limit=1.2; % 1871
         c_limit=1.2;  % 1950
         adj_nr=0;
         
          
          
        elseif  iso_nr==6 % PPA

         c_limit=0.05;
          adj_nr=10;
          
          
        
        elseif  iso_nr==9 && param_nr==1
         
         c_limit=0.38;
         
        elseif  iso_nr==9 && param_nr==2  % SIC
         
         c_limit=0.15*100;
         
         elseif  iso_nr==10  || iso_nr==14% PDO 
             
        adj_nr=25;    
       % c_limit=0.5;
       if param_nr==1
        c_limit=1.0; 
       elseif param_nr==2
        c_limit=0.1; 
        
       end

        
        elseif iso_nr==11 || iso_nr==12% SAM Fogt
    
          c_limit=0.5;     
        
         end
         
    end
    

s_all_c=s_all_c1;

if  iso_nr==9
    
p_level_rc=1; 

else
p_level_rc=0.05;   
end

%p_level_rc=0.01;

s_all_c(p_all>=p_level_rc)=NaN; 


%%%%%%%%%%%%%
[wa, wb]=size(s_all_c);

% for i=1:wa
%     for j=1:wb
% 
% %             if  isnan(s_all_c(i,j,:))==1
% %                s_all_c(i,j)=NaN;          
% %             end
% 
% 
%          if s_all_c(i,j)>=1
%                     s_all_c(i,j)=NaN;
%          end
%          
%          if s_all_c(i,j)<=-1
%                     s_all_c(i,j)=NaN;
%          end
% 
%          
%          
%          
% %          s_all_c(360,161)
%          
%     end
% end
%%%%%%%%%%%%%

% wr_c=find(s_all_c>100 | s_all_c<0.001);

%   lon1=-250;
%    lon2= 110;


fig('units','inches','width',14,'height',9,'font','Helvetica','fontsize',16,'border','on');



if   iso_nr==9  
    
        if param_nr==1
        lat2=10;
        adj_nr=25;
        elseif param_nr==2
        lat2=-50; 
        adj_nr=9.5;
        end
        
        
elseif  iso_nr==10  
       lat2=40;   
elseif  iso_nr==14  || iso_nr==15
    
        if param_nr==1
        lat1=5;
        lat2=-30;  
    
        elseif param_nr==2
        lat1=-50;    
        lat2=-90; 
        adj_nr=0;
        end
    
    
    
    
end




%%%%%%%%%%%
if param_nr==1
   axesm( proj,'MapLatLimit',[lat1+adj_nr lat2],'MapLonLimit',[lon1 lon2],'Grid','on','ParallelLabel','on','Frame','on',... %
       'MeridianLabel','on','FontWeight','bold','FontSize',14,...
       'mlabellocation',[0:30:180,0:-30:-300]); 
   
elseif param_nr==2 % stereo
    if proj_nr==1
      axesm( proj,'MapLatLimit',[lat1 lat2],'Grid','on','ParallelLabel','on','Frame','on',... %
       'MeridianLabel','on','FontWeight','bold','FontSize',18,...
       'mlabellocation',[0:30:180,0:-30:-180]); 
   
    elseif proj_nr==2 % mrec
       axesm( proj,'MapLatLimit',[lat1+10 lat2+5],'MapLonLimit',[lon1 lon2],'Grid','on','ParallelLabel','on','Frame','on',... %
       'MeridianLabel','on','FontWeight','bold','FontSize',18,...
       'mlabellocation',[0:60:180,0:-60:-280]);   
   
   
    end
   
   
end
   %Frame - around figure
   
   
%        axesm('MapProjection','mercator','FLatLimit',[-80 20],'MapLonLimit',[lon1 lon2],'ParallelLabel','on','MeridianLabel','on',...
%      'FontWeight','bold','FontSize',14) 
  
   set(gca,'box','off','XColor',[1,1,1],'YColor',[1,1,1]);  % turns off border box and axes  

% pcolor(era_long,era_lat,squeeze(b(1,:,:))'); shading flat; colorbar;

%%%%
%hSurf=surfm(double(era_lat),double(era_long),squeeze(b(1,:,:))');

% s_all_c_in=logical(s_all_c_in);
%  s_all_c(s_all_c_in)=NaN;

wr_cc=squeeze(s_all_c(:,:,1))';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SST_dataset==2 % ERSST
    
    
    [ab, ac]=size(wr_cc);
    
    for i=1:ab
        for j=1:ac
           
            
%            if j<=91 % 'New Guinea', ,'North and South America'
%                 wr_d=landmask(HadISST_lat(i),HadISST_lon(j),'Australia');
%            else
%                 wr_d=landmask(HadISST_lat(i),360-HadISST_lon(j),'Australia');
%            end
           
           
%            if wr_d==1

           if wr_cc(i,j)==0
             
           wr_cc(i,j)=NaN;
           
           end
            
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

  
end

%%%%%%%%%%%%%%%%%%%%





hSurf=surfm(double(HadISST_lat),double(HadISST_lon),wr_cc);

 hold on
 colormap(b2r(-c_limit,c_limit));
%%%%%%%%

% ylim([-90 0]);
br=colorbar;

%where the position arguments are [xposition yposition width height].

% inflate(h,1.5)
        pos=get(br,'Position');
%         pos(1)=pos(1)-0.018;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar   (+) --->>   +0.02;
        pos(1)=pos(1)-0.0;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar   (+) --->>   +0.02; 
            if param_nr==1 &&  (iso_nr==16 || iso_nr==17)
                
                pos(1)=pos(1)-0.015;  
                pos(2)=pos(2)+0.072;  
%      
                pos(3)=pos(3)- 0.0;  % widthcolorbar
                pos(4)=pos(4)-0.15;  % height colorbar  
                
            elseif   param_nr==1 
                pos(2)=pos(2)+0.095;  
%             elseif  param_nr==1 && iso_nr ==9
%             pos(2)=pos(2)+0.095;     
                pos(3)=pos(3)- 0.0;  % widthcolorbar
                pos(4)=pos(4)-0.15;  % height colorbar
            elseif param_nr==2 && iso_nr==6
                pos(1)=pos(1)+0.011; 
                pos(2)=pos(2)+0.041; 
                pos(3)=pos(3)- 0.0;  % widthcolorbar
                pos(4)=pos(4)-0.15;  % height colorbar
                
            elseif param_nr==2 && (iso_nr==9  ||  iso_nr==14)
                
                if proj_nr==1
                pos(1)=pos(1)+0.005; 
                pos(2)=pos(2)+0.07; 
                pos(3)=pos(3)- 0.0;  % widthcolorbar
                pos(4)=pos(4)-0.15;  % height colorbar
                
                elseif proj_nr==2 % merch iso
                pos(1)=pos(1)-0.015; 
                pos(2)=pos(2)+0.055; 
                pos(3)=pos(3)- 0.0;  % widthcolorbar
                pos(4)=pos(4)-0.15;  % height colorbar               
                
                end
                
                
            elseif  iso_nr==7 ||  iso_nr==15 
                pos(1)=pos(1)+0.014; 
                pos(2)=pos(2)+0.170;
                pos(3)=pos(3)- 0.0;  % widthcolorbar
                pos(4)=pos(4)-0.35;  % height colorbar
                
            elseif param_nr==2
                pos(2)=pos(2)+0.045;   
                pos(3)=pos(3)- 0.0;  % widthcolorbar
                pos(4)=pos(4)-0.15;  % height colorbar   
            end

        
        set(br,'Position',pos)




% colormap(b2r(-c_limit,c_limit));

color_alt=2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (iso_nr==7 || iso_nr==15 || iso_nr==14) && param_nr==2
    color_alt=4; 
elseif  iso_nr==9  
    color_alt=1; 
end


if color_alt==1
    colormap(brewermap(256,'*RdBu'));
elseif color_alt==2
    colormap(flipud(cbrewer('div','Spectral',10))); 
elseif color_alt==3
    colormap(flipud(cbrewer('div','PiYG',20)));    
elseif color_alt==4    
    colormap(cbrewer('div','BrBG',20));    
    
end



% colb_x=-0.02; colb_y=0;
% move([colb_x, colb_y],br)
%             inflate(1,0.76,h);

    if iso_nr==7 || iso_nr==8 || iso_nr==9  || iso_nr==10  || iso_nr==12 || iso_nr==14  || iso_nr==6 ||  iso_nr==15 ||  iso_nr==16 ||  iso_nr==17
        co_pa=[.1 .1 .1];
        
        line_w_nr=2;
        if param_nr==1
        line_w_nr=1.5;
        end
        
    else
        co_pa=[.9 .9 .9];
        line_w_nr=2;
    end

    h1= contourm( double(HadISST_lat),double(HadISST_lon),squeeze(p_all(:,:,1))', [0.05],'--','ShowText','off',...
    'Linecolor',co_pa,'LineWidth', line_w_nr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coast_nr=0;

if coast_nr==1
    load coast
    %%%%%%%%%%%%%
    % to be able to use coastline for continents in combination with ant
    % coastline and grounding line from bedmap
in_c=find(lat<-60);
lat_cr=lat;
lon_cr=long;
lat_cr(in_c)=NaN;
lon_cr(in_c)=NaN;

% if param_nr==1 &&  iso_nr==16
%         color_code=[.6 .6 .6]; 
if param_nr==1
       color_code=[.1 .1 .1];
       
elseif param_nr==2 && param_nr==2
       color_code=[.0 .0 .0];        
       
elseif param_nr==2
       color_code=[.6 .6 .6];  
end

%%%%%%%%%%%%%
 plot3m(lat_cr,lon_cr,'-','LineWidth', 2, 'color',color_code);
%plot3m(lat(in_c),long(in_c),'.k','MarkerSize', 1);
%%%%%%%%%%%%%


addpath C:\PHD\matlab_mapping
bedmap2('gl','LineWidth', 1.5, 'color',color_code);
bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end



%%%%%%%%%%%%%%%%%%%%

%letter='b';

% label_size=16;
label_size=20;

% if cor_nr==4
%     a1=axestext_c(0.04, +0.05, ['(',label_1,')'] );
% else
%     a1=axestext_c(0.1, -0.04, ['(',label_1,' ', label_2,')'] );
% % text(-1.9,-1.8, [label_1,' ', label_2],'FontWeight','bold','FontSize',label_size,'rotation',90); 
% end
% set(a1,'FontWeight','bold','FontSize',label_size,'rotation',90);


% if   iso_nr==9  
%         label_across=-70;
% elseif  iso_nr==10  
%          label_across=-70;
% end

label_vert=-300;
if param_nr==1 && (iso_nr==16 || iso_nr==17)
     label_across=-28;
elseif param_nr==1
     label_across=-70;
elseif param_nr==2 && (iso_nr==6 || iso_nr==7 || iso_nr==14 ||  iso_nr==15)
    label_across=-55; 

elseif param_nr==2 && iso_nr==9    
    label_across=-55;
    label_vert=-298;
elseif param_nr==2
     label_across=-85;   
end




% end



if  iso_nr==10  
mlabel('MLabelParallel',label_across,'PLabelLocation',[-75 -60 -45 -30 -15  0 15 30],'PLabelMeridian',label_vert) 
else
mlabel('MLabelParallel',label_across,'PLabelLocation',[-75 -60 -45 -30 -15  0 15 ],'PLabelMeridian',label_vert) 
end



%  % save figure

% if cor_nr==1
%     t4=title(['Trend ', name,' ', season]);
% elseif cor_nr==3
%    t4= title(['Regression ', iso,'/', name]);
% elseif cor_nr==4
    
    if  iso_nr==7
                 
         label_1='°C s.d.^-^1'; % per nomralized unit standard deviation
         
         
         if param_nr==1
             
         x_lab=0.920;y_lab=0.7;
         elseif param_nr==2
          label_1='SIC';   
          %x_lab=1.04;y_lab=0.41;
          x_lab=1.09;y_lab=0.45;
          
          
         end
         
         
         
        if PC_nr==2 
         pc_str='SAM';

        elseif PC_nr==3 && iso_nr==7
            pc_str='PSA1';

        elseif PC_nr==4 && iso_nr==7
            
            if fact==1
            pc_str='PSA2';
            elseif fact==-1
            pc_str='PSA2';
            end
        end
        
      elseif   iso_nr==6
        
          pc_str='PPA';
          label_1='SIC'; 
          x_lab=1.04;y_lab=0.43;
          letter='';
          PC_nr=0;
        
      elseif   iso_nr==8
          
            if nino_nr==4
            pc_str='Niño-4';
            elseif nino_nr==34
            pc_str='Niño-3.4';
            elseif nino_nr==12
            pc_str='Niño-1+2'; 
            end
        label_1='°C s.d.^-^1'; 
        x_lab=0.920;y_lab=0.70;
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
      elseif   iso_nr==16 % SST PCs tropical Pacific
          
       
            pc_str='SST PC';

        label_1='°C s.d.^-^1'; 
        x_lab=0.920;y_lab=0.70;
          
           if PC_nr==2
              pc_str='HadISST EOF1';

                if yr_s==1871
                variance_exp=54; %1871          
              
                elseif yr_s==1950
                %variance_exp=52; %1950
                variance_exp=62; %1950 narrower region 20S-20N
                
                elseif yr_s==1979
%                 variance_exp=54; %1979
                variance_exp=64; %1979

                end
                
           elseif  PC_nr==3
              pc_str='HadISST EOF2';
              
              if yr_s==1871
                variance_exp=10; %1871
              elseif yr_s==1950
%                 variance_exp=11; %1950
                variance_exp=12; %1950 narrower region 20S-20N
              elseif yr_s==1979
%                 variance_exp=12; %1979   
                variance_exp=14; %1979 
              end
              
           elseif  PC_nr==4
              pc_str='HadISST EOF3'; 
              if yr_s==1950
              variance_exp=6; %1950 narrower region 20S-20N
              
              elseif yr_s==1979
              variance_exp=4; %1950 narrower region 20S-20N
              
              end
              
           elseif  PC_nr==5
              pc_str='HadISST EOF4'; 
              if yr_s==1950
              variance_exp=3; %1950 narrower region 20S-20N
              
              elseif yr_s==1979
              variance_exp=3; %1950 narrower region 20S-20N
              
              end
              
              
              
           end
           
   a2=axestext_c(0.98, +1.0-0.001, [num2str(variance_exp),'%'] );
   set(a2,'FontWeight','bold','FontSize',label_size);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
      elseif   iso_nr==17 % ERSST SST PCs tropical Pacific
          
       
            pc_str='SST PC';

        label_1='°C s.d.^-^1'; 
        x_lab=0.920;y_lab=0.70;
          
           if PC_nr==2
              pc_str='ERSST EOF1';

                if yr_s==1871
                variance_exp=54; %1871          
              
                elseif yr_s==1950
                %variance_exp=52; %1950
                variance_exp=53; %1950 narrower region 20S-20N
                
                elseif yr_s==1979
%                 variance_exp=54; %1979
                variance_exp=64; %1979

                end
                
           elseif  PC_nr==3
              pc_str='ERSST EOF2';
              
              if yr_s==1871
%                 variance_exp=10; %1871
              elseif yr_s==1950
%                 variance_exp=11; %1950
                variance_exp=11; %1950 narrower region 20S-20N
              elseif yr_s==1979
%                 variance_exp=12; %1979   
%                 variance_exp=14; %1979 
              end
              
           elseif  PC_nr==4
              pc_str='ERSST EOF3'; 
              if yr_s==1950
              variance_exp=10; %1950 narrower region 20S-20N
              elseif yr_s==1979
%               variance_exp=4; %1950 narrower region 20S-20N
              
              end
           end
           
   a2=axestext_c(0.98, +1.0-0.001, [num2str(variance_exp),'%'] );
   set(a2,'FontWeight','bold','FontSize',label_size);                                   
           
   
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
      elseif   iso_nr==9 ||  iso_nr==10
            if  iso_nr==9
                pc_str='Trend';
                
                if param_nr==1
              	label_1='°C decade^-^1';
                x_lab=0.920;y_lab=0.74;
                
                elseif param_nr==2
                label_1='% decade^-^1';
                
                    alt_cp1=2;%%%%
                    
                    if alt_cp1==1
                    x_lab=0.920;y_lab=0.64;
                    elseif alt_cp1==2
                    x_lab=1.075;y_lab=0.34;
                    elseif alt_cp1==3
                    x_lab=1.025;y_lab=0.10;
              
                    
                    end
                    
                end
                
           
            elseif  iso_nr==10
                
                        if     index_nr==1
                        pc_str='PDO';  
                        elseif  index_nr==2
                        pc_str='IPO';  
                        end
                       
                        label_1='°C s.d.^-^1'; % per nomralized unit standard deviation
                        
      x_lab=0.920;y_lab=0.67;      
                        
          
            end
      
  
                   elseif  iso_nr==11
                        pc_str='SAM Fogt';  
                        label_1='°C s.d.^-^1'; % per nomralized unit standard deviation     
                     x_lab=0.920;y_lab=0.67;    
                     
                     
                    elseif  iso_nr==12
                        pc_str='SST Aus';  
                        label_1='°C s.d.^-^1'; % per nomralized unit standard deviation     
                        x_lab=0.920;y_lab=0.67;
                        
                    elseif  iso_nr==14
                            if el_nino_type==1
                            pc_str='CP';
                            elseif el_nino_type==2
                            pc_str='EP';
                            elseif el_nino_type==12
                            pc_str='CP alt 2';
                            end
                            
                            
                               if param_nr==1
                                label_1='°C s.d.^-^1'; % per nomralized unit standard deviation     
                                x_lab=0.920;y_lab=0.67;  
                        
                              elseif param_nr==2
                                label_1='SIC';   
                                %x_lab=1.04;y_lab=0.42;
                                
                                x_lab=1.1;y_lab=0.44;
                                
                               end
                              
                    elseif  iso_nr==15  
                        label_size=26;
                        if PC_nr==2
                        pc_str='PC1';
                                    if SIC_PCs_alt==1
                                    variance_exp=17;
                                    elseif SIC_PCs_alt==2
                                    variance_exp=17;
                                    elseif SIC_PCs_alt==5
                                    variance_exp=17;
                                    
                                    elseif SIC_PCs_alt==4
                                    variance_exp=15;
                                    
                                   end
                                    a2=axestext_c(0.95, +0.96, [num2str(variance_exp),'%'] );
                                    set(a2,'FontWeight','bold','FontSize',label_size);
                        
                                    
%                                             lon_b_c=[-150 -150 -140 -130 -120 -110 -100 -90 -80 -70 -60 -50 -40 -30 -30 -150];
%                                             lat_b_c=[ -90 -65   -65  -65  -65  -65  -65 -65 -65 -65 -65 -65 -65 -65 -90 -90];
%                                             plotm(lat_b_c,lon_b_c,'-','LineWidth',3,'Color',[1.0 0.0 .4])  
                                            
                                            
                                            lon_b_c=[-30 -40 -50 -60 -70 -80 -90 -100 -110 -120 -130 -140 -150   -150 -140 -130 -120 -110 -100 -90 -80 -70 -60 -50 -40 -30 -30 ];
                                            lat_b_c=[-75 -75 -75 -75 -75 -75 -75  -75  -75  -75  -75  -75  -75   -64   -64  -64  -64  -64  -64 -64 -64 -64 -64 -64 -64 -64 -75 ];
                                            plotm(lat_b_c,lon_b_c,'-','LineWidth',3,'Color',[1.0 0.0 .4])  
                                            
                                            
                                            
                                            
                        
                        
                        elseif PC_nr==3
                        pc_str='PC2';
                                     if SIC_PCs_alt==1
                                    variance_exp=11;
                                     elseif SIC_PCs_alt==2
                                    variance_exp=11; 
                                    elseif SIC_PCs_alt==5  %%%
                                    variance_exp=12;
                                    elseif SIC_PCs_alt==4
                                    variance_exp=10;
                                     end
                                    a2=axestext_c(0.95, +0.96, [num2str(variance_exp),'%'] );
                                    set(a2,'FontWeight','bold','FontSize',label_size);
                                    
                       elseif PC_nr==4
                        pc_str='PC3';
                                     if SIC_PCs_alt==1
                                    variance_exp=11;
                                     elseif SIC_PCs_alt==2
                                    variance_exp=8; 
                                    elseif SIC_PCs_alt==5 %%%%
                                    variance_exp=8;
                                    elseif SIC_PCs_alt==4
                                    variance_exp=10;
                                     end
                                    a2=axestext_c(0.95, +0.96, [num2str(variance_exp),'%'] );
                                    set(a2,'FontWeight','bold','FontSize',label_size);              
                        
                        
                        end
                               
                               if param_nr==2
%                                 label_1='SIC s.d.^-^1';   
                                label_1='SIC';   
%                                 x_lab=1.04; y_lab=0.41;
                               x_lab=1.09;y_lab=0.45;
                               end      
            
      end
    
%     if type_nr==2
%         season='monthly';
%     end

    
    
%    t4=title([ pc_str  ,' ',season  ]);


if type_nr==1
type_str=[];
elseif type_nr==2 && sea_nr==1
    season=[];
    type_str='monthly';
elseif type_nr==2
type_str='monthly';
end




      if param_nr==1  
      t4=title([ pc_str  ,' SST regression ',season,' ',type_str,' ',num2str(yr_s),'-',num2str(yr_e) ]);
      elseif param_nr==2 %&& ( PC_nr==2 || PC_nr==4 || PC_nr==999)
     t4=title([ pc_str  ,' SIC regression ',season,' ',type_str,' ',num2str(yr_s),'-',num2str(yr_e) ]);   
%      elseif param_nr==2 && (PC_nr==3 )%|| PC_nr==4)
%      t4=title([ pc_str  ,'(-1) SIC regression ',season,' ',type_str,' ',num2str(yr_s),'-',num2str(yr_e) ]); 

      end
      
 if type_nr==1
type_str='annaul';     
 end    
      
      
% set(t4,'FontWeight','bold','FontSize',20)

        pos=get(t4,'Position');
        
        if iso_nr==10 || iso_nr==9  
%         pos(2)=pos(2)-0.24;
        
              pos(2)=pos(2)-0.06;
%         elseif iso_nr==9 && param_nr==2  
%         pos(2)=pos(2)-0.15; % -0.35;
        elseif  (iso_nr==16 || iso_nr==17)  && param_nr==1
        pos(2)=pos(2)-0.06;

        elseif iso_nr==6 || (iso_nr==7 || iso_nr==14  || iso_nr==15 ) && param_nr==2  
        pos(2)=pos(2)-0.06;
%         elseif iso_nr==9 && param_nr==2 
%             
%         pos(2)=pos(2)-0.28;  
        
        else
        pos(2)=pos(2)-0.16;  
        end
        
        set(t4,'Position',pos,'FontSize',22, 'FontWeight', 'bold');    


        
        
  
        if iso_nr==10 
            letter_pos= [260 615 50 50]; 
        elseif (iso_nr==7 || iso_nr==8 || iso_nr==11 || iso_nr==12  || iso_nr==14 || iso_nr==6) && param_nr==1 
            letter_pos= [260 595 50 50];
        
        elseif  iso_nr==9 && param_nr==2
            letter_pos= [350 685 50 50];
        
        elseif  (iso_nr==16 || iso_nr==17) && param_nr==1
            letter_pos= [225 570 50 50];
        
        
        elseif  ( iso_nr==7 || iso_nr==14|| iso_nr==15 ) && param_nr==2 
            
            sic_alt=2;%%%%%%%%%%%%%%%%%%%%%%%
            
           if  sic_alt==1
            letter_pos= [270 545 50 50];
           elseif sic_alt==2
            letter_pos= [380 675 50 50];
           end
        
        else
         letter_pos= [270 560 50 50];   
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        letter_size=28;
        
        if param_nr==2
         letter_size=36;    
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     letter='d'

          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',letter_size ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
     

     
     
     if param_nr==2 && ( iso_nr==7 ||  iso_nr==15)
         
           a1=axestext_c(x_lab,y_lab, [' '] );
         
     else
         
         a1=axestext_c(x_lab,y_lab, ['(',label_1,')'] );
         
     end
     
     
     if (iso_nr==9 || iso_nr==14)
     
         if proj_nr==1
         label_size=26;
         elseif proj_nr==2
         label_size=20;    
         end
         
     end
     
% text(-1.9,-1.8, [label_1,' ', label_2],'FontWeight','bold','FontSize',label_size,'rotation',90); 

set(a1,'FontWeight','bold','FontSize',label_size,'rotation',90);
     
    if show_colorbar==1
    set(br, 'FontSize',18, 'FontWeight', 'bold'); 
    end    
    %%%%%%%%%%%%
    
     gridm('GLineStyle','--','Gcolor',[.1 .1 .1],'GLineWidth',1.5,...
    'MLineLimit',[lat1 lat2],'Galtitude', .02)


% %%%%%%%%%%%%areas of PSA1 and PSA2 overlap%%%%%%%%%%%%%%%%%%%%%%%%
% folder_c='C:\PHD\matlab_storage_of_output_files\';
% load([folder_c,'ind_psa_sum']);
% 
%           color_code=[1 .7 .1];                  
%     h2= contourm( double(HadISST_lat),double(HadISST_lon),ind_psa_sum',[1,1],'-k','ShowText','off',...
%     'Linecolor',color_code,'LineWidth', 2.5);   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSA2_box=0;
eRS_box=1;


if PC_nr==4 && PSA2_box==1 % SST box south of Aus
    
        lon_b_c=[120 120 155 155 120];
        lat_b_c=[ -45 -36 -36 -45 -45];
        plotm(lat_b_c,lon_b_c,'-','LineWidth',3,'Color',[1 .7 .1])  
        
        
elseif iso_nr==9 && param_nr==2 && eRS_box==1 
        
        lon_b_c=[-140 -140 -166 -166 -140];
        lat_b_c=[ -74 -70 -70 -74 -74];
        plotm(lat_b_c,lon_b_c,'-','LineWidth',3,'Color',[0.9 .1 .9])       
        
    
end





    %%%%%%%%%%%%%

if (PC_nr==2 || PC_nr==3 || PC_nr==4) && iso_nr==7 && param_nr==1
    
filename=['HadISST_regression_',name,'_',pc_str,'_',season,'_',type_str];


elseif iso_nr==11 || iso_nr==12  || iso_nr==14 || iso_nr==15 

filename=['HadISST_regression_',name,'_',pc_str,'_',season,'_',type_str,'_',num2str(yr_s),'-',num2str(yr_e)];
elseif SST_dataset==1
    if iso_nr==9
        iso='trend';
    end
% filename=['HadISST_regression_',name,'_',iso,'_',season];

filename=['HadISST_regression_',name,'_',pc_str,'_',season,'_',type_str,'_',num2str(yr_s),'-',num2str(yr_e)];

elseif SST_dataset==2
filename=['ERSST_regression_',name,'_',pc_str,'_',season,'_',type_str,'_',num2str(yr_s),'-',num2str(yr_e)];

end
   % num2str(round(date_annual(start_t))),'_',num2str(round(date_annual(end_t))),proj];
% 
% % filedir ='rice_isotopes\annual\';
filedir ='C:\PHD\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);

 % save as png 
orient landscape
% increase quality of picture by increase to -r500

quality_level_nr=3;%%%%%%%%%%%%%%%%%%%%%%%%%

if quality_level_nr==1
   quality_level_str= '-r110';
elseif quality_level_nr==2
  % quality_level_str= '-r150';
   quality_level_str= '-r190';
elseif quality_level_nr==3   
   quality_level_str= '-r240';
end

  export_fig('-png','-nocrop','-painters', '-depsc','-opengl', quality_level_str, savefilename_c); % PNG-nocrop'

%% for saving PSA pattern surfaces (step 1)-= one more step=- 
  
% save regression surface  
folder_c='C:\PHD\matlab_storage_of_output_files\';
% if sea_nr==1
    savefilename =[folder_c,filename,'_',num2str(p_level_rc),'_c.mat'];
% elseif sea_nr>=2    
%      savefilename =[folder_c,filename,'_',num2str(p_level_rc),'_',season,'.mat'];   
% end

 %save(savefilename,'s_all_c'); 
 
 save(savefilename,'s_all_c','s_all'); 
 
 
%% Where do PSA1 and PSA2 overlap
% PSA2
p_level_rc=0.05;
%p_level_rc=0.01;
 type_nr=2; % (2) monthly
 
 sea_nr=1;
 
 folder_c='C:\PHD\matlab_storage_of_output_files\';
 %%%%%%%%%%%%%%%
    if sea_nr==1
        season='annual';
    elseif sea_nr==2
        season='DJF';
    elseif sea_nr==3
        season='MAM';
    elseif sea_nr==4
        season='JJA';
    elseif sea_nr==5
        season='SON';
    end
%%%%%%%%%%%%%%
if type_nr==1
type_str=[];

elseif type_nr==2
type_str='monthly_';

end


%%%%%%%%%%%%%%%%%%%
% if sea_nr==1
%     load([folder_c,'HadISST_regression_SST_PSA2_monthly','_',num2str(p_level_rc),'.mat']);
if sea_nr>=3        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      load([folder_c,'HadISST_regression_SST_PSA2','_', season,'_',type_str,num2str(p_level_rc),'.mat']);  
elseif sea_nr==2   % DJF
      load([folder_c,'HadISST_regression_SST_PSA1','_', season,'_',type_str,num2str(p_level_rc),'.mat']); 
   %   load([folder_c,'HadISST_regression_SST_PSA2','_', season,'_',type_str,num2str(p_level_rc),'.mat']); 
elseif sea_nr==1 % annual
      season=[];
    load([folder_c,'HadISST_regression_SST_PSA2','_', season,'_',type_str,num2str(p_level_rc),'.mat']); 
      
end
    

s_all_c_psa2=s_all_c;  % this is the p-level surface as rs with lower p-levels are masked above

%%%%%%%%%%%%
 ind_psa2_nonnan=~isnan(s_all_c_psa2);
%%%%%%%%%%

s_all_c_psa2(find(s_all_c_psa2<0))=-1;  % Negative values -1
s_all_c_psa2(find(s_all_c_psa2>0))=1;  % Negative values -1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSA1

% if sea_nr==1
%     load([folder_c,'HadISST_regression _SST_PSA1_monthly','_',num2str(p_level_rc),'.mat']);
if sea_nr>=3
      load([folder_c,'HadISST_regression_SST_PSA1','_', season,'_',type_str,num2str(p_level_rc),'.mat']);  
elseif sea_nr==2 % DJF
        load([folder_c,'HadISST_regression_SST_SAM','_', season,'_',type_str,num2str(p_level_rc),'.mat']);  
      %   load([folder_c,'HadISST_regression_SST_PSA1','_', season,'_',type_str,num2str(p_level_rc),'.mat']);         
elseif sea_nr==1
              season=[];
    load([folder_c,'HadISST_regression_SST_PSA1','_', season,'_',type_str,num2str(p_level_rc),'.mat']); 
end

s_all_c_psa1=s_all_c;

%%%%%%%%%%%%%%%%%
 ind_psa1_nonnan=~isnan(s_all_c_psa1);
%%%%%%%%%%%%%%%%%

s_all_c_psa1(find(s_all_c_psa1<0))=-1;  % Negative values -1
s_all_c_psa1(find(s_all_c_psa1>0))=1;  % Negative values -1





% [~,ind_test_c]=ismember(s_all_c_psa1,s_all_c_psa2); % is not member of 

%%%%%%%%%%%%%%
 ind_psa_sum_nonactive=ind_psa1_nonnan+ind_psa2_nonnan;
%%%%%%%%%%%%%

 ind_psa_sum=s_all_c_psa1+s_all_c_psa2;
 ind_psa_sum_c=ind_psa_sum;
 
ind_nan= isnan(ind_psa_sum_c);
 
ind_psa_sum_c(ind_nan)=9999; 
 
%%

if sea_nr==1
filename=['ind_psa_sum_',num2str(p_level_rc),'.mat'];
elseif sea_nr>=2
filename=['ind_psa_sum_',season,'_',num2str(p_level_rc),'.mat'];
end

 savefilename =[folder_c,filename]; 
% save(savefilename,'ind_psa_sum'); 
 save(savefilename,'ind_psa_sum_c','ind_psa_sum_nonactive'); 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% IPO SAM version (use)
% Where do SAM and PSA1 overlap, SST and SIC
% PSA2
p_level_rc=0.05;
%p_level_rc=0.01;
 type_nr=2; % (2) monthly
 
 yr_s=1979;yr_e=2014;
 
 sign_nr=1; % (1) SAM*(-1)+PSA1, (2) SAM +PSA1
 
 if sign_nr==1
     fact_nr=-1;
 elseif sign_nr==2
     fact_nr=1;
 end
 
 
 sea_nr=1;
 
 folder_c='C:\PHD\matlab_storage_of_output_files\';
 %%%%%%%%%%%%%%%

        season='annual';

%%%%%%%%%%%%%%
if type_nr==1
type_str=[];

elseif type_nr==2
type_str='monthly_';

end

param_nr=2;% (1-SST/ 2 SIC)


if param_nr==1
param_str='SST';    
elseif param_nr==2
    
alt_pc=2;%(1) sam psa1 (2) pcs SIC %%%%%%%    
param_str='SIC';

end

%%%%%%%%%%%%%%%%%%%
% if sea_nr==1
%     load([folder_c,'HadISST_regression_SST_PSA2_monthly','_',num2str(p_level_rc),'.mat']);
% if sea_nr>=3        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       load([folder_c,'HadISST_regression_SST_PSA2','_', season,'_',type_str,num2str(p_level_rc),'.mat']);  
% elseif sea_nr==2   % DJF
%       load([folder_c,'HadISST_regression_SST_PSA1','_', season,'_',type_str,num2str(p_level_rc),'.mat']); 
%    %   load([folder_c,'HadISST_regression_SST_PSA2','_', season,'_',type_str,num2str(p_level_rc),'.mat']); 
% elseif sea_nr==1 % annual
      season=[];
      
      if param_nr==1
        load([folder_c,'HadISST_regression_',param_str,'_SAM','_', season,'_',type_str,num2str(p_level_rc),'.mat']);
      elseif param_nr==2
          if alt_pc==1
        load([folder_c,'HadISST_regression_',param_str,'_SAM','_', season,'_',type_str,num2str(yr_s),'-',num2str(yr_e),'_',num2str(p_level_rc),'.mat']);
          elseif alt_pc==2
        load([folder_c,'HadISST_regression_',param_str,'_PC1','_', season,'_',type_str,num2str(yr_s),'-',num2str(yr_e),'_',num2str(p_level_rc),'_c.mat']);  
          end
       end
      
% end
%     

s_all_c_sam=s_all_c;  % this is the p-level surface as rs with lower p-levels are masked above

%%%%%%%%%%%%
 ind_sam_nonnan=~isnan(s_all_c_sam);
%%%%%%%%%%

s_all_c_sam(find(s_all_c_sam<0))=-1;  % Negative values -1
s_all_c_sam(find(s_all_c_sam>0))=1;  %
s_all_c_sam(find(s_all_c_sam==0))=NaN; % % remove zeros outside of the sea ice zone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSA1

% % if sea_nr==1
% %     load([folder_c,'HadISST_regression _SST_PSA1_monthly','_',num2str(p_level_rc),'.mat']);
% if sea_nr>=3
%       load([folder_c,'HadISST_regression_SST_PSA1','_', season,'_',type_str,num2str(p_level_rc),'.mat']);  
% elseif sea_nr==2 % DJF
%         load([folder_c,'HadISST_regression_SST_SAM','_', season,'_',type_str,num2str(p_level_rc),'.mat']);  
%       %   load([folder_c,'HadISST_regression_SST_PSA1','_', season,'_',type_str,num2str(p_level_rc),'.mat']);         
% elseif sea_nr==1
              season=[];
              
              if param_nr==1
                load([folder_c,'HadISST_regression_SST_PSA1','_', season,'_',type_str,num2str(p_level_rc),'.mat']); 
              elseif param_nr==2
                if alt_pc==1
                load([folder_c,'HadISST_regression_',param_str,'_PSA1','_', season,'_',type_str,num2str(yr_s),'-',num2str(yr_e),'_',num2str(p_level_rc),'.mat']);
                elseif alt_pc==2
                load([folder_c,'HadISST_regression_',param_str,'_PC2','_', season,'_',type_str,num2str(yr_s),'-',num2str(yr_e),'_',num2str(p_level_rc),'_c.mat']);     
                end
              end
    
    
    
% end

s_all_c_psa1=s_all_c;

%%%%%%%%%%%%%%%%%
 ind_psa1_nonnan=~isnan(s_all_c_psa1);
%%%%%%%%%%%%%%%%%

s_all_c_psa1(find(s_all_c_psa1<0))=-1;  % Negative values -1
s_all_c_psa1(find(s_all_c_psa1>0))=1;  % Negative values -1
s_all_c_psa1(find(s_all_c_psa1==0))=NaN;




% [~,ind_test_c]=ismember(s_all_c_psa1,s_all_c_psa2); % is not member of 

%%%%%%%%%%%%%%
 ind_psa_sum_nonactive=ind_sam_nonnan+ind_psa1_nonnan;
%%%%%%%%%%%%%
    IPO_nr=2; % (1) IPO+, (2) IPO -  % just a name for SIC
    
    if IPO_nr==1
        fac_i=-1;
        ipo_str='ipo_p_';
    elseif IPO_nr==2
         fac_i=1;
         ipo_str='ipo_n_';
    end

 ind_psa_sum=s_all_c_sam*( fac_i)+s_all_c_psa1; % multiply by -1 negative corr
 
 
 ind_psa_sum_c=ind_psa_sum;
 
ind_nan= isnan(ind_psa_sum_c);
 
ind_psa_sum_c(ind_nan)=9999; 
 
%% save

if sea_nr==1
    if param_nr==1
    filename=['ind_sam_psa1_sum_',ipo_str,num2str(p_level_rc),'.mat'];
    elseif param_nr==2
        if alt_pc==1    
        filename=['ind_',param_str,'sam_psa1_sum_',ipo_str,num2str(p_level_rc),'.mat'];
        elseif alt_pc==2 
        filename=['ind_',param_str,'pc1_pc2_sum_',ipo_str,num2str(p_level_rc),'.mat'];
        end
    end

elseif sea_nr>=2
filename=['ind_sam_psa1_sum_',season,'_',num2str(p_level_rc),'_c.mat'];
end

 savefilename =[folder_c,filename]; 
% save(savefilename,'ind_psa_sum'); 
 save(savefilename,'ind_psa_sum_c','ind_psa_sum_nonactive'); 
 
 