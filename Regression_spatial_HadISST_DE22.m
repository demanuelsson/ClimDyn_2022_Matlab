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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('[\bRead in data ]\b')
tic
% if iso_nr==7 %***************
% iso='PCs';
% elseif iso_nr==9 %***************
% iso=[]; % just trend
% elseif iso_nr==14 %***************
% iso='CP EP index';
% elseif iso_nr==15 %***************
% iso='SIC PCs';
% end
 
show_maximum_point=0;
 
show_title=1; % (1/0)
area_2_box=0 ; % 0/1 for RICE SST corr Box and text for Area 2
  
 
  lock_scalebar=1;
  
 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SST_dataset=1; % (1) HadISST 

if SST_dataset==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HadISST data
% version 1.1 Rayner et al. 2003


% addpath C:\PHD\HadISST_data\ % Ienovo
%    ncdisp('HadISST_ice_c.nc');
%    ncdisp('HadISST_sst.nc');

 HadISST_time=ncread([data_dir,'HadISST_sst_c.nc'],'time'); %units         = 'days since 1870-1-1 0:0:0'
 % missing values  missing_value = -1e+30 ( looks more like -1000)\
 %monthly values
 
 HadISST_time_c=HadISST_time+datenum(1870,1,1);  % to 2015 month 7
 
 date_vec=datevec(double(HadISST_time_c)); 
 yyyy=date_vec(:,1); 
 mm=date_vec(:,2);
 HadISST_year_num=yyyy+(mm/12-1/12); 

 HadISST_time_bnds =ncread([data_dir,'HadISST_sst_c.nc'],'time_bnds');
 HadISST_lat=ncread([data_dir,'HadISST_sst_c.nc'],'latitude');
 HadISST_lon=ncread([data_dir,'HadISST_sst_c.nc'],'longitude');
 
 
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

   HadISST_ice=ncread([data_dir,'HadISST_ice_c.nc'],'sic'); % 360x180x1728 lon x lat x time
   M=HadISST_ice;
   
 elseif strcmp(name,'SST')==1 
   M=ncread([data_dir,'HadISST_sst_c.nc'],'sst'); % 360x180x1728 lon x lat x time
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

% %%%%%%%%%%%%%%%%%%%%%%%
% mm_in=mm(HadISST_start:HadISST_count);  % seasonal indexing. Only used for monthly reg [not used]
% 
%  if sea_nr==1
%     mm_in_c=[HadISST_start:HadISST_count]';
%     
% %  elseif sea_nr==2  % DJF
% %      mm_in_c=[find(mm_in==12);  find(mm_in==1); find(mm_in==2)];
% %      mm_in_c=sort( mm_in_c);
% %      mm_in_c=mm_in_c(3:end-1); % dont use JF for first year as there is no D
% %      
% % elseif sea_nr==3  % MAM
% %      mm_in_c=[find(mm_in==3);  find(mm_in==4); find(mm_in==5)];
% %      mm_in_c=sort( mm_in_c);
% %      mm_in_c=mm_in_c(4:end);
% %      
% % elseif sea_nr==4  % JJA
% %      mm_in_c=[find(mm_in==6);  find(mm_in==7); find(mm_in==8)];
% %      mm_in_c=sort( mm_in_c);
% %      mm_in_c=mm_in_c(4:end);
% %      
% %  elseif sea_nr==5  % SON
% %      mm_in_c=[find(mm_in==9);  find(mm_in==10); find(mm_in==11)];
% %      mm_in_c=sort( mm_in_c);
% %      mm_in_c=mm_in_c(4:end);
% 
% elseif sea_nr==7  % AMJJASON
%      mm_in_c=[find(mm_in==4);  find(mm_in==5); find(mm_in==6); find(mm_in==7); find(mm_in==8); find(mm_in==9); find(mm_in==10); find(mm_in==11)];
%      mm_in_c=sort( mm_in_c);
%      mm_in_c=mm_in_c(1:end);
%       
%  end
% %%%%%%%%%%%%%%%%%%%
%%
 tic

[nlon, nlat, ntim]= size(M(:,:,:)); %
HadISST_c = permute(M(:,:,:), [3 2 1]);


%%% -=(1)=-
% Remove monthly climatology for all cells 

 HadISST_c2 = double(reshape(HadISST_c, ntim, nlat*nlon));  % One time series (column) for each grid point
 
 %  lat 1, 2, 3, 4,..........33...lat 1, 2, 3, 4.....33       until 33X240=7920
 %  
 %1
 %time
 
 [HadISST_c3,clim_z500] = annave(HadISST_c2);   % Removes seasonal cycle

 HadISST_c4= reshape(HadISST_c3,  ntim, nlat, nlon );

 
%%% -=(2)=- 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Since the grid size decreases as you move towards the pole, 
% weight each grid box (i.e., multiply the time series at each grid box) 
% by the square root of the cosine of latitude (the weights are based on 
% the square root of the cosine so that the covariance matrix is weighted by the cosine of latitude).

%  z = cosweight(z, lat);
% weighted to account for the decrease in area towards the pole.
%era_z500_c2=cosweight(era_z500_c,era_lat_c); % Original UoW function    (time x lat x lon) 
HadISST_c5=cosweight(HadISST_c4,HadISST_lat);
%%%%%%%%%%%

% Back to old format again
HadISST_c6= reshape(HadISST_c5,  ntim, nlat, nlon );
M = permute(HadISST_c6, [3 2 1]);  %
toc
%%

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
toc
%%
fprintf('[\bSeasonal anomaly fields ]\b')
tic

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

          HadISST_M_season_amjjas=nanmean(dummy1(3:8,:),1);  %3. seasonal 
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
% %     %AMJJASON
elseif strcmp(season,'AMJJASON')==1 % [Used]
          HadISST_M_season_amjjason=nanmean(dummy1(4:11,:),1);  % 4.seasonal 
          HadISST_M_season_amjjason=HadISST_M_season_amjjason(1:end);
          HadISST_M_season_reg_amjjason(i,j,:)=HadISST_M_season_amjjason; % means for regression          
          HadISST_M_season_dummy_amjjason=detrend(HadISST_M_season_amjjason); 
          HadISST_M_season_detrend_amjjason(i,j,:)=HadISST_M_season_dummy_amjjason;            
              
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
fprintf('[\bLoad Z500 PCs ]\b')
tic
if iso_nr==7  % ZPCs

load([filedir,'ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2011_AMJJASON_c52_annual_mean_varimax_c2.mat']); % 2011
        if type_nr==1 
            time_c=MA_PCs_save(:,1);
            yr_2=find( time_c==yr_e);
           
        elseif type_nr==2
            time_c=MA_PCs_save_monthly(:,1);
            yr_2=find( time_c==yr_e)+11;
        end
            
     yr_1=find( time_c==yr_s); 
  
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
           % PC_nr=2; % (2) SAM (3) PSA1 (4) PSA2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
             label_2='';
             iso='z500';
   
            if PC_nr==2 %SAM. PC1 (but first coloumn in MA is time)
                                          
                if sea_nr==1  || sea_nr==7
                    if param_nr==1
                        fact=1;    
                    elseif param_nr==2  %SIC
%                     fact=1;%%%%%%%%%%%%%%%%%change polarity SAM-%%%%%%%%%%%%%%%%%%%%%%%
                        % fact=-1; % to corespond to SIC PC regression fig.
                        
                         fact=1; % to corespond to SIC PC regression fig. 2014
                    end  
                end
                
            
                
          
                
                
            elseif PC_nr==3
                
                        
                
                if sea_nr==1  || sea_nr==7
                   if param_nr==1 

                    fact=1;   %PSA1+
                   elseif param_nr==2

                     fact=1;   %PSA1   % to corespond to SIC PC regression fig. 2011
                   end 
                end
                 

                 
            elseif PC_nr==4
                 
                                      % Changed sign convetion            
                
                % use -1 because wave train eminates from Aus positive SSTA?
                % fact=1; %PSA2-  
                
                 if sea_nr==1 || sea_nr==7
                     if param_nr==1

                        fact=1; %PSA2+ 
                     elseif param_nr==2

                         fact=1; % to corespond to SIC PC regression fig.2011
                         
                     end
                 end
                
     
                
                % use as anticyclone in the Ross Sea likely to be associated with enriched isotopes. 
            end
            
            
            
         if type_nr==1   && sea_nr==1  % annual
            y=MA_PCs_save((yr_1:yr_2),PC_nr)*fact;
            
              
         elseif  type_nr==1  && sea_nr>=2  % seasonal one value per year
             
          
            y_c2=MA_PCs_save(:,PC_nr)*fact;
            y= y_c2; 
             
            %%%%%%%%%%%%%%%%%%%%%%% 
            
         elseif type_nr==2 && sea_nr==1 %  monthly
            y=MA_PCs_save_monthly((yr_1:yr_2),PC_nr)*fact;
            
           elseif type_nr==2 && sea_nr>=2 %  monthly seasonal 
               
             y=MA_PCs_save_monthly((mm_in_c),PC_nr)*fact;   
            
         end                
end
toc
%% regression of deseasonal 
fprintf('[\bRegression ... ]\b')
tic
method_nr=2;  % Use 2 regstats it works fine

if type_nr==1 % Interannual, one value yearly for seasonal [used]
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
    elseif sea_nr==7
    X=HadISST_M_season_reg_amjjason(:,:,1:end); 
    end     

% elseif type_nr==2 % for monthly reg [not used]
%     
%     if SST_dataset==1 % HadISST
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if  set_yr==1 &&  param_nr==1
%         load('HadISST_sst_c4.mat') % from HadISST_cosweight_clim.m 1979-
%         ntim=432;
%     
%     
%         elseif (set_yr==6 &&  param_nr==1)  || (iso_nr==16 && param_nr==1) % 1950-2014
%             
%         yr_s_nr=1;
%         
%             if yr_s_nr==1
%             load('C:\PHD\matlab_lib\Data\HadISST_sst_1950_2014_clim.mat')
%             ntim=780;
%     
%             elseif yr_s_nr==2
%             % 1871-2014
%             load('C:\PHD\matlab_storage_of_output_files\clim_SST_1871_2014.mat');
%             ntim=1728;
%             HadISST_sst_c4=era_z500_c4;
%         
%             end
%             
%         
%         elseif (set_yr==1 || set_yr==8  || set_yr==9) &&  param_nr==2
%         load('C:\PHD\matlab_storage_of_output_files\HadISST_clim_sic_1979_2014.mat');
%         ntim=432;    
%             
% 
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     elseif SST_dataset==2 % ERSST
%       
%             load('C:\PHD\matlab_storage_of_output_files\clim_SST_ERSST_1950_2014.mat');
%             ntim=780;
%             HadISST_sst_c4=era_z500_c4;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     end
%     
%     
%     if SST_dataset==1
%         nlat=180;
%         nlon=360;
%         
%     elseif SST_dataset==2    
%         nlat=89;
%         nlon=180;
%         
%     end
%     
%     HadISST_sst_c5 = reshape(HadISST_sst_c4,ntim, nlat, nlon); 
%     HadISST_sst_c6=permute(HadISST_sst_c5,[3 2 1]);
%   %  X=HadISST_sst_c6(:,:,(yr_1:yr_2));
%   
%         if  sea_nr==1
%             
%             if set_yr==1 || set_yr==10
%                 
%              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% %              if yr_e==2014   
% %              X=HadISST_sst_c6(:,:,(1:end));
% %              elseif yr_e==2013
% %              X=HadISST_sst_c6(:,:,(1:end-12));    
% %              end
%              
%              
%              
%              if yr_s==1871 || yr_s==1950
%              X=HadISST_sst_c6(:,:,(1:end));
%              elseif yr_s==1979 && param_nr==1
%              X=HadISST_sst_c6(:,:,(349:end));
%              elseif yr_s==1979 && param_nr==2
%              X=HadISST_sst_c6(:,:,(1:end)); 
%              
%              
%              end
%              
%              
%              %%%%%%%%%%%%%%%
%              
%             elseif set_yr==8
%                 if yr_e==2014
%                 X=HadISST_sst_c6(:,:,(253:end));   % 2000-2014 IPO-
%                 elseif yr_e==2013
%                 X=HadISST_sst_c6(:,:,(253:end-12));   % 2000-2013 IPO-
%                 end
%             elseif set_yr==9
%              X=HadISST_sst_c6(:,:,(1:252));   % 2000-2014 IPO-
%             
%             end
%         elseif  sea_nr>=2
% %             999
%             X=HadISST_sst_c6(:,:,(mm_in_c));
%         end
        
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
toc

%% Wilks test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
fprintf('[\bWilks test ]\b')
contour_alt_nr=2; % (1) p regular local, (2) with Wilks 2006, 2016 test
if contour_alt_nr==2
 

p_in= squeeze(p_all(:,:,1));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alt 2
% from Github 
% find the locations where the null hypothesis 
% could be rejected locally
% alp_c=.05;
alp_c=.1;
rej_null = find(p_in<alp_c); % changes from grid to vector

no_rej_null = length(rej_null); 
%sort from smallest to largest the p values that
%would suggest the null hypothesis could be rejected
%locally
[rej_null_sort, rej_null_ind] = sort(p_in(rej_null)); 
[rn_x,rn_y] = ind2sub([360,180],rej_null(rej_null_ind));
p_local_rej = ones(360,180);
p_fdr_rej = ones(360,180);
no_grid_points = 360*180;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Figure
%%%%%%%%%%%%%%%%%%%%%%%
tic
fprintf('[\bFigure ]\b')
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

% [c_max, c_min, c_limit]=Ma_max_min(s_all_c1);
% 
% c_max_c=c_max;
% 
% c_limit_c=c_limit;

    if lock_scalebar==1
         
         
        if  iso_nr==7 || iso_nr==15 % PCs
         % c_limit=0.5;
%                 if param_nr==1
%                 c_limit=1.0;
%                 elseif param_nr==2
%                 c_limit=0.1;
%                 end
                
                
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

fig('units','inches','width',14,'height',9,'font','Helvetica','fontsize',16,'border','on');
%%%%%%%%%%%
if param_nr==1
   axesm( proj,'MapLatLimit',[lat1+adj_nr lat2],'MapLonLimit',[lon1 lon2],'Grid','on','ParallelLabel','on','Frame','on',... %
       'MeridianLabel','on','FontWeight','bold','FontSize',14,...
       'mlabellocation',[0:60:360]);
   
elseif param_nr==2 % stereo
    if proj_nr==1
      axesm( proj,'MapLatLimit',[lat1 lat2],'Grid','on','ParallelLabel','on','Frame','on',... %
       'MeridianLabel','on','FontWeight','bold','FontSize',18,...
       'mlabellocation',[0:30:180,0:-30:-180]); 
   
    elseif proj_nr==2 % merc
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

%The position arguments are [xposition yposition width height].

        pos=get(br,'Position');
                pos(1)=pos(1)+x_move_colorbar; 
                pos(2)=pos(2)+y_move_colorbar;
                pos(3)=pos(3)+adj_width;  % widthcolorbar
                pos(4)=pos(4)+adj_height;  % height colorbar
        set(br,'Position',pos)
        br.TickLabels = strrep(br.TickLabels, '-', '–');% changes to minus symbol

        if color_alt==1
            colormap(brewermap(256,'*RdBu'));
        elseif color_alt==2
            colormap(flipud(cbrewer2('div','Spectral',10)));
        elseif color_alt==3
            colormap(flipud(cbrewer2('div','PiYG',20)));              
        elseif color_alt==4
            colormap(flipud(cbrewer2('div','RdYlGn',20)));          
        elseif color_alt==5    
             colormap(cbrewer2('div','BrBG',20));     
        elseif color_alt==6    
            colormap(flipud(cbrewer2('div','PuOr',20)));
        elseif color_alt==7    
            colormap(flipud(cbrewer2('div','RdGy',20)));      
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% significance contour
%     h2= contourm( double(HadISST_lat),double(HadISST_lon),squeeze(p_all(:,:,1))', [0.05],'-','ShowText','off',...
%     'Linecolor',[0.6 0.8 0.8],'LineWidth', line_w_nr);
    % Wilks test
    h1= contourm( double(HadISST_lat),double(HadISST_lon),p_fdr_rej', [0.05],'-','ShowText','off',...
    'Linecolor',co_pa,'LineWidth', line_w_nr);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
bedmap2('gl','LineWidth', 1.5, 'color',color_code);
bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end
%%%%%%%%%%%%%%%%%%%%


% label_size=16;
label_size=20;

% if cor_nr==4
%     a1=axestext_c(0.04, +0.05, ['(',label_1,')'] );
% else
%     a1=axestext_c(0.1, -0.04, ['(',label_1,' ', label_2,')'] );
% % text(-1.9,-1.8, [label_1,' ', label_2],'FontWeight','bold','FontSize',label_size,'rotation',90); 
% end
% set(a1,'FontWeight','bold','FontSize',label_size,'rotation',90);




label_vert=-300;
label_across=-64.9;
mlabel('MLabelParallel',label_across,'PLabelLocation',[-75 -60 -45 -30 -15  0 15 ],'PLabelMeridian',label_vert) 
% 'MLabelLocation',[ -120 -60  0 60 120 180],




    
    if  iso_nr==7
                 
         label_1='°C s.d.^-^1'; % per nomralized unit standard deviation
         
         
         if param_nr==1

         elseif param_nr==2
          label_1='SIC';   
     
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
        
        set(t4,'Position',pos,'FontSize',16, 'FontWeight', 'bold');    


        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         letter_size=28;
%         
%         if param_nr==2
%          letter_size=36;    
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',letter_size ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
     

     
     
     if param_nr==2 && ( iso_nr==7 ||  iso_nr==15)
         
           a1=axestext_c(x_lab,y_lab, [' '] );
         
     else
         
         a1=axestext_c(x_lab,y_lab, ['(',label_1,')'] );
         
     end
     
     

     
% text(-1.9,-1.8, [label_1,' ', label_2],'FontWeight','bold','FontSize',label_size,'rotation',90); 
label_size2=18;
set(a1,'FontWeight','bold','FontSize',label_size2,'rotation',90);
     
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
    
filename=['HadISST_regression_',name,'_',pc_str,'_',season,'_',type_str,'_',num2str(yr_s),'-',num2str(yr_e)];


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


savefilename_c=strcat(filedir,'figures\',filename);

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
cd('G:\My Drive\ISO_CFA\matlab')
  export_fig('-png','-nocrop','-painters', '-depsc','-opengl', quality_level_str, savefilename_c); % PNG-nocrop'
cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab')
toc
%% for saving PSA pattern surfaces (step 1)-= one more step=- 
  
% save regression surface  
% folder_c='C:\PHD\matlab_storage_of_output_files\';
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
 
 