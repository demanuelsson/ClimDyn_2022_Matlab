%%%%%%%%%%%%%%%%%%%%%%%%% 
% Regression ERA-Interim
%    
% Regression_spatial_ERA_I_z500_SAT_DE21.m
% Daniel Emanuelsson 2021
% Matlab 2020b
% Github version 1  [Uploaded]
% Github version 2  23-08-2020
%%%%%%%%%%
% subfunctions
% * fig.m fileexchange
% * export_fig.m fileexchange
% Wilks test
% saves surfaces, that is ploted further down
%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
clear
close all

     for    PC_nr=2:3 % (2) PC1 (SAM), (3) PC2 (PSA1)     (4) PC3 PSA2, (5) PC4


% site='RICE';

downscale_nr=1; % 0/1
seasonaly_nr=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cor_nr=4;  % (1) just trend, (4) PCs from EOF z500, (6) PCs from EOF 2mT
% change dataset too
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type_nr=1; %(0) if seasonaly(1) annual use for PCs (2) monthly, for trends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
era_name_nr=1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1. Z500 **** (Fig 2)
% 3. v850 **** (Fig. 11c) trend
% 4. u850
% 5. 2mT SAT ***** (Fig. 3)

if era_name_nr==1
    name='z500';
elseif era_name_nr==3 % uesd in Fig. 11c
    name='v850';
elseif era_name_nr==4    
    name='u850';    
elseif era_name_nr==5
    name='2mT'; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        yr_nr_s=1;%%%%%%%%%%%
        
        if yr_nr_s==1
            yr_s=1979;
            yr_e=2011;
        elseif yr_nr_s==2 % uesd in Fig. 11c
            yr_s=2000;
          %  yr_e=2013; 
            yr_e=2011; 
        end
        
    

        sea_nr=7; %%%%%% (1) Annual (4) JJA (5) SON (6) MAMJJASON (7) AMJJASON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
          elseif sea_nr==6
            season='MAMJJASON';
          elseif sea_nr==7
            season='AMJJASON';    
          end


%%%%%%%%%%Figure format
 figure_format=2; % (1) - EPS, (2) - PNG 
%%%%%%%%%%
proj='stereo';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in ERA-interim data
% addpath C:\PHD\ERA_interim
  addpath C:\Users\benman\matlab_lib\Data\ERA_interim

name_c='ERA_int_monthly_z500_2.nc';
%   ncdisp(name_c)

era_time=ncread(name_c,'time');
era_long=ncread(name_c,'longitude');
era_lat=ncread(name_c,'latitude');
% era_sst=ncread(name_c,'sst');
clear name_c

% era_date=((double(era_time)./24)./365.2400)+1900;

dayssince111=era_time/24;
datevalue=dayssince111+datenum(1900,1,1);
 date_vec=datevec(double(datevalue)); 
 yyyy=date_vec(:,1); 
 mm=date_vec(:,2);
 era_year_num=yyyy+(mm/12-1/12); 
 
%%%%%%%%%%%%%%%%%%%%%%%
 

if strcmp(name,'2mT')==1
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of daily means
    % surface 2mT
    
    name_c='ERA_int_monthly_2m_T.nc';
    %   ncdisp(name_c)
    era_T=ncread(name_c,'t2m'); % 2m temp
    era_T=era_T- 273.15;
    letter='a';
    lat1=-90;
    lat2=-30; 
    label_1='°C s.d.^-^1'; 
 elseif strcmp(name,'z500')==1 %&& strcmp(season,'annual')==1
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of dail means
    % geopotentail 500 hPa
    
    name_c='ERA_int_monthly_z500_2.nc'; 
    era_z500=ncread(name_c,'z'); % GPH z500 %%%%%% file _2 with 2014 and 2015 too ......
    era_z500=era_z500/9.80665;
    
    if PC_nr==2
    letter='a';
    elseif PC_nr==3
    letter='b';
    end
    
    label_1='m';
  
        if strcmp( proj,'stereo')==1
            lat1=-90;
            lat2=0;
        elseif strcmp( proj,'mercator')==1
            lat1=-90;
            lat2=40;
            lon1=-260;
            lon2=-30;
        end

  
 elseif strcmp(name,'v850')==1
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of dail means
    % V wind, pressure level 850 hPa 
    
    name_c='ERA_int_monthly_v850.nc'; 
%     era_v850=ncread(name_c,'v') ; % v850 'm s**-1'
    letter='c';
    lat1=-90; 
    lat2=-50;  
    if type_nr==1
    label_1='m s^-^1';
    end
    
  elseif strcmp(name,'u850')==1   
    
    name_c='ERA_int_monthly_u850.nc'; 
%     era_u850=ncread(name_c,'u'); % v850 'm s**-1'
    letter='c';
    lat1=-90; 
    lat2=-45;  
    if type_nr==1
    label_1='m s^-^1';
    end
 end


era_count=find(era_year_num==yr_e)+11;  % data in months
%era_start=find(era_year_num==1980-1);

%if strcmp(season,'annual')==1 % annual 
era_start=find(era_year_num==yr_s);


mm_in=mm(era_start:era_count);  % seasonal index

 if sea_nr==1
    mm_in_c=[era_start:era_count]';
 elseif sea_nr==2  % DJF
     mm_in_c=[find(mm_in==12);  find(mm_in==1); find(mm_in==2)];
     mm_in_c=sort( mm_in_c);
   %  mm_in_c=mm_in_c(3:end-1); % dont use JF for first year as there is no D
      mm_in_c=mm_in_c(1:end); % dont use JF for first year as there is no D    
elseif sea_nr==3  % MAM
     mm_in_c=[find(mm_in==3);  find(mm_in==4); find(mm_in==5)];
     mm_in_c=sort( mm_in_c);
     mm_in_c=mm_in_c(1:end);
     
elseif sea_nr==4  % JJA
     mm_in_c=[find(mm_in==6);  find(mm_in==7); find(mm_in==8)];
     mm_in_c=sort(mm_in_c);
     mm_in_c=mm_in_c(1:end);
     
 elseif sea_nr==5  % SON
     mm_in_c=[find(mm_in==9);  find(mm_in==10); find(mm_in==11)];
     mm_in_c=sort( mm_in_c);
     mm_in_c=mm_in_c(1:end);
 elseif sea_nr==6  % MAMJJASON    
     mm_in_c=[find(mm_in==3); find(mm_in==4); find(mm_in==5); find(mm_in==6); find(mm_in==7); find(mm_in==8); find(mm_in==9);  find(mm_in==10); find(mm_in==11)];
     mm_in_c=sort( mm_in_c);
     mm_in_c=mm_in_c(1:end); 
 elseif sea_nr==7  % AMJJASON    
     mm_in_c=[find(mm_in==4); find(mm_in==5); find(mm_in==6); find(mm_in==7); find(mm_in==8); find(mm_in==9);  find(mm_in==10); find(mm_in==11)];
     mm_in_c=sort( mm_in_c);
     mm_in_c=mm_in_c(1:end);      
     
 end
%%%%%%%%%%%%%%%%%%%   
toc
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %ERA_M_annual_reg
 tic
 if cor_nr==1  
         iso=[]; 
     
        if      type_nr==1 % annual
 
            wr_yr=[1979:2014]';
            yr_1=find(wr_yr==yr_s);       
            yr_2=find(wr_yr==yr_e);
            
                if sea_nr==1 % annual
                y=wr_yr(yr_1:yr_2);
                elseif sea_nr>=2  % Seasonal (1980-
                y=wr_yr(yr_1+1:yr_2);
                end
                
        elseif   type_nr==2 % monthly
                
                y=[1:era_count]';            
%                 if sea_nr==1
                y=y(mm_in_c);
%                 elseif sea_nr>=2              
        end
                
                
            label_2='yr^-^1';
                if era_name_nr==1
                label_1='m decade^-^1';
                elseif era_name_nr==3 || era_name_nr==4
                label_1='m s^-^1 decade^-^1';
                end
        
                letter='c';
          
  elseif cor_nr==4        
         

      if sea_nr==1 && type_nr==2
             
      elseif sea_nr==1 && type_nr==1
          
      elseif sea_nr==7    
       %  load('C:\Users\benman\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2011_MAMJJASON_c52_annual_mean_varimax_c2.mat');
         load('C:\Users\Machine\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2011_AMJJASON_c52_annual_mean_varimax_c2.mat');  %  % AMJJASON 
      end
               
                  label_1='m s.d.^-^1';
    
    
        if type_nr==1   %%%%%%%%% Annual 
            time_c=MA_PCs_save(:,1);
        elseif seasonaly_nr==1 && sea_nr>=2   
            time_c=MA_PCs_save_seasonaly(:,1);
        elseif type_nr==2 %%%%%%%%%%% Monthly
            time_c=MA_PCs_save_monthly(:,1);
        end
        
        
        
            yr_1=find( time_c==yr_s);       
            yr_2=find( time_c==yr_e);
    
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %     PC_nr=2; % (2) SAM (3) PSA1 (4) PSA2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (1) is just time
%         for    PC_nr=2:4  
            
            
            label_2='';
            iso='z500';
                                      paper_nr=2; % (1) iso (2) EOF 
    
        if PC_nr==2  % SAM
            
                fact_nr=1;  %2011
                 if sea_nr==1 || sea_nr==6 || sea_nr==7
                                     
                    if era_name_nr==1

                            variance_exp=num2str(round(reof_out2.rexpv(1)*10)/10); % varimax
                            letter='a';
                    elseif era_name_nr==3
                            letter='c';
                    elseif era_name_nr==5
                            %letter='d';
                            letter='g';
                    end  
                end
                
             
            elseif PC_nr==3 % PSA1

               fact_nr=1; %2011 
                if sea_nr==1 || sea_nr==6 || sea_nr==7

                        variance_exp=num2str(round(reof_out2.rexpv(2)*10)/10); % varimax
                        
                        if era_name_nr==1
                            letter='b';
                         elseif era_name_nr==5
                          %  letter='e';
                            letter='h';
                         end
                end
                
                
    
              
            elseif PC_nr==4 % PSA2
               %  fact_nr=-1;  %2014  To get sign as in Kidson 1988 paper
                 fact_nr=1;% 2011
                 if sea_nr==1
                     %variance_exp=10;
                    % variance_exp=8.2;
                    % variance_exp=9; % varimax
                    variance_exp=9.9; % varimax
                        if era_name_nr==1
                            letter='c';
                        elseif era_name_nr==5
                            % letter='f';
                            letter='i';

                        end 
                 end
        end
            
            
            
         if type_nr==2 %  monthly % works for all seasons
%                 if sea_nr==1
               % y=MA_PCs_save_monthly((yr_1:yr_2),PC_nr)*fact;%%%%%%%%%%%%%%%%%%%%%%%%%%%
               if sea_nr==1
              y=MA_PCs_save_monthly((mm_in_c),PC_nr)*fact_nr;
              
               elseif seasonaly_nr==1
               y=MA_PCs_save_seasonaly(:,PC_nr)*fact_nr;
               else
            %  y=MA_PCs_save_monthly(:,PC_nr)*fact_nr;  
            %    y=MA_PCs_save_monthly(:,PC_nr)*fact_nr;
              y=MA_PCs_save_monthly(1:length(mm_in_c),PC_nr)*fact_nr;     % seasonal
               end
         elseif type_nr==1
             
             y=MA_PCs_save((yr_1:yr_2),PC_nr)*fact_nr;
               
          end
  
 end
 toc
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
% regression of deseasonal 
%%%%%%%%%%%%%%%%%
tic

% method_nr=2;  % Use (2) regstats. It works fine

  
    
if  era_name_nr==1  %%%%% Monthly anomalies data
    if downscale_nr==0
        
% cosweight and monthly clim removed global z500 not downscaled (run start of eof code)
     load('C:\Users\Machine\matlab_lib\Data\era_z500_c4.mat');
     ntim=444;
     nlat=241;
     nlon=480;
     
    elseif downscale_nr==1
    load('C:\Users\Machine\matlab_storage_of_output_files\clim_z500_ERA-Interim_downscale_0.75_1979_2012.mat'); % downscaled
    
    %eof_code_run_c4_annual_seasonal.m
    
    ntim=408;
    nlat=90;
    nlon=360; 
    
    end

%     
%     ntim=408;
%     nlat=47;
%     nlon=240;
    

    
    era_z500_c6_reg_c2 = reshape(era_z500_c4,ntim, nlat,nlon); 
    era_z500_c7_reg_c2=permute(era_z500_c6_reg_c2,[3 2 1]);
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    if seasonaly_nr==1 && sea_nr>2
        
    era_z500_c4=era_z500_c4(mm_in_c,:);
      
      seasonal_str='seasonaly';
       
        [ac1 ac2 ]=size(era_z500_c4);
        
        Ma_save_season=NaN(ac1/3,ac2);
        Ma_save_season_t=NaN(ac1/3,1);
        kc=1;
        
            for ib=1:3:ac1
        
           dummy_c1 =mean(era_z500_c4((ib:ib+2),:));
           
%            dummy_t=mean(era_year_num_c3(ib:ib+2));
           
           Ma_save_season(kc,:)=dummy_c1;
%            Ma_save_season_t(kc,:)=dummy_t;
           
           kc=kc+1;
           
           
            end
        
      era_z500_c4=Ma_save_season; 
           
          ntim_c=33;
          era_z500_c6_reg_c2 = reshape(era_z500_c4,ntim_c, nlat,nlon); 
          era_z500_c7_reg_c2=permute(era_z500_c6_reg_c2,[3 2 1]); 
          X=era_z500_c7_reg_c2;
          
    elseif type_nr==1 %%%%%%%%%%%%%%%%% annual
        
%         era_z500_c4=era_z500_c4(mm_in_c,:);
      
      seasonal_str='annualy';
       
        [ac1 ac2 ]=size(era_z500_c4);
        
        Ma_save_season=NaN(ac1/12,ac2);
%         Ma_save_season_t=NaN(ac1/12,1);
        kc=1;
        
            for ib=1:12:ac1
        
           dummy_c1 =mean(era_z500_c4((ib:ib+11),:));
           
%            dummy_t=mean(era_year_num_c3(ib:ib+2));
           
           Ma_save_season(kc,:)=dummy_c1;
%            Ma_save_season_t(kc,:)=dummy_t;
           
           kc=kc+1;
           
           
            end
        
      era_z500_c4=Ma_save_season((yr_1:yr_2),:); 
           
          %ntim_c=34;
          ntim_c=size([yr_1:yr_2]');
          era_z500_c6_reg_c2 = reshape(era_z500_c4,ntim_c(1), nlat,nlon); 
          era_z500_c7_reg_c2=permute(era_z500_c6_reg_c2,[3 2 1]); 
          X=era_z500_c7_reg_c2;      
      
      
    elseif type_nr==2
        
        
    %X=era_z500_c7_reg_c2(:,:,(yr_1:yr_2));%%%%%%%%%%%%%%%%%%%%%%%%
     X=era_z500_c7_reg_c2(:,:,(mm_in_c));
     % time_check=era_year_num(mm_in_c);
    
        
        
     
    end
     
  %%%%%%%%%%      v850
     elseif type_nr==2  && era_name_nr==3%%%%% v850 Monthly data  
         
     filedir ='C:\Users\Machine\matlab_storage_of_output_files\';
     if downscale_nr==0
     
%       load([filedir,'ERA_interim_v850_clim.mat']); % cosweight and monthly clim removed global z500 not downscaled (run start of eof code)
     X=ERA_interim_v850_clim(:,:,(mm_in_c));  

     elseif downscale_nr==1

      load([filedir,'clim_v850_ERA-Interim_downscale_0.75_1979_2012.mat']); 
      
    %eof_code_run_c4_annual_seasonal.m
    
%         ntim=444; 
%         nlat=241;
%         nlon=480;

    ntim=408;
    nlat=90;
    nlon=360;


%     
         era_z500_c6_reg_c2 = reshape(era_z500_c4,ntim, nlat,nlon); 
         era_z500_c7_reg_c2=permute(era_z500_c6_reg_c2,[3 2 1]);
    
   % X=era_z500_c7_reg_c2(:,:,(yr_1:yr_2));%%%%%%%%%%%%%%%%%%%%%%%%
    X=era_z500_c7_reg_c2(:,:,(mm_in_c));
     end 
    
     
        
        
  %%%%%%%%%%  
  
  elseif type_nr==2  && era_name_nr==4%%%%% u850 Monthly data  

       filedir ='C:\Users\Machine\matlab_storage_of_output_files\';
       
%          ntim=432; % already done for this one
%          nlat=241;
%          nlon=480;
       
       
       if downscale_nr==0
       load([filedir,'clim_u850_ERA-Interim_1979_2014.mat']);
       
       
       elseif  downscale_nr==1
        load([filedir,'clim_u850_ERA-Interim_downscale_0.75_1979_2012.mat']);
                % downscaled      
        ntim=408;
        nlat=90;
        nlon=360;
       end
       
         era_z500_c6_reg_c2 = reshape(era_z500_c4,ntim, nlat,nlon); 
         era_z500_c7_reg_c2=permute(era_z500_c6_reg_c2,[3 2 1]);
         
       X=era_z500_c7_reg_c2(:,:,(mm_in_c)); 
  
   
elseif era_name_nr==5%%%%% 2mT Monthly data    

            filedir ='C:\Users\Machine\matlab_storage_of_output_files\';
            if downscale_nr==0
            load([filedir,'era_2mT_c8.mat']); % cosweight and monthly clim removed global z500 not downscaled (run start of eof code)     
            
        ntim=444;
        nlat=241;
        nlon=480;
            era_z500_c6_reg_c2= reshape(era_z500_c8,ntim, nlat,nlon); 
            era_z500_c7_reg_c2= permute(era_z500_c6_reg_c2,[3 2 1]);
        
        
            elseif downscale_nr==1
        load([filedir,'clim_2mT_ERA-Interim_downscale_0.75_1979_2012.mat']);
        ntim=408;
        nlat=60;
        nlon=360;
        
         era_z500_c6_reg_c2 = reshape(era_z500_c4,ntim, nlat,nlon); 
         era_z500_c7_reg_c2=permute(era_z500_c6_reg_c2,[3 2 1]);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
         if seasonaly_nr==1 && sea_nr>2
             
            era_z500_c4=era_z500_c4(mm_in_c,:);
      
             
          seasonal_str='seasonaly';
       
        [ac1 ac2 ]=size(era_z500_c4);
        
        Ma_save_season=NaN(ac1/3,ac2);
        Ma_save_season_t=NaN(ac1/3,1);
        kc=1;
        
            for ib=1:3:ac1
        
           dummy_c1 =mean(era_z500_c4((ib:ib+2),:));
           
%            dummy_t=mean(era_year_num_c3(ib:ib+2));
           
           Ma_save_season(kc,:)=dummy_c1;
%            Ma_save_season_t(kc,:)=dummy_t;
           
           kc=kc+1;
                   
            end
        
          era_z500_c4=Ma_save_season; 
           
          ntim_c=33;
          era_z500_c6_reg_c2 = reshape(era_z500_c4,ntim_c, nlat,nlon); 
          era_z500_c7_reg_c2=permute(era_z500_c6_reg_c2,[3 2 1]);
      
      
      
          X=era_z500_c7_reg_c2;    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
         elseif type_nr==1 %%%%%%%%%%%%%%%%% annual
        
%         era_z500_c4=era_z500_c4(mm_in_c,:);
      
      seasonal_str='annualy';
       
        [ac1 ac2 ]=size(era_z500_c4);
        
        Ma_save_season=NaN(ac1/12,ac2);
%         Ma_save_season_t=NaN(ac1/12,1);
        kc=1;
        
            for ib=1:12:ac1
        
           dummy_c1 =mean(era_z500_c4((ib:ib+11),:));
           
%            dummy_t=mean(era_year_num_c3(ib:ib+2));
           
           Ma_save_season(kc,:)=dummy_c1;
%            Ma_save_season_t(kc,:)=dummy_t;
           
           kc=kc+1; 
       
          
            end
            
          era_z500_c4=Ma_save_season((yr_1:yr_2),:); 
           
          %ntim_c=34;
          ntim_c=size([yr_1:yr_2]');
          era_z500_c6_reg_c2 = reshape(era_z500_c4,ntim_c(1), nlat,nlon); 
          era_z500_c7_reg_c2=permute(era_z500_c6_reg_c2,[3 2 1]); 
          X=era_z500_c7_reg_c2; 
            
            
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
         elseif type_nr==2
         X=era_z500_c7_reg_c2(:,:,(mm_in_c)); 
         end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            
            end
        
 
    
    %X=era_z500_c7_reg_c2(:,:,(yr_1:yr_2));%%%%%%%%%%%%%%%%%%%%%%%%
%             X=era_z500_c7_reg_c2(:,:,(mm_in_c));      
     
    
end
    
[A B C]=size(X);
% b=zeros(2, A,B);bINT=zeros(2,2,A,B);r=zeros(324,A,B);rINT=zeros(324,2,A,B);
b=zeros(2, A,B);
bINT=zeros(2,2,A,B);
r=zeros(C,A,B); % residual
rINT=zeros(C,2,A,B);
STATS=zeros(1,4,A,B);  % 3 p-value

 p_all=zeros(A,B,1);  s_all=zeros(A,B,1);  % just one layer needed?

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Regression regstats function
tic
 
% for i=1:A
%     for j=1:B
        
parfor i=1:A
    for j=1:B        
                        
            s=regstats(squeeze(X(i,j,:)),y,'linear','all');
            s.tstat;
            s_all(i,j,:)=s.tstat.beta(2); % beta- regresion coeficient 
            p_iso_rec=s.tstat.pval;
            p_all(i,j,:)= p_iso_rec(2); % The second element is the p-value for a significantly non-zero slope.
                       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
         [a, MSGID] = lastwarn(); % Comment/uncomment if needed
         warning('off', MSGID)     
            
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wilks test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contour_alt_nr=2; % (1) p regular local, (2) with wilks 2006, 2016 test

if contour_alt_nr==2
p_in=p_all;

% from github 
%find the locations where the null hypothesis 
%could be rejected locally
%alp_c=.05;
alp_c=.1;
rej_null = find(p_in<alp_c); % changes from grid to vector

no_rej_null = length(rej_null); 
%sort from smallest to largest the p values that
%would suggest the null hypothesis could be rejected
%locally
[rej_null_sort, rej_null_ind] = sort(p_in(rej_null)); 
[rn_x,rn_y] = ind2sub([360,nlat],rej_null(rej_null_ind));
p_local_rej = ones(360,nlat);
p_fdr_rej = ones(360,nlat);
no_grid_points = 360*nlat;
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
   fdr_thres = alp_c*((no_rej_null+1)-rej_null_ind(i))/no_grid_points; % remove *2
    
    
    
    %find local p values that are less than fdr threshold
    if(pi_c<fdr_thres)
        p_fdr_rej(rn_x(i),rn_y(i)) = pi_c;
        n_fdr = n_fdr+1;
    end 
end


 
end

 toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% turn off warning
% w = warning('query','last')
% id = w.identifier
% warning('off',id)
toc
%  s downscaled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
corr_label=1; %(1/0) %%%%%%%%%%% show colorbar and label
lock_scalebar=1; % (1/0)
%%%%


if type_nr==2 && cor_nr==1%%% TREND only    Monthly
    
    s_all_c=s_all;
%     s_all_c(p_all>=0.1)=NaN; 
    s_all_c=s_all_c*10*12; % per decade. %%%%%%  Agree with Bromwich et al. 2012
    

elseif   type_nr==2 || type_nr==1
    
   s_all_c=s_all;  % unit per s.d.   
end
%%%%%%%%%%%%%%%%%%%%%%%

[c_max, c_min, c_limit]=Ma_max_min(s_all_c);
c_max_c=c_max;
c_limit_c=c_limit;

    if lock_scalebar==1
             
         if    era_name_nr==5 && (cor_nr==4 || cor_nr==6)
%                  c_limit= 1.7879;% 2mT EOF1 2014 max                  
            %       c_limit= 1.51;% 2mT EOF1 2014
                   c_limit= 1.01;% 2mT EOF1 2012
                                
         elseif    era_name_nr==1 && cor_nr==4 && (sea_nr==1 || sea_nr==7)
                        % c_limit= 47.5260;
                         
                       % c_limit= 30; 
                        c_limit= 22.01; 
         elseif    era_name_nr==3 && cor_nr==1 && yr_nr_s==2         
                         
                        c_limit= 1.00; % 1.6590
         end
      
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig

fig('units','inches','width',10,'height',9,'font','Helvetica','fontsize',16,'border','on');

   axesm( proj,'MapLatLimit',[lat1 lat2],'Grid','on','ParallelLabel','on','Frame','on',... %
       'MeridianLabel','on','FontWeight','bold','FontSize',18,...  % font size of degree labels
       'mlabellocation',[0:30:180,0:-30:-180]); 

  
set(gca,'box','off','XColor',[1,1,1],'YColor',[1,1,1]);  % turns off border box and axes  
hSurf=surfm(double(era_lat),double(era_long),squeeze(s_all_c(:,:,1))');
hold on
 
        if corr_label==1
            br=colorbar('FontWeight', 'bold', 'FontSize',18);   
            pos=get(br,'Position');  
        
            if cor_nr==4 && era_name_nr==1
            
            pos(1)=pos(1)+0.04;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar 
            pos(2)=pos(2)+0.05;      
            pos(3)=pos(3)+0.014;  % widthcolorbar
            pos(4)=pos(4)-0.1;  % height colorbar 
            
            
            elseif cor_nr==4 && era_name_nr==5 || cor_nr==6 && era_name_nr==5
                
            
            pos(1)=pos(1)+0.04;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar 
            pos(2)=pos(2)+0.05;      
            pos(3)=pos(3)+ 0.014;  % widthcolorbar
            pos(4)=pos(4)-0.1;  % height colorbar    
                
            

            elseif cor_nr==1 && era_name_nr==3
            
            pos(1)=pos(1)+0.02;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar 
            pos(2)=pos(2)+0.09;  
            pos(3)=pos(3)+ 0.000;  % widthcolorbar
            pos(4)=pos(4)-0.195;  % height colorbar  
                    
            else
            pos(1)=pos(1)+0.06;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar 
            pos(2)=pos(2)+0.15;      
            pos(3)=pos(3)-0.0045;  % widthcolorbar
            pos(4)=pos(4)-0.3;  % height colorbar    
        
            end
        
        set(br,'Position',pos)     
        
        end
        
    
      colormap(b2r(-c_limit,c_limit));
      
      
      if era_name_nr==1 % z500 %%%%%%%%
       % color_alt=4;
        color_alt=6;
      elseif era_name_nr==5
        color_alt=1;
      else
        color_alt=1;    
      end
      
      
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
        end
      
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % p=0.05 significance level contour      
      
    if cor_nr==4 || cor_nr==5   || cor_nr==6 || cor_nr==3 || cor_nr==1    
        co_pa=[.1 .1 .1];
        
    else
        co_pa=[.9 .9 .9];
    end


    if contour_alt_nr==1
    h1= contourm( double(era_lat),double(era_long),squeeze(p_all(:,:,1))', [0.05],'--','ShowText','off',...
    'Linecolor',co_pa,'LineWidth', 2);

    elseif contour_alt_nr==2 %wilks

    h1= contourm( double(era_lat),double(era_long),p_fdr_rej', [0.05],'--','ShowText','off',...
    'Linecolor',co_pa,'LineWidth', 2);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coast_nr=1;%%%%(1/0) on/off %%%%% Coast lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if coast_nr==1

    color_nr=1;
    if color_nr==1
       color_code=[.1 .1 .1];   
    elseif color_nr==2
       color_code=[.6 .6 .6];  
    end

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
plot3m(lat_cr,lon_cr,'-k','LineWidth', 1.5,'color',color_code);
%plot3m(lat(in_c),long(in_c),'.k','MarkerSize', 1);

%addpath C:\PHD\matlab_mapping
addpath C:\Users\benman\matlab_mapping
bedmap2('gl','k','LineWidth', 1.5,'color',color_code);
bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% label above colorbar
label_size=22;
  
if cor_nr==1 && era_name_nr==3 %
    a1=axestext_c(1.101, +0.33, ['(',label_1,')'] );    
elseif cor_nr==4 || cor_nr==5 || cor_nr==1 || cor_nr==6
    a1=axestext_c(0.04, +0.028, ['(',label_1,')'] );    
else
    a1=axestext_c(0.1, -0.04, ['(',label_1,' ', label_2,')'] );
end

set(a1,'FontWeight','bold','FontSize',label_size,'rotation',90);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variance explained by EOF
        if ((cor_nr==4 && era_name_nr==1) || (cor_nr==6 && era_name_nr==5)) && sea_nr==1  
            a2=axestext_c(0.95, +0.98, [num2str(variance_exp),'%'] );
            set(a2,'FontWeight','bold','FontSize',label_size+4);           
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lat long labels. Long  at this lat and placement of lat         
    if cor_nr==6 || cor_nr==4 && era_name_nr==5
        mlab_par=-28;
    elseif cor_nr==1 && era_name_nr==3    
        mlab_par=-64; 
    else
        mlab_par=-10;
    end
        
mlabel('MLabelParallel',mlab_par,'PLabelLocation',[-75 -60 -45 -30 -15  0],'PLabelMeridian',100) 


%%%%%%%%%%%%%%%%%%%%%%
%  Title

if cor_nr==1    
        if type_nr==2
            type_str='Monthly';
        else
            type_str='Annual';    
        end
        
        t4=title(['Trend ', name,' ', type_str,' ',num2str(yr_s),'-',num2str(yr_e)]);   
            
elseif cor_nr==3
        t4= title(['Regression ', iso,'/', name]);
   
elseif cor_nr==4 || cor_nr==5  % EOFs 
    
        if type_nr==2
            type_str='Monthly';
        else
            type_str='Annual';    
        end
    
    if PC_nr==2
        pc_str='SAM';
    elseif PC_nr==3
        pc_str='PSA1';
     elseif PC_nr==4
        pc_str='PSA2';        
    end 
    
    if sea_nr==1
        season='Annual';
    end
     
            if cor_nr==4 && sea_nr==1 &&  era_name_nr==1
                t4=title([ pc_str,' ','EOF',num2str(PC_nr-1)]);   
%             elseif cor_nr==4 && sea_nr==1 &&  era_name_nr==2
%                 t4=title([ pc_str,' ','EOF',num2str(PC_nr-1),' ',type_str]);   
             
             elseif cor_nr==4 && sea_nr==1 &&   era_name_nr==5
                 
%                     if PC_nr<=3
                        t4=title([ pc_str,' ','PC',num2str(PC_nr-1),iso,'Annual']);   
%                     elseif PC_nr==4
%                         t4=title([ pc_str,' ','PC',num2str(PC_nr-1),', SAT ',type_str,' (-1)']);       
%                     end
                        
             
            elseif cor_nr==4 && sea_nr>=1 &&  ( era_name_nr==1 || era_name_nr==5 )
%             t4=title([ pc_str,' ',name,' ',season]);  
            
%             elseif cor_nr==4 && sea_nr==1 &&  era_name_nr==2
%               t4=title([ pc_str,' ',season,' ', name]);  
                
            elseif cor_nr==5
             t4=title([ pc_str]);      
            end
            
 elseif cor_nr==6           
  
        if type_nr==2
        type_str='Monthly';
        else
        type_str='Annual';    
        end
        
            if PC_nr==2
                pc_str='PC1';
                %pc_str='EOF1';
                
            elseif PC_nr==3
              pc_str='PC2';
              %  pc_str='EOF2';
            elseif PC_nr==4
                pc_str='PC3';
                 %    pc_str='EOF3';
            elseif PC_nr==5     
                pc_str='PC4'; 
            end

                    t4=title(['EOF',num2str(PC_nr-1),', SAT ',type_str]);  
       
 
end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Letter
        letter_pos= [170 710 50 62];
        letter_size=38;
    

          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',letter_size ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
     
     %%%%%%%%%%%%%%%%%
     % colorbar fontsize 
     if cor_nr==4
        scalebar_fontsize=22;    
     else
        scalebar_fontsize=16;
     end
     
      set(br, 'FontSize',scalebar_fontsize, 'FontWeight', 'bold');
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Variance explained by EOF
        if cor_nr==4 && era_name_nr==1 
            a2=axestext_c(0.95, +0.98, [num2str(variance_exp),'%'] );
            set(a2,'FontWeight','bold','FontSize',label_size+4);
        elseif cor_nr==6 && era_name_nr==2 && PC_nr<=3    
            a2=axestext_c(0.9, +0.90, [num2str(variance_exp),'%'] );
            set(a2,'FontWeight','bold','FontSize',label_size+4);            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fig
        str_end='_c';
     
     
        if cor_nr==4 || cor_nr==6      
            filename=['ERA_interim_regression_',name,'_',pc_str,'_',season,...
                '_',type_str,'_',num2str(yr_s),'-',num2str(yr_e),str_end];
        
        elseif cor_nr==1 
            filename=['ERA_interim_regression_',name,'_',season,'_',type_str,'_',num2str(yr_s),'-',num2str(yr_e)];
        
        else
            filename=['ERA_interim_regression_',name,'_',iso,'_',season];
        end
        

filedir ='C:\Users\Machine\matlab_storage_of_output_files\figures\';

if corr_label==1
savefilename_c=strcat(filedir,filename,'_crop');
else
savefilename_c=strcat(filedir,filename);

end


% save as png 
orient landscape
cd('G:\My Drive\ISO_CFA\matlab')
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r240',savefilename_c);
export_fig('-pdf','-painters', '-depsc','-opengl', '-r240',savefilename_c); % PNG-nocrop'
cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_fields=1;

if save_fields==1

 filedir ='C:\Users\Machine\matlab_storage_of_output_files\';
 filename=strcat(filedir,filename);
savefilename_c=strcat(filename,'.mat');
save(savefilename_c,'p_all','s_all','p_fdr_rej','era_lat','era_long');

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Save EOF surfaces (use the next one)
% 
% filedir ='C:\PHD\matlab_storage_of_output_files\';
% 
%       filename_c=['ERA_interim_regression_surface_',name,'_',pc_str,'_',season,'_',type_str,'_',num2str(yr_s),'_',num2str(yr_e)];
%       savefilename_c=strcat(filedir,filename_c);
% 
%       %save
% 
%     save(savefilename_c,'s_all','p_all'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% For saving regression pattern surfaces -= one more step=- (use)
 %  s_all_c=s_all_c1;

 save_con=1;
 if  save_con==1
    p_level_rc=0.05;
    %p_level_rc=0.01;
if contour_alt_nr==1
 s_all_c(p_all>=p_level_rc)=NaN;
elseif contour_alt_nr==2
 s_all_c(p_fdr_rej>=p_level_rc)=NaN; % including wilks test
end

%%%%%%%%%%%%%
% [wa, wb]=size(s_all_c); 
% save regression surface  

folder_c='C:\Users\Machine\matlab_storage_of_output_files\';

%savefilename =[folder_c,filename,'_',num2str(p_level_rc),'_c2.mat'];
savefilename =[filename,'_',num2str(p_level_rc),'_c6.mat'];
%_c2 rotate

save(savefilename,'s_all_c'); 
 end


     end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %% Alt 2 save surface -= Loads saved surfaces from above =- make contours for overlay on corr plots, Z500, SAT corr
 tic
p_level_rc=0.05;
era_name_nr=1; % 1, 5
% season='annual';
% type_str='Monthly';
%season='Annual';
%season='MAMJJASON';
season='AMJJASON';
type_str='Annual';% one value per year

yr_s=1979;
yr_e=2011;
%  
  folder_c='C:\Users\Machine\matlab_storage_of_output_files\';

    if era_name_nr==5
%         load([folder_c,'ERA_interim_regression_2mT_PC1_',season,'_',type_str,'_',...
%             num2str(yr_s),'-',num2str(yr_e),'_c_',num2str(p_level_rc),'_c5.mat']);
        load([folder_c,'ERA_interim_regression_2mT_SAM_',season,'_',type_str,'_',...
            num2str(yr_s),'-',num2str(yr_e),'_c_',num2str(p_level_rc),'_c5.mat']);

    elseif era_name_nr==1
        load([folder_c,'ERA_interim_regression_z500_SAM_',season,'_',type_str,'_',...
            num2str(yr_s),'-',num2str(yr_e),'_c_',num2str(p_level_rc),'_c6.mat']);
    end      

    

s_all_c_pc1=s_all_c;  % this is the p-level surface as rs with lower p-levels are masked above

%%%%%%%%%%%%
 ind_pc1_nonnan=~isnan(s_all_c_pc1);
%%%%%%%%%%

s_all_c_pc1(find(s_all_c_pc1<0))=-1; % Negative values -1
s_all_c_pc1(find(s_all_c_pc1>0))=1;  % positive values 1


    if era_name_nr==5
%             load([folder_c,'ERA_interim_regression_2mT_PC2_',season,'_',type_str,'_',...
%                 num2str(yr_s),'-',num2str(yr_e),'_c_',num2str(p_level_rc),'_c5.mat']);
            load([folder_c,'ERA_interim_regression_2mT_PSA1_',season,'_',type_str,'_',...
                num2str(yr_s),'-',num2str(yr_e),'_c_',num2str(p_level_rc),'_c6.mat']);
            
    elseif era_name_nr==1
            load([folder_c,'ERA_interim_regression_z500_PSA1_',season,'_',type_str,'_',...
            num2str(yr_s),'-',num2str(yr_e),'_c_',num2str(p_level_rc),'_c6.mat']);
    end

s_all_c_pc3=s_all_c;

%%%%%%%%%%%%%%%%%
 ind_pc3_nonnan=~isnan(s_all_c_pc3);
%%%%%%%%%%%%%%%%%

s_all_c_pc3(find(s_all_c_pc3<0))=-1;  % Negative values -1
s_all_c_pc3(find(s_all_c_pc3>0))=1;  % positive values 1



% [~,ind_test_c]=ismember(s_all_c_psa1,s_all_c_psa2); % is not member of 

%%%%%%%%%%%%%%
 ind_psa_sum_nonactive=ind_pc1_nonnan+ind_pc3_nonnan;
%%%%%%%%%%%%%

 ind_psa_sum=s_all_c_pc1+s_all_c_pc3;
 ind_psa_sum_c=ind_psa_sum;
 
ind_nan= isnan(ind_psa_sum_c); 
ind_psa_sum_c(ind_nan)=9999; 

    if era_name_nr==5
        filename=['pc_2mT_ind_psa_sum_',num2str(p_level_rc),'_c5.mat'];
    elseif era_name_nr==1
        filename=['pc_z500_ind_psa_sum_',num2str(p_level_rc),'_',season,'_c5.mat'];
    end


 savefilename =[folder_c,filename]; 
% save(savefilename,'ind_psa_sum'); 
 save(savefilename,'ind_psa_sum_c','ind_psa_sum_nonactive'); 
 toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot from saved fields (using first section now, as it runs faster on the new PC)
% 
clear
close all
era_name_nr=1; % (1) Z500 (2) 2mT
PC_nr=1; % 1,2,3,11,12,13

cor_nr=4; % (4) Z500 PCs,(6) SAT PCs
yr_nr_s=1;
corr_label=1;


sea_nr=1;
if sea_nr==1
season='annual';
end


%type_str='monthly';
type_str='Annual';
proj='stereo';
lat1=-90;
yr_s=1979;
yr_e=2011;

if era_name_nr==1
   era_name_str='Z500';
   title_str=era_name_str;
   lat2=-20;
   label_1='m s.d.^-^1';
   if PC_nr==1
    PC_nr_str='SAM';
    letter='a';
   % variance_exp=23;
     variance_exp=40; % rotate
%    variance_exp=25; % varimax
%    variance_exp=30; % varimax
   elseif PC_nr==2
    PC_nr_str='PSA1';
    letter='b';
    variance_exp=11;
    %variance_exp=14;
    % variance_exp=12; % varimax
   %  variance_exp=20; % varimax
     
   elseif PC_nr==3    
    PC_nr_str='PSA2';
    letter='c';
    %variance_exp=10;
    variance_exp=8;
   %  variance_exp=9; % varimax
   %   variance_exp=10; % varimax
   end   
   
   
   
elseif era_name_nr==2
   era_name_str='2mT';
   lat2=-30;
   label_1='°C s.d.^-^1';
   title_str='SAT';
   if PC_nr==1
    PC_nr_str='PC1';
    letter='d';
    %variance_exp=20;
    %variance_exp=17;
    % variance_exp=22; % varimax
     variance_exp=30; % 
     
   elseif PC_nr==2
    PC_nr_str='PC2';
  %  letter='e';
    letter='f';
    %variance_exp=13;
    %variance_exp=11;
    % variance_exp=12; % varimax
     variance_exp=12;
   elseif PC_nr==3    
    PC_nr_str='PC3';
 %   letter='f';
    letter='e';
    %variance_exp=8;
   % variance_exp=6;
%      variance_exp=8; % varimax
     variance_exp=7; 
   elseif PC_nr==11
     PC_nr_str='SAM'; %? name
    letter='g';
      
   elseif PC_nr==12
      PC_nr_str='PSA1';
      letter='h';
   elseif PC_nr==13    
     PC_nr_str='PSA2';
     letter='i';
   end
   
   
   
   
end


folder_c='C:\Users\Machine\matlab_storage_of_output_files\';
%load([folder_c,'ERA_interim_regression_',era_name_str,'_',PC_nr_str,'_annual_Monthly_1979-2011_c.mat']);
load([folder_c,'ERA_interim_regression_',era_name_str,'_',PC_nr_str,'_Annual_Annual_1979-2011_c.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
corr_label=1; %(1/0) %%%%%%%%%%% show colorbar and label
lock_scalebar=1; % (1/0)
%%%%
type_nr=2;

if type_nr==2 && cor_nr==1%%% TREND only    Monthly
    
    s_all_c=s_all;
%     s_all_c(p_all>=0.1)=NaN; 
    s_all_c=s_all_c*10*12; % per decade. %%%%%%  Agree with Bromwich et al. 2012
    

elseif   type_nr==2 || type_nr==1
    
   s_all_c=s_all;  % unit per s.d.   
end
%%%%%%%%%%%%%%%%%%%%%%%

[c_max, c_min, c_limit]=Ma_max_min(s_all_c);
c_max_c=c_max;
c_limit_c=c_limit;

    if lock_scalebar==1
             
         if    era_name_nr==2 && (cor_nr==4 || cor_nr==6)
%                  c_limit= 1.7879;% 2mT EOF1 2014 max                  
            %       c_limit= 1.51;% 2mT EOF1 2014
                   %c_limit= 1.01;% 2mT EOF1 2012
                   c_limit= 0.75;%             
         elseif    era_name_nr==1 && cor_nr==4 && sea_nr==1
                        % c_limit= 50.01;
                      %   c_limit= 30.01;
                         c_limit= 25.01;
                         
         elseif    era_name_nr==3 && cor_nr==1 && yr_nr_s==2         
                         
                        c_limit= 1.00; % 1.6590
         end
      
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig

fig('units','inches','width',10,'height',9,'font','Helvetica','fontsize',16,'border','on');

   axesm( proj,'MapLatLimit',[lat1 lat2],'Grid','on','ParallelLabel','on','Frame','on',... %
       'MeridianLabel','on','FontWeight','bold','FontSize',18,...  % font size of degree labels
       'mlabellocation',[0:30:180,0:-30:-180]); 

  
set(gca,'box','off','XColor',[1,1,1],'YColor',[1,1,1]);  % turns off border box and axes  
hSurf=surfm(double(era_lat),double(era_long),squeeze(s_all_c(:,:,1))');
hold on
 
        if corr_label==1
            br=colorbar('FontWeight', 'bold', 'FontSize',18);   
            pos=get(br,'Position');  
        
            if cor_nr==4 && era_name_nr==1
            
            pos(1)=pos(1)+0.04;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar 
            pos(2)=pos(2)+0.1;      
            pos(3)=pos(3)+0.014;  % widthcolorbar
            pos(4)=pos(4)-0.195;  % height colorbar 
            
            
            elseif cor_nr==4 && era_name_nr==2 || cor_nr==6 && era_name_nr==2
                
            
            pos(1)=pos(1)-0.03;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar 
            pos(2)=pos(2)+0.15;      
            pos(3)=pos(3)+ 0.0;  % widthcolorbar
            pos(4)=pos(4)-0.315;  % height colorbar    
                
        
            end
        
        set(br,'Position',pos)     
        
        end
        
    
      colormap(b2r(-c_limit,c_limit));
      
      
      if era_name_nr==1 % z500 %%%%%%%%
        color_alt=6;
      elseif era_name_nr==5
        color_alt=1;
      else
        color_alt=1;    
      end
      
      
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
        elseif  color_alt==7   
        load('C:\Users\Machine\matlab_storage_of_output_files\cmap_colormap_save_c2.mat')
        colormap(cmap_colormap_save_c2)     
            
             
        end
      
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % p=0.05 significance level contour      
      
    if cor_nr==4 || cor_nr==5   || cor_nr==6 || cor_nr==3 || cor_nr==1    
        co_pa=[.1 .1 .1];
        
    else
        co_pa=[.9 .9 .9];
    end


 contour_alt_nr=2;   
    
if contour_alt_nr==1
    h1= contourm( double(era_lat),double(era_long),squeeze(p_all(:,:,1))', [0.05],'--','ShowText','off',...
    'Linecolor',co_pa,'LineWidth', 2);
elseif contour_alt_nr==2 % Wilks test
    h1= contourm( double(era_lat),double(era_long),p_fdr_rej', [0.05],'--','ShowText','off',...
    'Linecolor',co_pa,'LineWidth', 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coast_nr=1;%%%%(1/0) on/off %%%%% Coast lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if coast_nr==1

    color_nr=1;
    if color_nr==1
       color_code=[.1 .1 .1];   
    elseif color_nr==2
       color_code=[.6 .6 .6];  
    end

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
plot3m(lat_cr,lon_cr,'-k','LineWidth', 1.5,'color',color_code);
%plot3m(lat(in_c),long(in_c),'.k','MarkerSize', 1);


%addpath C:\Users\benman\matlab_mapping
bedmap2('gl','k','LineWidth', 1.5,'color',color_code);
bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% label above colorbar
label_size=22;
  
if cor_nr==1 && era_name_nr==3 %
    a1=axestext_c(1.101, +0.33, ['(',label_1,')'] );    
elseif cor_nr==4 && era_name_nr==1
    a1=axestext_c(0.04, +0.028, ['(',label_1,')'] ); 
elseif   era_name_nr==2  
    a1=axestext_c(0.15, +0.10, ['(',label_1,')'] );  
else
    a1=axestext_c(0.1, -0.04, ['(',label_1,' ', label_2,')'] );
end

set(a1,'FontWeight','bold','FontSize',label_size,'rotation',90);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variance explained by EOF
        if cor_nr==4 && era_name_nr==1 
            a2=axestext_c(0.95, +0.98, [num2str(variance_exp),'%'] );
            set(a2,'FontWeight','bold','FontSize',label_size+4);
        elseif cor_nr==6 && era_name_nr==2 && PC_nr<=3    
            a2=axestext_c(0.9, +0.90, [num2str(variance_exp),'%'] );
            set(a2,'FontWeight','bold','FontSize',label_size+4);            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lat long labels. Long  at this lat and placement of lat         
    if era_name_nr==1
        mlab_par=-28;
    elseif era_name_nr==2    
        mlab_par=-38; 
     %   mlab_par=-42; 
    else
        mlab_par=-10;
    end
        
mlabel('MLabelParallel',mlab_par,'PLabelLocation',[-75 -60 -45 -30 -15  0],'PLabelMeridian',100) 


%%%%%%%%%%%%%%%%%%%%%%
%  Title
show_title=1;
if show_title==1 && PC_nr==3 && era_name_nr==1 || show_title==1 && PC_nr==2 && era_name_nr==2 || show_title==1 && PC_nr==13 && era_name_nr==2
t4=title(title_str);       
pos=get(t4,'Position');
               
% move title        
    if cor_nr==1
        pos(1)=pos(1)+0.130;
        pos(2)=pos(2)-0.030;
    elseif cor_nr==4 && era_name_nr==1
        pos(1)=pos(1)+1.65;
        pos(2)=pos(2)-0.30;
    elseif cor_nr==4 && era_name_nr==2 || cor_nr==6 && era_name_nr==2 
        pos(1)=pos(1)+1.36;
        pos(2)=pos(2)-0.59;      
        

    end
        set(t4,'Position',pos,'FontSize',26, 'FontWeight', 'bold');  

end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Letter
   if era_name_nr==1
        letter_pos= [170 710 50 62];
   elseif era_name_nr==2
        letter_pos= [230 650 50 62];
   end
        letter_size=38;
    

          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',letter_size ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
     
     %%%%%%%%%%%%%%%%%
     % colorbar fontsize 
     if cor_nr==4
        scalebar_fontsize=22;    
     else
        scalebar_fontsize=16;
     end
     
      set(br, 'FontSize',scalebar_fontsize, 'FontWeight', 'bold'); 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
     % save fig
        str_end='_c';
     
    
        if cor_nr==4 || cor_nr==6      
            filename=['ERA_interim_regression_', era_name_str,'_',PC_nr_str,'_',season,...
                '_',type_str,'_',num2str(yr_s),'-',num2str(yr_e),str_end];
        
        end
        

filedir ='C:\Users\Machine\matlab_storage_of_output_files\figures\';

if corr_label==1
savefilename_c=strcat(filedir,filename,'_crop');
else
savefilename_c=strcat(filedir,filename);

end


% save as png 
orient landscape
cd('G:\My Drive\ISO_CFA\matlab')
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r240',savefilename_c); % PNG-nocrop'
export_fig('-pdf','-painters', '-depsc','-opengl', '-r240',savefilename_c); % PNG-nocrop'
cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab')
toc
