%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DE HadISST_rice_corr_annual_seasonal_DE20_c8.m
% Daniel Emanuelsson, 2021
% 
% GitHub version 1  - 04-09-2018 
% GitHub version 2  - 23-08-2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%
% * corrcoef_df.m   UoW Steig
% * annave.m        UoW online archive Atmospheric Science
% * star_coord_WA.m DE
% * cosweight.m     UoW online archive Atmospheric Science
% * fig.m           fileexchange

% Wilks test
%%%%%%%%%%%%%%%%%%%%%%
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

 eval('Config_corr_SIC')
 %eval('Config_corr_SIC_suppl_m')
 
end
toc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

% Input

% param_nr=2;% (1-SST/ 2 SIC)

% SST_dataset=1; % (1-HadISST, 2-ERSST) 

rcontour_psa=0; %(0)(1) show where PSA pattern overlap interferance (2) shows where PSA2 is significant

% coast_nr=1;  % On/Off
    
% if param_nr==2
%     name='SIC';
% end
%  


  
%%%%%%%%%%%%%%%%%%%%%
% sea_nr=3;
if sea_nr==1    
    season='annual';
elseif sea_nr==2
    season='MAMJJASON';
elseif sea_nr==3    
    season='AMJJASON';
end
%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%%%%%%%%%%%%%
 show_title=4; % (1)-long (2)-short (3) use for Interval SIC figure  (4) r (x,y)
  
        
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Int_in=0; %%% turn (1/0) on/off, turn on for Fig 10a, b 
%  
%  interval_nr=2; %(2, 3, 4)%%%%%%%%%%%%%%%%%%% 2,3 in use for both SST and SIC
%  
%  if  Int_in==1
%      
%        alt_pc=2;
%     
%     if  interval_nr==2 %%%%%%%%%%%%%%%%%% used change name, in-phase
% 
%         yr_s=1994;
%         yr_e=2001;
%         letter='b';
%         show_title=3;
%                 
%     elseif  interval_nr==3 %%%%%%%%%%%%%%%%%% used change name, not in-phase
% 
%         yr_s=1979; % 
%         yr_e=1993;
%         letter='a';
%         show_title=3;
%         
%     elseif  interval_nr==4 % out-of-phase
%       
%         yr_s=2005; % 
%         yr_e=2011;
%         letter='c';
%         show_title=3;       
%     end
%     
%  else
%     
%        alt_pc=0; % dont show stippling
% 
%  end
 
%%%%%%%%%%%%%%%%
 
 show_colorbar=1;
 iso='d18O';
 show_maximum_point=0;
 p_level_rc=0.05; % p level choose 0.1 for figures 0.05 for SST time series
 
%%%%%%%%%%%%%%%%%%%%%%%  
% if  strcmp(name,'SIC')==1
%   proj='stereo';
%   lat1=-90;
%   lat2=-50;
%   box_use=0;  % (0/1)  
% end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HadISST data
% HadISST data
% version 1.1 Rayner et al. 2003
%
%           standard_name = 'sea_ice_area_fraction'
%           long_name     = 'Monthly 1 degree resolution sea ice concentration'

%addpath C:\Users\Machine\matlab_lib\Data\HadISST_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    name_c=[data_dir,'HadISST_ice_c.nc'];
    HadISST_time=ncread(name_c,'time'); %units         = 'days since 1870-1-1 0:0:0'
    % missing values  missing_value = -1e+30 ( looks more like -1000)\
    %monthly values
 
    HadISST_time_c=HadISST_time+datenum(1870,1,1);
 
    date_vec=datevec(double(HadISST_time_c)); 
    yyyy=date_vec(:,1); 
    mm=date_vec(:,2);
    HadISST_year_num=yyyy+(mm/12-1/12); 

    HadISST_time_bnds =ncread(name_c,'time_bnds');
    HadISST_lat=ncread(name_c,'latitude');
    HadISST_lon=ncread(name_c,'longitude');
    
    
         %   name='ice'; 
        HadISST_ice=ncread(name_c,'sic'); % 360x180x1728 lon x lat x time
        M=HadISST_ice;   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


HadISST_year_num(HadISST_start);
HadISST_year_num(HadISST_count);

% HadISST_year_num()-HadISST_year_num()
% HadISST_start:HadISST_count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ===================================================
 HadISST_lat_c=HadISST_lat;
 [nlon, nlat, ntim]= size(M);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  era_z500_c = reshape(era_z500, ntim, nlat, nlon);  %

HadISST_z500_c = permute(M, [3 2 1]);  %  z500_eof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 %%% -=(1)=-
% Remove monthly climatology for all cells 

 HadISST_z500_c2 = double(reshape(HadISST_z500_c, ntim, nlat*nlon));  % One time series (column) for each grid point
 
 %  lat 1, 2, 3, 4,..........33...lat 1, 2, 3, 4.....33       until 33x240=7920
[HadISST_z500_c3,clim_z500] = annave(HadISST_z500_c2);   % Removes seasonal cycle  

HadISST_z500_c4= reshape(HadISST_z500_c3,  ntim, nlat, nlon );

%%%%%%%%%%%%%%%
%%% -=(2)=-
%  z = cosweight(z, lat);

HadISST_z500_c5=cosweight(HadISST_z500_c4,HadISST_lat_c); % Original UoW function    (time x lat x lon) 

%%%%%%%%%%%


% Back to old format again
HadISST_z500_c6= reshape(HadISST_z500_c5,  ntim, nlat, nlon );

M = permute(HadISST_z500_c6, [3 2 1]);  %  z500_eof
 
%==========================================================
% Annual and seasonal means
% (360x180x1728) (Long,Lat,month)
[m n t]=size(M(:,:,HadISST_start:HadISST_count));% time(months) % 1894-2011
for i= 1:m
    for j= 1:n

        dummy1=reshape(squeeze(M(i,j,HadISST_start:HadISST_count)), 12, div);

   
if strcmp(season,'annual')==1
         HadISST_M_annual_temp=nanmean(dummy1,1); % annual HadISST-values

         
%          HadISST_M_annual_temp=anomaly(HadISST_M_annual_temp); % try
%          normalizing
         HadISST_M_annual_dummy=detrend(HadISST_M_annual_temp(1:end));        
         HadISST_M_annual(i,j,:)=HadISST_M_annual_temp(1:end);
         HadISST_M_annual_detrend(i,j,:)= HadISST_M_annual_dummy; % detrended
         
elseif strcmp(season,'MAMJJASON')==1         
         
          HadISST_M_season_mamjjason=nanmean(dummy1(3:11,:),1);  %3. seasonal (JJA)
          HadISST_M_season_mamjjason=HadISST_M_season_mamjjason(1:end);
          HadISST_M_season_reg_mamjjason(i,j,:)=HadISST_M_season_mamjjason; % means for regression         
          HadISST_M_season_dummy_mamjjason=detrend(HadISST_M_season_mamjjason); 
          HadISST_M_season_detrend_mamjjason(i,j,:)=HadISST_M_season_dummy_mamjjason;
          
elseif strcmp(season,'AMJJASON')==1            
          HadISST_M_season_amjjason=nanmean(dummy1(4:11,:),1);  %3. seasonal (JJA)
          HadISST_M_season_amjjason=HadISST_M_season_amjjason(1:end);
          HadISST_M_season_reg_amjjason(i,j,:)=HadISST_M_season_amjjason; % means for regression         
          HadISST_M_season_dummy_amjjason=detrend(HadISST_M_season_amjjason); 
          HadISST_M_season_detrend_amjjason(i,j,:)=HadISST_M_season_dummy_amjjason;
    
    
end

%   clear       dummy1 HadISST_M_annual_dummy HadISST_M_annual_temp

    end
end
toc
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if strcmp(iso,'dD')==1 || strcmp(iso,'d18O')==1
    
% isotopes record
%date_annual=[1900:2011]';

%    iso_alt_nr=3;
        
        if iso_alt_nr==1
             load([filedir,'RICE_combined_Deep_1213B_c23.mat']);
             iso_season_str='annual';
        elseif iso_alt_nr==2
            load([filedir,'RI_combined_Deep_1213B_annual_means_no_summer_c24.mat']); % MAMJJASON 
%             load('C:\Users\benman\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_clim_c24.mat'); % clim removed, no difference
            iso_season_str='ExtendedWinter';
        elseif iso_alt_nr==3    
            load([filedir,'RI_combined_Deep_1213B_annual_means_no_summer_AMJJASON_c24.mat']); % AMJJASON 
%             load('C:\Users\benman\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_AMJJASON_clim_c24.mat'); % clim removed, no difference
            iso_season_str='ExtendedWinter';           
        end

date_annual=MA_save(:,1);    
start_t=find(date_annual==yr_s);
end_t=find(date_annual==yr_e); 

    stacked_record_annual_Ma=MA_save;
%    X=detrend(stacked_record_annual_Ma((start_t:end_t),2));
    X=detrend(stacked_record_annual_Ma((start_t:end_t),3)); % d18O
end 
   
if strcmp(season,'annual')==1
    %season='Annual';
    Y=HadISST_M_annual_detrend; 
 
elseif strcmp(season,'MAMJJASON')==1   
   Y=HadISST_M_season_detrend_mamjjason;
   
elseif strcmp(season,'AMJJASON')==1   
   Y=HadISST_M_season_detrend_amjjason;    
end


[A B C]=size(Y);
R=zeros(2,2,A,B);P=zeros(2,2,A,B);
for i=1:A
    for j=1:B
        Q=squeeze(Y(i,j,:));
        
            
        [R(:,:,i,j) P(:,:,i,j)]=corrcoef_df(X,Q, 'rows','pairwise');
        

    end
end
% % for RI
% R(:,:,find(floor(HadISST_lon)==199),find(floor(HadISST_lat)==-80))
% P(:,:,find(floor(HadISST_lon)==199),find(floor(HadISST_lat)==-80))
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Wilks test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contour_alt_nr=2; % (1), regular p (2) including Wilks test

if contour_alt_nr==2

p_in=squeeze(P(1,2,:,:));
%p_fdr_rej=fdr_wilks_test(p_in); % made for 360x180 grid 


% Rw=squeeze(R(1,2,:,:));

p_in=p_in(:,(146:167)); % SH sea ice edge

% A_nan=isnan(p_in);
% A_nan_nr=sum(sum(A_nan));


%  wr_c=nanmean(HadISST_M_annual_detrend,3);
 wr_c=nanmean(Y,3);
 
 wr_c=wr_c(:,(146:167));
 wr_c2=wr_c; 
 wr_c2(isnan(wr_c))=0; % turn nans to zero
 
A_ind2=find(wr_c2~=0);
A_nan_nr2=length(A_ind2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from Github 
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
[rn_x,rn_y] = ind2sub([360,17],rej_null(rej_null_ind));
p_local_rej = ones(360,22);
p_fdr_rej = ones(360,22);
%no_grid_points = 360*180;
%no_grid_points = A_nan_nr2; % fewer grids in sea ice data  
no_grid_points = A_nan_nr2;
n_fdr = 0; 

%loop over all the cases where it the null hypothesis
%is assumed could be rejected locally
for i = 1:no_rej_null
    pi_c = p_in(rej_null(rej_null_ind(i)));
    p_local_rej(rn_x(i),rn_y(i))= pi_c;
    % compute fdr threshold
    %fdr_thres = 2*alp_c*i/no_grid_points;
      fdr_thres = alp_c*((no_rej_null+1)-rej_null_ind(i))/no_grid_points; % remove 2, eq. 3 Wilks 2016
    
    %find local p values that are less than fdr threshold
    if(pi_c<fdr_thres)
        p_fdr_rej(rn_x(i),rn_y(i)) = pi_c;
        n_fdr = n_fdr+1;
    end 
end
 toc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Figure  
 tic
  lock_scalebar=1;%%%%%%%%%%%%%%%%%%%%%%%
   
  if Int_in==1
    RSAS_SIC_box=0;    
  elseif Int_in==0    
    RSAS_SIC_box=0; % (1) area defined using RICE corr pattern (2) ADP following Yuan 2004
  end
  
 %%%%%%%

  f2=fig('units','inches','width',10.5,'height',10.5,'font','Helvetica','fontsize',16,'border','on'); 

  axesm( 'stereo','MapLatLimit',[lat1 lat2],'Grid','on',...
      'Frame','on','ParallelLabel','on',...
       'MeridianLabel','on','FontWeight','bold','FontSize',18, 'mlabellocation',[0:-30:-360, 0:30:180]); 
 
  set(gca,'box','off','XColor',[1,1,1],'YColor',[1,1,1]);  % turns off border box and axes 

  hold on

  set(gca,...
     'linewidth',2,...
     'FontWeight','bold' );
      
 gridm('GLineStyle','--','Gcolor',[.5 .5 .5],'GLineWidth',1.5,... %%%%%%%%%% Grid
    'MLineLimit',[lat1 lat2],'Galtitude', .02)


Rc=R;
Rc5=R;
% Rc2=R;

% Rc(P>=p_level_rc)=NaN; 
% Rc2(P>=0.01)=NaN; % Mask out above or equal R's
Rc5(P>=0.1)=NaN; 


c_min=squeeze(R(1,:,:));
c_min=c_min(2,:);
c_min=min(c_min);

c_max=squeeze(R(1,:,:)); 
c_max=c_max(2,:);
c_max=max(c_max);

c_limit=max(abs(c_min),abs(c_max));
c_limit_c=c_limit;


if strcmp(name,'SIC')==1

 X_r=squeeze(Rc(1,2,:,:))';

X_r(:,(100:360))=NaN;  %!!!!!!!!!! Look ONLY in RS ABS region,  nan out around east antarcica 
X_r((1:90),:)=NaN; % southern Hem
 
%  X_r_max=nanmax(nanmax(X_r));   
  X_r_max=nanmin(nanmin(X_r));  
 
X_r_long=nanmin(X_r); % 1 long
xr=find(X_r_long==X_r_max); 
xr_lon=HadISST_lon(xr);

X_r_lat=nanmin(X_r,[],2); %lat
xr=find(X_r_lat==X_r_max); 
xr_lat=HadISST_lat(xr);
 
   
end


if lock_scalebar==1 
   c_limit=0.701; 
    
end

      mlabel('MLabelParallel',-54  ,'PLabelLocation', [-75 -60 -45 -30],'PLabelMeridian',100) 
      
%%%%%%%%%%%%%%%% Fill in seam  % only one side
% HadISST_lat_c=HadISST_lat;
HadISST_lon_c=HadISST_lon;
HadISST_lon_c(361,:)=179.99;

Rc6=squeeze(Rc5(1,2,:,:))';
Rc6(:,361)=Rc6(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hSurf=surfm(double(HadISST_lat),double(HadISST_lon_c),squeeze(Rc6));
colormap(b2r(-c_limit,c_limit));
%colormap(jet);
% colormap(brewermap(256,'*RdBu'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Int_in==1 % color for Fig. 10
   color_alt=1;   
else
    color_alt=6;    
end


if color_alt==1
    colormap(brewermap(256,'*RdBu'));
elseif color_alt==2
    colormap(flipud(cbrewer2('div','Spectral',10)));
elseif color_alt==3
    colormap(flipud(cbrewer2('div','PiYG',20)));      
elseif color_alt==4
    colormap(flipud(cbrewer2('div','RdGy',20)));      
elseif color_alt==5
    colormap(flipud(cbrewer2('div','RdBu',20)));    
elseif color_alt==6
    colormap(cbrewer2('div','BrBG',20));
elseif color_alt==7
    colormap(flipud(cbrewer2('div','RdYlBu',20)));    
end

%%%%%%%%%%%%%%%%%%%%%%
% hSurfAx=(gca);
cRange= caxis;
% contourm function keeps projection defined above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if contour_alt_nr==1
    h1= contourm( double(HadISST_lat),double(HadISST_lon),squeeze(P(1,2,:,:))', [0.05],'-k','ShowText','off',...
   'LineWidth', 2);

 elseif contour_alt_nr==2 % Wilks test

    h1= contourm( double(HadISST_lat(146:167)),double(HadISST_lon),p_fdr_rej', [0.05],'-k','ShowText','off',...
   'LineWidth', 2);
 end


% Then reset the color axis according to the range determined above:
caxis(cRange); 
alpha(1) 
%%%%%%%%%%%%%%%%%%%%%%%%
if RSAS_SIC_box==1    
        
        lon_b_c=[-120 -120 -150 -150 -120];
        lat_b_c=[ -73 -60 -60 -73 -73];
        plotm(lat_b_c,lon_b_c,'-','LineWidth',4,'Color',[1 .7 .1]) 
        
 elseif RSAS_SIC_box==2         % ADP
        
        lon_b_c=[-130 -130 -150 -150 -130];
        lat_b_c=[ -70 -60 -60 -70 -70];
        plotm(lat_b_c,lon_b_c,'-','LineWidth',4,'Color',[1 .2 .2])  % 1 .7 .1
        
        lon_b_c=[-40 -40 -60 -60 -40];
        lat_b_c=[ -70 -60 -60 -70 -70];
        plotm(lat_b_c,lon_b_c,'-','LineWidth',4,'Color',[1 .2 .2])         
        
end        
%%%%%%%%%%%%%%%%%%%%%%%%

hSurfAx=(gca);
cRange= caxis;

if show_colorbar==1
  %h=colorbar('SouthOutside'); 
  h=colorbar('EastOutside');  
%   set(h, 'Position', [.76 .145 .015 .75], 'FontSize',18, 'FontWeight', 'bold');  
end

    
if yr_e==2009 % suppl. mat

    a1=axestext_c(1.1,+0.6, ['Correlation'] ,'rotation',-90,'FontSize',22, 'FontWeight', 'bold');    %

else
    
    a1=axestext_c(1.09,+0.61, ['Correlation'] ,'rotation',-90,'FontSize',24, 'FontWeight', 'bold');    %

end 
    
    
    
    
    % modifiy colorbar
    pos=get(h,'Position');
    
%     if yr_e==2009
    
        
%         pos(1)=pos(1)+0.024;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar   (+) --->>   +0.02;
        pos(1)=pos(1)+0.015; % 0.035
        pos(2)=pos(2)+0.11;  
        pos(3)=pos(3)+ 0.012;  % widthcolorbar
        pos(4)=pos(4)-0.22;  % height colorbar
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        
%     else
%         
%         pos(1)=pos(1)+0.028; % 0.035
%         pos(2)=pos(2)+0.14;  
%         pos(3)=pos(3)- 0.0052;  % widthcolorbar
%         pos(4)=pos(4)-0.31;  % height colorbar       
%         
%         
%         
%     end
        
        set(h,'Position',pos)
        
h.TickLabels = strrep(h.TickLabels, '-', '–');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

if strcmp(iso,'dD')==1
    iso_str='{\delta}D';
elseif strcmp(iso,'d18O')==1
    iso_str='{\delta}^1^8O';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% letter
TextBox = uicontrol('style','text');
if  strcmp(proj,'merc')==1          
            
  letter_position=  [180 500 50 55];% [200 635 50 50];
  
    %set(TextBox,'String',letter,'position',[130 820 50 50],'FontWeight','bold','FontSize',28 ); % x position y position xsize ysize  
    set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',28 ); % x position y position xsize ysize            
            
elseif  strcmp(proj,'stereo')==1
    
       if strcmp(name,'SIC')==1 
 %      letter_position=  [200 850 60 60];%
       letter_position=  [200 810 60 60];% 
       else
       letter_position=   [330 840 50 50];
       end
    
        
       
       if yr_e==2009
            set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',44 ); % x position y position xsize ysize 
       else
           set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',45 ); % x position y position xsize ysize
       end
            
            
            %285
end
            set(TextBox,'foregroundcolor', [0 0 0], ...
            'backgroundcolor', [1 1 1]);
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 coast_nr=1;  % On/Off
if coast_nr==1
    
    load coast
    %%%%%%
    % to be able to use coastline for continents in combination with ant
    % coastline and grounding line from bedmap
    in_c=find(lat<-60);
    lat_cr=lat;
    lon_cr=long;
    lat_cr(in_c)=NaN;
    lon_cr(in_c)=NaN;

    %%%%%%%%%%%%%

%    addpath C:\Users\benman\matlab_mapping
        if strcmp(name,'SIC')==1
            color_code=[.4 .4 .4];    
            plot3m(lat_cr,lon_cr,'-','LineWidth', 2,'color',color_code);
            bedmap2('gl','LineWidth', 2,'color',color_code);
        end

    bedmap2('gl','LineWidth', 2, 'color',color_code );
    bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if strcmp(name,'SIC')==1 && show_maximum_point==1     
    plotm(double(xr_lat),double(xr_lon),'.','color',[.4,.8,.2],'MarkerSize',36);  % max corr point

         TextBox = uicontrol('style','text');
         set(TextBox,'String',[' r = ',num2str(round(X_r_max*1000)/1000)],...
             'position',[310 130 170 50],'FontWeight','bold','FontSize',22 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);    
end

%%%%%%%%%%%%%%
% marker for ice-core site
if  strcmp(site,'RICE')==1
        site_coor=star_coord_WA(site); % RICE marker
        % Circle Alt
        plotm(site_coor(1),site_coor(2),'.','MarkerSize',36,'MarkerFaceColor','none','MarkerEdgeColor',[.9,.1,.9]); %,'LineWidth',2.5); % [.4,.6,.1]
        %%%% Dot Alt
        %plotm(site_coor(1),site_coor(2),'.','MarkerSize',40,'MarkerFaceColor','none','MarkerEdgeColor',[.9,.1,.9],'LineWidth',2.5);   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title
if show_title==1 || show_title==2 || show_title==3 || show_title==4
                                    SST_dataset_str='HadISST';   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if show_title==1                       
                        hp1=  title([site,' ',iso_str,' ',name,' ',SST_dataset_str,' ',season,...
                            ' ', num2str(round(date_annual(start_t))),' - ',num2str(round(date_annual(end_t)))],'FontSize',20, 'FontWeight', 'bold');
        elseif show_title==2
                        hp1=  title([iso_str,', ',name],'FontSize',28, 'FontWeight', 'bold');                                                
        elseif show_title==3
                        hp1=  title([iso_str,' ',name,' ',SST_dataset_str,' ',...
                            ' ', num2str(round(date_annual(start_t))),' - ',num2str(round(date_annual(end_t)))],'FontSize',20, 'FontWeight', 'bold');
        elseif show_title==4
                        hp1=  title(['r(',iso_str,', ',name,')'],'FontSize',28, 'FontWeight', 'bold');  
        end

 % move title       
        if strcmp(proj,'stereo')==1  
            pos=get(hp1,'Position');
            pos(2)=pos(2)-0.07;
            set(hp1,'Position',pos,'FontSize',30, 'FontWeight', 'bold');    
%         elseif strcmp(proj,'mec')==1
%             pos=get(hp1,'Position');
%             pos(2)=pos(2)-0.02;
%         	set(hp1,'Position',pos,'FontSize',18, 'FontWeight', 'bold');    
        end
       
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if show_colorbar==1
        set(h, 'FontSize',20, 'FontWeight', 'bold'); 
    end
    
%%%%%%%%%%
% areas of EOF1 and EOF2 overlap%%%%%%%%%%%%%
% load contours                        
% if rcontour_psa==1
%             alt_pc=2 ;%%%%%% (1) z500-SIC reg pattern (2) PC-SIC patterns %%%%%%%%%
%             p_level_rcontour_psa=0.05; % 0.05 or 0.01
% 
%     if sea_nr==1
%                   if Int_in==0
%                      contour_set=1;%%%%%%%%%%%%%%%
%                   elseif interval_nr==2 || interval_nr==3 || interval_nr==4 
%                      contour_set=2;    
%                   end
%                                          
%     % load surfaces    
%         if contour_set==1
%             %PSA
%             filename_con=['ind_psa_sum_',num2str(p_level_rcontour_psa),'.mat'];
%         elseif contour_set==2
%         % SAM PSA1
%       
%                                   
%             if interval_nr==2 || interval_nr==4 
% %                 if param_nr==1 % SST
% %                     filename_con=['ind_sam_psa1_sum_ipo_p_',num2str(p_level_rcontour_psa),'.mat'];
%                  if param_nr==2 % SIC
%                   
%                     if alt_pc==1
%                         filename_con=['ind_SICsam_psa1_sum_ipo_p_',num2str(p_level_rcontour_psa),'.mat'];
%                     elseif alt_pc==2 
%                         filename_con=['ind_SICpc1_pc2_sum_ipo_p_',num2str(p_level_rcontour_psa),'.mat'];   %%%%%%%%%%%%%% 
%                     end
%                     
%                 end
%             elseif interval_nr==3 % IPO-
% %                 if param_nr==1 % SST
% %                     filename_con=['ind_sam_psa1_sum_ipo_n_',num2str(p_level_rcontour_psa),'.mat'];                
%                 if param_nr==2 % SIC
%                     if alt_pc==1
%                         filename_con=['ind_SICsam_psa1_sum_ipo_n_',num2str(p_level_rcontour_psa),'.mat'];
% %                       filename_con=['ind_SICpc1_pc2_sum_ipo_p_',num2str(p_level_rcontour_psa),'.mat']; % same as above change colorcode instead
%                     elseif alt_pc==2
%                         filename_con=['ind_SICpc1_pc2_sum_ipo_n_',num2str(p_level_rcontour_psa),'.mat'];    %%%%%%%%%%%%%%%%%%%
%                     end
%                 end
%             end
%         end 
%     end
% 
% folder_c='C:\PHD\matlab_storage_of_output_files\';
% load([folder_c,filename_con]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% %%%% overlap stippling strengthened/canceled out corr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if param_nr==2 && alt_pc==2 % Using  PC SIC regression patterns 
%     
%  %%%%%%%%%%%%%%%%%%%%  
%    active_regions_stipple=0; % (1/0) on/off turn off for checking things turn on for final fig (slow)
%     
% if active_regions_stipple==1    
%     
%    color_code2=[0 0.8 .2]; 
%    color_code=[1 .7 .1];   % Yellow 
%      [c1,c2,c3]=size(ind_psa_sum_c);
%  
%    for i=1:c1  %360 lon
%      for  k=140:166% c2; % lat
%    
%         
%         ind_cp=  find(ind_psa_sum_c(i,k)==2 | ind_psa_sum_c(i,k)==-2); %%%%%%% 2  -2 Yellow pattern reinforce one another
%         
%         if isempty(ind_cp)==1
%             ind_cp=9999; % assign if empty
%         end
%         
%          land_ind1 = landmask(double(HadISST_lat(k)), double(HadISST_lon(i)),'Antarctica'); % check if lat long falls on land
%          land_ind2 = landmask(double(HadISST_lat(k)), double(HadISST_lon(i)),'North and South America'); % check if lat long falls on land
% 
%          
%          
%          % plot colored stippling
%          % Reinforce
%          if ind_cp==1 && land_ind1==0 &&  land_ind2==0    
%             s1=plotm(double(HadISST_lat(k)), double(HadISST_lon(i)),'color',color_code,'Marker','.','MarkerSize',8);  
%          end
% 
%          %%%%%%%%%%%%%%%%
%          ind_cp2=  find(ind_psa_sum_c(i,k)==0); %%%%%%% 0 green pattern cancel one another
%          
%                  if isempty(ind_cp2)==1        
%                     ind_cp2=9999;
%                  end
%                  
%                  
%          if ind_cp2==1 && land_ind1==0 &&  land_ind2==0    
%          s2=plotm(double(HadISST_lat(k)), double(HadISST_lon(i)),'color',color_code2,'Marker','.','MarkerSize',8);  
%          end        
%          
%          
%      end
%    end
%    
%    
% end
%     
%  clear c1 c2 c3  
%     %%%%%%%%%%%%% 
% end
%  
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
% if   sea_nr==1 &&  param_nr==2 && alt_pc==2 
%         
%            color_code=[0.8 0.2 .8];   % not an active region
%            [C11, h11]= contourm( double(HadISST_lat),double(HadISST_lon),ind_psa_sum_nonactive',[0,0],'-','LineStyle','none');
%            [cp_11x, cp_11y, cp_11z]=C2xyz (C11);
% 
% ind_psa_sum_nonactive_c=ind_psa_sum_nonactive;
% ind_psa_sum_nonactive_c(ind_psa_sum_nonactive_c==0)=NaN;
% ind_psa_sum_nonactive_c(ind_psa_sum_nonactive_c==2)=1;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% stipple_nr=0; % (0/1) non-active areas
%  
%     if  stipple_nr==1
%  
%     [c1,c2,c3]=size(ind_psa_sum_nonactive_c);
%  
%         for i=1:c1  %360 lon
%             for  k=140:166% c2; % lat
%          
%          
%             ind_cp=isnan(ind_psa_sum_nonactive_c(i,k));
%          
%             land_ind1 = landmask(double(HadISST_lat(k)), double(HadISST_lon(i)),'Antarctica'); % check if lat long falls on land
%             land_ind2 = landmask(double(HadISST_lat(k)), double(HadISST_lon(i)),'North and South America'); % check if lat long falls on land
%          
%             if ind_cp==1 && land_ind1==0 &&  land_ind2==0
%                 s1=plotm(double(HadISST_lat(k)), double(HadISST_lon(i)),'k','Marker','.','MarkerSize',6);  
%             end
%      
%      
%             end
%         end
%  
%     end
%         
% end
% 
% %%%%%%%%%%%%alt v.%%%%%%%%%%%
% elseif rcontour_psa==2
%     p_level_rcontour_psa=0.05; % 0.05 or 0.01
% 
%     if sea_nr==1
%         filename_con2=['HadISST_regression _SST_PSA2_monthly_',num2str(p_level_rcontour_psa),'.mat'];
%         filename_con1=['HadISST_regression _SST_PSA1_monthly_',num2str(p_level_rcontour_psa),'.mat'];
%         
%     end
%     
%     folder_c='C:\PHD\matlab_storage_of_output_files\';
%     load([folder_c,filename_con2]);
% 
%     s_all_c_psa2=s_all_c;  % this is the p-level surface as rs with lower p-levels are masked above
%     ind_psa2_nonnan=~isnan(s_all_c_psa2);
%     
%     color_code=[1 .7 .1];                  
%     h2= contourm( double(HadISST_lat),double(HadISST_lon),double(ind_psa2_nonnan)',[1,1],'-k','ShowText','off',...
%     'Linecolor',color_code,'LineWidth', 2.5); 
% 
% end
hold off
%%%%%%%%%%%%%%%%%%%%
% save figure
% filedir ='C:\Users\Machine\matlab_storage_of_output_files\figures\'; 
    if rcontour_psa==1
                filename=strcat(filedir,'figures\','HadISST_',name,'_',iso,'_',season,'_',proj,'_',num2str(yr_s),'_',num2str(yr_e),'_wPSA_pattern');
    else
                filename=strcat(filedir,'figures\','HadISST_',name,'_',iso,'_',season,'_',proj,'_',num2str(yr_s),'_',num2str(yr_e));
    end

savefilename_c=filename;

 % save as png 
orient landscape

cd('G:\My Drive\ISO_CFA\matlab')

export_fig('-pdf','-painters', '-depsc','-opengl', '-r190',savefilename_c);  %func needs ghostscript ,'-nocrop' 
export_fig('-png','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG  110  '-nocrop',

cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab')  
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%