%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Composite differences
% composite_SAM_PSA1_in_phase_ERA_I_z500_2mT_HadISST_SIC_DE21.m
% Read in ERA-interim data
% also HadISST SIC data at the end 
%%%%%%%%%%%%%%%%%%%%%%%%%
% Daniel Emanuelsson, 2021
% MATLAB 2020b
% t-test
% FDR-test Wilks 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    eval('Config_inPhase')
end
toc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_c=[data_dir,'ERA_int_monthly_z500_2.nc'];
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
eof_type=3;

if eof_type==1

eof_str='';
elseif eof_type==2
%eof_str='rotate';
eof_str='varimax';

elseif eof_type==3
    eof_str='varimax';
    %load('C:\Users\benman\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2011_MAMJJASON_c52_annual_mean_varimax_c2');
    load('C:\Users\Machine\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2011_AMJJASON_c52_annual_mean_varimax_c2');
end


time_c=MA_PCs_save(:,1);
            
            
era_count=find(era_year_num==yr_e)+11;  % data in months
era_start=find(era_year_num==yr_s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mm_in_c=[era_start:era_count]';
yr_1=find( time_c==yr_s);       
yr_2=find( time_c==yr_e);    
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             PC_nr=2; % (2) SAM (3) PSA1 (4) PSA2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (1) is just time
            
            
%          if PC_nr==2  % SAM
            
                if eof_type==2
               % fact_nr=1;  %2014
                fact_nr=-1;  %2011 eof  
                elseif eof_type==3
                fact_nr=1;
                end
%             elseif PC_nr==3 % PSA1

%                  if eof_type==2
%                     fact_nr=1; %2014 and 2011 eof
%                  elseif eof_type==3
%                     fact_nr=-1; 
%                  end
%               
%             elseif PC_nr==4 % PSA2
%                  fact_nr=-1;  %2014  To get sign as in Kidson 1988 paper

%         end            
            
            
            
            label_2='';
            iso='z500';
            
            
    %   y=MA_PCs_save_monthly((mm_in_c),PC_nr)*fact_nr;  % monthly
       SAMi=MA_PCs_save((1:33),PC_nr)*fact_nr;  
       
      %plot(time_c,nanmoving_average(y,5),'k-')
      f1=fig('units','inches','width',18,'height',5,'font','Helvetica','fontsize',16,'border','on');
      h1= plot(MA_PCs_save((1:33),1),SAMi(1:33),'k-','LineWidth',2.5);% annual
      hold on
      box on
         fact_nr2=1; 
     %  y2=MA_PCs_save_monthly((mm_in_c),3)*1;  
       PSA1i=MA_PCs_save((1:33),3)*fact_nr2;
       
%        plot(time_c,nanmoving_average(y2,5),'r-')
       h2=plot(MA_PCs_save((1:33),1),PSA1i(1:33),'r-','LineWidth',2.5);
       

% threshold 
       %perc_std=1;
       perc_std=0.5;
       
       SAM_lim=perc_std*nanstd(SAMi);
       PSA1_lim=perc_std*nanstd(PSA1i);
       
    lrc3=hline([0],':');
%   set(lrc3,'Color',[0.8,0.5,0.1],'LineWidth',1.5);  
  set(lrc3,'Color',[0.1,0.8,0.1],'LineWidth',1.5);     
     
  lrc3=hline([SAM_lim,-SAM_lim],':');
%   set(lrc3,'Color',[0.8,0.5,0.1],'LineWidth',1.5);  
  set(lrc3,'Color',[0.1,0.1,0.1],'LineWidth',1.5);
  
  
  
  lrc3=hline([PSA1_lim,-PSA1_lim],':');
%   set(lrc3,'Color',[0.8,0.5,0.1],'LineWidth',1.5);  
  set(lrc3,'Color',[0.9,0.0,0.0],'LineWidth',1.5);  
  
  
                    set(gca,...
          'linewidth',2,...
          'FontWeight','bold' );
         
%correlation
%  factor_nr
[r,p]=corrcoef_df(SAMi(11:21),PSA1i(11:21)); % 1989-2000
rp=[r(2),p(2)]
p_level(p(2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAM index limits
       
SAMi_c_p=SAMi;
SAMi_c_p(SAMi<SAM_lim)=0;
SAMi_c_p(SAMi_c_p>0)=1;

SAMi_c_n=SAMi;
SAMi_c_n(SAMi>-SAM_lim)=0;
SAMi_c_n(SAMi_c_n<0)=-1;

SAMi_c_c=SAMi_c_p+SAMi_c_n; % combined


% PSA1 index limits
PSA1i_c_p=PSA1i;
PSA1i_c_p(PSA1i<PSA1_lim)=0;
PSA1i_c_p(PSA1i_c_p>0)=1;

PSA1i_c_n=PSA1i;
PSA1i_c_n(PSA1i>-PSA1_lim)=0;
PSA1i_c_n(PSA1i_c_n<0)=-1;

PSA1i_c_c=PSA1i_c_p+PSA1i_c_n; % combined  


ind_c1=find((SAMi_c_c==1) & (PSA1i_c_c==-1)); % SAM+, PSA1-

MA_PCs_save(ind_c1,1) % years

plot(MA_PCs_save(ind_c1,1),0,'k*','MarkerSize',14,'linewidth',1.2 ) % years


ind_c2=find((SAMi_c_c==-1) & (PSA1i_c_c==1)); % SAM-, PSA1+

MA_PCs_save(ind_c2,1) % years

plot(MA_PCs_save(ind_c2,1),0,'ko','MarkerSize',14,'linewidth',1.2 )


set(gca,'XLim',[yr_s-1  yr_e+1])  
xlabel('Year','FontWeight','bold','FontSize', 20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% save index

end_str='_c1'; 
name_str=['SAM_PSA_index',end_str];

%  filedir= 'C:\Users\Machine\matlab_storage_of_output_files\';
  savefilename =[filedir, name_str '.mat'];
 
 
 save_nr=1; % (1/0)

if save_nr==1
    % save(savefilename,'MA_save'); 
    save(savefilename,'ind_c1','ind_c2');
end


   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% letter
TextBox = uicontrol('style','text');
set(TextBox,'String',letter,'position',letter_position,'FontSize',32 ); % ,'FontWeight','bold'
                 set(TextBox,'foregroundcolor', [0 0 0], ...
            'backgroundcolor', [1 1 1]);
%%%%%%%%%%%%
% legend  
hs1='SAM';
hs2='PSA1';
l2=legend([h1 h2], hs1 , hs2);  
%set(l2,'FontSize',20, 'FontWeight', 'bold', 'location', 'NorthEastOutside'); 
 set(l2,'FontSize',14, 'FontWeight', 'bold', 'location', 'NorthWest'); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x y width height] for the figure, increase margins
    pos1 = get(gca, 'Position');
    pos1(1) = 0.2;
    pos1(3) = 0.6;
    set(gca, 'Position', pos1) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change to correct minus symbol
ax=gca;
set(ax,'ticklabelinterpreter','none')  %or 'tex' but not 'latex'
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
xticklabels(ax, strrep(xticklabels(ax),'-','–'));    
    
    
%%%%%%%%%%%%%%%%%%%%
% save figure

%  filedir ='C:\Users\Machine\matlab_storage_of_output_files\figures\';  
 filename=strcat(filedir,'figures\','SAM_PSA1_in-phase_years_',eof_str,'c');
 savefilename_c=filename;

% save as png 
orient landscape
% increase quality of picture by increase to -r500
% if figure_format==1
cd('G:\My Drive\ISO_CFA\matlab')
export_fig('-pdf','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c); % EPS works  func needs ghostscript 
% elseif figure_format==2
export_fig('-png','-painters', '-depsc','-nocrop','-opengl', '-r190',savefilename_c); % PNG  110
cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab')       

%% loop ERA-interim composite
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% era_name_nr=1; % 1 or 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Z500 *
% 2. 2mT *

if era_name_nr==1
name='z500';

elseif era_name_nr==2
name='2mT'; 

end

yr_s=1979;
yr_e=2011;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sea_nr=3;                    %%%%%%%%%%          Season
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



%       Size:       480x241x421
%        Dimensions: longitude,latitude,time

 addpath C:\Users\Machine\matlab_lib\Data\ERA_interim
% ncdisp('ERA_int_monthly_z500.nc') 
% ncdisp('ERA_int_monthly_2m_T.nc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. 0 Read in ERA-interim data

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
    % Monthly means of dail means
    % surface 2mT 
    name_c='ERA_int_monthly_2m_T_c.nc';
    era_T=ncread(name_c,'t2m'); % 2m temp
    era_T=era_T- 273.15;
    letter='c';
    lat1=-90; 
    lat2=-30;  
 
 
 elseif strcmp(name,'z500')==1 % && strcmp(season,'annual')==1
    % ECMWF ERA-Inerim (Dee et al. 2011)
    % Monthly means of dail means
    % geopotentail 500 hPa  
     
     
    name_c='ERA_int_monthly_z500_2.nc'; 
    era_z500=ncread(name_c,'z'); 
    era_z500=era_z500/9.80665;
    %   letter='';
    letter='b'; 
  
        if strcmp( proj,'stereo')==1
        lat1=-90;
        lat2=-10;
        end
 end
%
% era_count=396; % that overlaps with RICE record 1979-2011
era_count=find(era_year_num==yr_e)+11;  % data in months
%era_start=find(era_year_num==1980-1);

% annual 
era_start=find(era_year_num==yr_s);
div=ceil((era_count-era_start)/12);
    if strcmp(name,'z500')==1
        M=era_z500;
    elseif strcmp(name,'2mT')==1
        M=era_T;
    end
%%%%
%% resize grid
 tic
 reanalysis_nr=1;
 if reanalysis_nr==1 % Only for ERA-Interim and ERA-20C(higher resolution)
 
 [a b c]=size(M);
 
 
 %scale=0.5;
 scale=0.75; % correspond to hadisst
 %scale=1;
 
    for i=1:c
        
    era_z500_v(:,:,i)=resizem(M(:,:,i),scale);
    
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
    
    M=era_z500_v;
    era_lat=era_lat_v;
    era_long=era_lon_v;
 
 downscale_str=['downscale_',num2str(scale)];   
    
 end
 toc
 %%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% ===================================================
 M=M(:,(91:end),:);% SH
 era_lat_c=era_lat;

 %[nlon, nlat, ntim]= size(M);
 [nlon, nlat, ntim]= size(M); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  era_z500_c = reshape(era_z500, ntim, nlat, nlon);  %  z500_eof

era_z500_c = permute(M, [3 2 1]);  %  z500_eof

  
 %%% -=(1)=-
% Remove monthly climatology for all cells 

 era_z500_c2 = double(reshape(era_z500_c, ntim, nlat*nlon));  % One time series (column) for each grid point
 
 %  lat 1, 2, 3, 4,..........33...lat 1, 2, 3, 4.....33       until 33X240=7920
 
 
  [era_z500_c3,clim_z500] = annave(era_z500_c2);   % checked, Removes
%  seasonal cycle


 era_z500_c4= reshape(era_z500_c3,  ntim, nlat, nlon );

 
%%% -=(2)=- 
%  z = cosweight(z, lat);
% weighted to account for the decrease in area towards the pole.

%era_z500_c2=cosweight(era_z500_c,era_lat_c); % Original UoW function    (time x lat x lon) 
era_z500_c5=cosweight(era_z500_c4,era_lat_c(91:end));

% Back to old format again
era_z500_c6= reshape(era_z500_c5,  ntim, nlat, nlon );

M = permute(era_z500_c6, [3 2 1]);  %  z500_eof
 


% S_lon=266;S_lat=227;%check
% M ( S_lon, S_lat-20,2) 
%  M=cosweight_c(M,era_lat); % UoW function
%  M ( S_lon, S_lat-20,2) % check should changes value before and after
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Annual seasonal means

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
          ERA_M_annual_dummy=ERA_M_annual(1:end); % change to 1 to get 1979 annual  %%%%%%%%%%%%%%%%%%
          ERA_M_annual_c1(i,j,:)=ERA_M_annual_dummy;            
    elseif  strcmp(season,'MAMJJASON')==1
        
          ERA_M_season_mamjjason=nanmean(dummy1(3:11,:),1);  
          ERA_M_season_mamjjason=ERA_M_season_mamjjason(1:end);
          ERA_M_season_reg_mamjjason(i,j,:)= ERA_M_season_mamjjason; % for regression
          ERA_M_season_dummy_mamjjason=ERA_M_season_mamjjason; 
          ERA_M_season_mamjjason_c1(i,j,:)=ERA_M_season_dummy_mamjjason; 
          
   elseif  strcmp(season,'AMJJASON')==1
          ERA_M_season_amjjason=nanmean(dummy1(4:11,:),1);  
          ERA_M_season_amjjason=ERA_M_season_amjjason(1:end);
          ERA_M_season_reg_amjjason(i,j,:)= ERA_M_season_amjjason; % for regression
          ERA_M_season_dummy_amjjason=ERA_M_season_amjjason; 
          ERA_M_season_amjjason_c1(i,j,:)=ERA_M_season_dummy_amjjason; 
            
    end  
    
    end
end

 if  strcmp(season,'annual')==1
    Y=ERA_M_annual_c1;
 elseif  strcmp(season,'MAMJJASON')==1
    Y=ERA_M_season_mamjjason_c1;
 elseif  strcmp(season,'AMJJASON')==1
    Y=ERA_M_season_amjjason_c1;    
    
 end
%%%%%%%%%
% composite difference
 Xc=Y(:,:,ind_c1); % select using index SAM+ PSA-
 Xc2=nanmean( Xc,3);% composite surface
 
 
 Xcc=Y(:,:,ind_c2); % select using index SAM- PSA+
 Xcc2=nanmean( Xcc,3);% composite surface
 
 % minus
 Xccc=Xcc2-Xc2;
 
 
 clear Xcc2 Xc2
 toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %% t-test significance level, composites
 tic
 t_test_P_store=NaN(m,n);
 
 for i= 1:m
    for j= 1:n
%  if eof_type==2
%   x_wc=squeeze(Xc(i,j,(1:4)));
%  else
  x_wc=squeeze(Xc(i,j,:));
%  end
 y_wc=squeeze(Xcc(i,j,:));
 
 %[H t_test_P]=ttest(x_wc, y_wc);
 [H t_test_P]=ttest2(x_wc, y_wc);
 t_test_P_store(i,j)=t_test_P;
 
    end
 end
 
 clear  i j x_wc y_wc t_test_P
 toc
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Wilks test 
 tic 
 contour_alt_nr=2; % (1) t-test p (2) and Wilks test
if  contour_alt_nr==2 
 
% Wilks (2006, 2016) test
p_in=t_test_P_store;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
   % fdr_thres = 2*alp_c*i/no_grid_points;
%     fdr_thres = 2*alp_c*((no_grid_points+1)-rej_null(rej_null_ind(i)))/no_grid_points; % modified otherwise the sort/rank isn't taken into acount
   %     fdr_thres = 2*alp_c*((no_rej_null+1)-rej_null_ind(i))/no_grid_points; % modified otherwise the sort/rank isn't taken into acount
      fdr_thres = alp_c*((no_rej_null+1)-rej_null_ind(i))/no_grid_points; % remove 2*     eq. 2 Wilks 2016
 
    
    
    %find local p values that are less than fdr threshold
    if(pi_c<fdr_thres)
        p_fdr_rej(rn_x(i),rn_y(i)) = pi_c;
        n_fdr = n_fdr+1;
    end 
end

end
 
 toc
 %% Figure ERA-I
 %%%%%%%
  proj='stereo'; 

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

if era_name_nr==1
  c_limit=60;
  color_alt=9;
elseif  era_name_nr==2
  c_limit=2.0;
  color_alt=1;
end

hSurf=surfm(double(era_lat(91:end)),double(era_long),Xccc');
colormap(b2r(-c_limit,c_limit));

      

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
elseif color_alt==8    
    colormap(flipud(cbrewer2('div','RdYlGn',20)));    
elseif color_alt==9    
    colormap(flipud(cbrewer2('div','PuOr',20))); 
end

%%%%%%%%%%%%%%%%%%%%%%
% hSurfAx=(gca);
cRange= caxis;


 h=colorbar('EastOutside');
        pos=get(h,'Position');
  if era_name_nr==1
 

        pos(1)=pos(1)+0.0; % 0.035
        pos(2)=pos(2)+0.09;  
        pos(3)=pos(3)+ 0.012;  % widthcolorbar
        pos(4)=pos(4)-0.18;  % height colorbar
  elseif era_name_nr==2
        pos(1)=pos(1)-0.035; % 0.035
        pos(2)=pos(2)+0.09;  
        pos(3)=pos(3)+ 0.012;  % widthcolorbar
        pos(4)=pos(4)-0.16;  % height colorbar
        
  end
        
  set(h,'Position',pos)
  %%%%%%%%%%%%%%%
  h.TickLabels = strrep(h.TickLabels, '-', '–');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 
   if era_name_nr==1
    a1=axestext_c(1.05,+0.58, ['Z500 (m)'] ,'rotation',-90,'FontSize',22, 'FontWeight', 'bold');
     mlabel('MLabelParallel',-16  ,'PLabelLocation', [-75 -60 -45 -30],'PLabelMeridian',100)
       letter_position=   [200 760 60 80];
  elseif era_name_nr==2
    a1=axestext_c(1.035,+0.57, ['SAT (°C)'] ,'rotation',-90,'FontSize',22, 'FontWeight', 'bold'); 
     mlabel('MLabelParallel',-36  ,'PLabelLocation', [-75 -60 -45 -30],'PLabelMeridian',100)
       letter_position=   [200 740 60 80];
   end             
  
   
   
  
    
if contour_alt_nr==1
      h1= contourm(double(era_lat),double(era_long),t_test_P_store', [0.05],'-k','ShowText','off',...
   'LineWidth', 2);
elseif contour_alt_nr==2 % Wilks test

      h1= contourm(double(era_lat(91:end)),double(era_long),  p_fdr_rej', [0.05],'-k','ShowText','off',...    % wilks test
   'LineWidth', 2);  
end
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%
% letter

%    letter='a';

TextBox = uicontrol('style','text');

    letter_position=   [200 760 60 80];
   
     set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',40 ); 
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
%         if strcmp(name,'SIC')==1
            color_code=[.4 .4 .4];    
            plot3m(lat_cr,lon_cr,'-','LineWidth', 2,'color',color_code);
            bedmap2('gl','LineWidth', 2,'color',color_code);
%         end

    bedmap2('gl','LineWidth', 2, 'color',color_code );
    bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code)
end      
        
   hp1=  title(['(SAM−, PSA1+) minus (SAM+, PSA1−)',...
                 ],'FontSize',16, 'FontWeight', 'bold');
                                              
            pos=get(hp1,'Position');
           % pos(2)=pos(2)-0.03;
            pos(2)=pos(2)-0.3;
            set(hp1,'Position',pos,'FontSize',28, 'FontWeight', 'bold');
            
     set(h, 'FontSize',20, 'FontWeight', 'bold'); 
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

hold off
%%%%%%%%%%%%%%%%%%%%
% save figure
% filedir ='C:\Users\Machine\matlab_storage_of_output_files\figures\'; 
%     if rcontour_psa==1
%                 filename=strcat(filedir,'ndsic_',name,'_',iso,'_',season,'_',proj,'_',num2str(yr_s),'_',num2str(yr_e),'_wPSA_pattern');
    
                filename=strcat(filedir,'figures\','compositie_in_phase_SAM_PSA1_',name,...
                    '_',eof_str,'_',proj);
                


savefilename_c=filename;

 % save as png 
orient landscape

cd('G:\My Drive\ISO_CFA\matlab')
export_fig('-pdf','-painters', '-depsc','-opengl', '-r190',savefilename_c); % func needs ghostscript
export_fig('-png','-painters', '-depsc','-opengl', '-r190',savefilename_c); % PNG  ,'-nocrop'
cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab')  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Sea Ice HadISST
tic
clearvars -except ind_c1 ind_c2 eof_type contour_alt_nr filedir

%%%%%%%%%%%%%%%%%%%%%%
%  clear all
 close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input

param_nr=2;% (1-SST/ 2 SIC)

% SST_dataset=1; % (1-HadISST, 2-ERSST) 

rcontour_psa=0; %(0)(1) show where PSA pattern overlap interferance (2) shows where PSA2 is significant

%coast_nr=0;  % On/Off
    
if param_nr==2
    name='SIC';
end
 
 site='RICE';
  yr_s=1979;
  %yr_e=2009; 
  yr_e=2011; %%%%%%%%% RICE 2009 %%% suppl. material
%%%%%%%%%%%%%%%%%%%%%
sea_nr=1;
season='annual'; 
%%%%%%%%%%%%%%%%%%%%%%


  
if  strcmp(name,'SIC')==1
  proj='stereo';
  lat1=-90;
  lat2=-50;
  box_use=0;  % (0/1)  
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HadISST data
% HadISST data
% version 1.1 Rayner et al. 2003
%
%           standard_name = 'sea_ice_area_fraction'
%           long_name     = 'Monthly 1 degree resolution sea ice concentration'


addpath  C:\Users\Machine\matlab_lib\Data\HadISST_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    name_c='HadISST_ice_c.nc';
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
% % cosweight
% S_lon=17;S_lat=170;%check
% M ( S_lon, S_lat-5,2) 
%  M=cosweight_c(M,HadISST_lat); % UoW function
%  M ( S_lon, S_lat-5,2)

HadISST_year_num(HadISST_start);
HadISST_year_num(HadISST_count);

% HadISST_year_num()-HadISST_year_num()
% HadISST_start:HadISST_count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ===================================================
% seasonality 
 
 HadISST_lat_c=HadISST_lat;
 [nlon, nlat, ntim]= size(M);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% -=(1)=-
%  z = cosweight(z, lat);
% weighted to account for the decrease in area towards the pole.

%  era_z500_c = reshape(era_z500, ntim, nlat, nlon);  %  z500_eof

HadISST_z500_c = permute(M, [3 2 1]);  %  z500_eof


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 %%% -=(1)=-
% Remove monthly climatology for all cells 

 HadISST_z500_c2 = double(reshape(HadISST_z500_c, ntim, nlat*nlon));  % One time series (column) for each grid point
 
 %  lat 1, 2, 3, 4,..........33...lat 1, 2, 3, 4.....33    
 
% Removes seasonal cycle
 [HadISST_z500_c3,clim_z500] = annave(HadISST_z500_c2);   % checked
 
 
 HadISST_z500_c4= reshape(HadISST_z500_c3,  ntim, nlat, nlon );
 
  %%% -=(2)=-
 HadISST_z500_c5=cosweight(HadISST_z500_c4,HadISST_lat_c); % Original UoW function    (time x lat x lon) 

% Back to old format again
HadISST_z500_c6= reshape(HadISST_z500_c5,  ntim, nlat, nlon );

M = permute(HadISST_z500_c6, [3 2 1]);  %  z500_eof
 
 %==========================================================
% Annual seasonal means
% run from start if year is changed above
% (360x180x1728) (Long,Lat,month)
[m n t]=size(M(:,:,HadISST_start:HadISST_count));% time(months) % 1894-2011
for i= 1:m
    for j= 1:n

        dummy1=reshape(squeeze(M(i,j,HadISST_start:HadISST_count)), 12, div);

   
if strcmp(season,'annual')==1
         HadISST_M_annual_temp=nanmean(dummy1,1); % annual HadISST-values

         
%          HadISST_M_annual_temp=anomaly(HadISST_M_annual_temp); % try
%          normalizing
        % HadISST_M_annual_dummy=detrend(HadISST_M_annual_temp(1:end));
          HadISST_M_annual_dummy=HadISST_M_annual_temp(1:end);  
          
         HadISST_M_annual(i,j,:)=HadISST_M_annual_temp(1:end);
         HadISST_M_annual_c1(i,j,:)= HadISST_M_annual_dummy; % detrended
         
elseif strcmp(season,'MAMJJASON')==1
    
          HadISST_M_season_mamjjason=nanmean(dummy1(3:11,:),1);  %3. seasonal (JJA)
          HadISST_M_season_mamjjason=HadISST_M_season_mamjjason(1:end);
          HadISST_M_season_reg_mamjjason(i,j,:)=HadISST_M_season_mamjjason; % means for regression         
          HadISST_M_season_dummy_mamjjason=HadISST_M_season_mamjjason; 
          HadISST_M_season_mamjjason_c1(i,j,:)=HadISST_M_season_dummy_mamjjason;
          
elseif strcmp(season,'AMJJASON')==1
    
          HadISST_M_season_amjjason=nanmean(dummy1(3:11,:),1);  %3. seasonal (JJA)
          HadISST_M_season_amjjason=HadISST_M_season_amjjason(1:end);
          HadISST_M_season_reg_amjjason(i,j,:)=HadISST_M_season_amjjason; % means for regression         
          HadISST_M_season_dummy_amjjason=HadISST_M_season_amjjason; 
          HadISST_M_season_amjjason_c1(i,j,:)=HadISST_M_season_dummy_amjjason;          
  
end

%   clear       dummy1 HadISST_M_annual_dummy HadISST_M_annual_temp

    end
end

if strcmp(season,'annual')==1
    Y=HadISST_M_annual_c1;
elseif strcmp(season,'MAMJJASON')==1    
    Y=HadISST_M_season_mamjjason_c1;
elseif strcmp(season,'AMJJASON')==1    
    Y=HadISST_M_season_amjjason_c1;
end



% composite difference
 Xc=Y(:,:,ind_c1); % select using index SAM+ PSA-
 Xc2=nanmean( Xc,3);% composite surface
 

 Xcc=Y(:,:,ind_c2); % select using index SAM- PSA+
 Xcc2=nanmean(Xcc,3);% composite surface
 
 % minus
 Xccc=Xcc2-Xc2;

toc

 %% t-test
 tic
 
 t_test_P_store=NaN(m,n);
 
 for i= 1:m
    for j= 1:n
%  if eof_type==1
 x_wc=squeeze(Xc(i,j,:));
%  else
%  x_wc=squeeze(Xc(i,j,(1:4)));    
% %  end
 y_wc=squeeze(Xcc(i,j,:));
 
 [H t_test_P]=ttest2(x_wc, y_wc);
 
 t_test_P_store(i,j)=t_test_P;
 
    end
 end
 
 clear  i j x_wc y_wc t_test_P
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if contour_alt_nr==2 % Wilks test
 
p_in=t_test_P_store;
%p_in=p_in(:,(91:180)); % SH
p_in=p_in(:,(146:167)); % SH sea ice edge

% A_nan=isnan(p_in);
% A_nan_nr=sum(sum(A_nan));


% find active sea ice cells
% %%%%%%%%%%%%%%%%%%%%%
% load('HadISST_M_annual_detrend'); % go w. 0.1
% wr_c=nanmean(HadISST_M_annual_detrend,3);
%  
% wr_c=wr_c(:,(146:167));
% wr_c2=wr_c; 
% wr_c2(isnan(wr_c))=0; % turn nans to zero
%  
% A_ind2=find(wr_c2~=0);
% A_nan_nr2=length(A_ind2);
% %%%%%%%%%


Xccc_c=Xccc;
Xccc_c(isnan(Xccc))=0; % turn nans to zero

A_ind=find(Xccc_c(:,(146:167))~=0);

A_ind_nr=length(A_ind);

%A_nan_nr2=180*360-A_nan_nr;
%A_nan_nr2=90*360-A_nan_nr;
A_nan_nr2=22*360-A_ind_nr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%      fdr_thres = 2*alp_c*i/no_grid_points;
%     fdr_thres = 2*alp_c*((no_rej_null+1)-rej_null_ind(i))/no_grid_points; % modified otherwise the sort/rank isn't taken into acount
     fdr_thres = alp_c*((no_rej_null+1)-rej_null_ind(i))/no_grid_points;% remove 2*
    
    %find local p values that are less than fdr threshold
    if(pi_c<fdr_thres)
        p_fdr_rej(rn_x(i),rn_y(i)) = pi_c;
        n_fdr = n_fdr+1;
    end 
end

% p_fdr_rej=fdr_wilks_test(p_in); % Wilks test made for 360x180 grid
 end
 toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig SIC

  proj='stereo';
  lat1=-90;
  lat2=-50;
  
 % c_limit=0.2;
   c_limit=0.1501;

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

        mlabel('MLabelParallel',-54  ,'PLabelLocation', [-75 -60 -45 -30],'PLabelMeridian',100) 
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hSurf=surfm(double(HadISST_lat),double(HadISST_lon),Xccc');
colormap(b2r(-c_limit,c_limit));

  color_alt=6;    

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

  h=colorbar('EastOutside');  
%   set(h, 'Position', [.76 .145 .015 .75], 'FontSize',18, 'FontWeight', 'bold');  
     % modifiy colorbar
    pos=get(h,'Position');
    
 a1=axestext_c(1.07,+0.525, ['SIC'] ,'rotation',-90,'FontSize',22, 'FontWeight', 'bold');         
          
% %         pos(1)=pos(1)+0.024;      %%%%%%%%%%%%%%%%%%%% x-position of colorbar   (+) --->>   +0.02;
        pos(1)=pos(1)+0.024; % 0.035
        pos(2)=pos(2)+0.09;  
        pos(3)=pos(3)+ 0.012;  % widthcolorbar
        pos(4)=pos(4)-0.18;  % height colorbar         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
     set(h,'Position',pos)
     set(h, 'FontSize',20, 'FontWeight', 'bold'); 
     
     h.TickLabels = strrep(h.TickLabels, '-', '–');
%%%%%%%%%%%%%

if contour_alt_nr==1
      h1= contourm(double(HadISST_lat),double(HadISST_lon),t_test_P_store', [0.05],'-k','ShowText','off',...
   'LineWidth', 2,'Color',[.99 .0 .4]);
elseif contour_alt_nr==2
%       h1= contourm(double(HadISST_lat),double(HadISST_lon),p_fdr_rej', [0.05],'-k','ShowText','off',...
%    'LineWidth', 2);
      h1= contourm(double(HadISST_lat(146:167)),double(HadISST_lon),p_fdr_rej', [0.05],'-k','ShowText','off',...
   'LineWidth', 2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% letter
letter='d';
TextBox = uicontrol('style','text');

   letter_position=  [180 800 80 80];%
   set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',40 ); % x position y position xsize ysize

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

%       addpath C:\Users\benman\matlab_mapping
      
        if strcmp(name,'SIC')==1
            color_code=[.4 .4 .4];    
            plot3m(lat_cr,lon_cr,'-','LineWidth', 2,'color',color_code);
            bedmap2('gl','LineWidth', 2,'color',color_code);
        end

    
    bedmap2('gl','LineWidth', 2, 'color',color_code );
    bedmap2('patchshelves','LineWidth', 1.0,'facecolor','none','frame','on','edgecolor',color_code)
end 
  
   hp1=  title(['(SAM−,PSA+) minus (SAM+,PSA−) '],'FontSize',28, 'FontWeight', 'bold');    

            pos=get(hp1,'Position');
            pos(2)=pos(2)-0.05;
            set(hp1,'Position',pos,'FontSize',28, 'FontWeight', 'bold'); 
 
%%%%%%%%%%%%%%%%%%%%
% save figure
eof_str='varimax';
%  filedir ='C:\Users\benman\matlab_storage_of_output_files\figures\'; 
filename=strcat(filedir,'figures\','HadISST_compositie_in_phase_SAM_PSA1_',name,'_',eof_str,'_',season,'_',proj);
savefilename_c=filename;

 % save as png 
orient landscape
cd('G:\My Drive\ISO_CFA\matlab')
  export_fig('-pdf','-painters', '-depsc','-opengl', '-r190',savefilename_c); %func needs ghostscript 

  export_fig('-png','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG  ,'-nocrop'
cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          