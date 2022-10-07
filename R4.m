
% Comparison to NOAA SAM index
% https://stateoftheocean.osmc.noaa.gov/atm/sam.php


% ncdisp('C:\Users\Machine\matlab_lib\Data\samNOAA.nc')% Accessed on 4/10/22

SAM_x2=ncread('C:\Users\Machine\matlab_lib\Data\samNOAA.nc','SAM');
SAM_time=ncread('C:\Users\Machine\matlab_lib\Data\samNOAA.nc','TSAXIS');
SAM_time_c=SAM_time+datenum(1900,1,1);  
 
 date_vec=datevec(SAM_time_c); 
 yyyy=date_vec(:,1); 
 mm=date_vec(:,2);
 SAM_year_num=yyyy+(mm/12-1/12);
 
%  mean(SAM_x2(12+4:12+11));

 month_ind=repmat(1:12,1,43)';
 ind_c=find(month_ind>=4 & month_ind<=11);
 
% SAM_NOAA_eW=nan(43,1);
 SAM_NOAA_eW=[];
 
 for i=1:8:length(ind_c)
     
 wr_c=SAM_x2(ind_c(i:i+7));
 SAM_c2=nanmean(wr_c);
 SAM_NOAA_eW=[SAM_NOAA_eW;SAM_c2];
  end
 



    load('SAM_NOAA_2.mat') % old file
    title_str1='SAM NOAA';
    %SAM_x=nanmean(SAM_NOAA_2,2);% annual means
    SAM_x=nanmean(SAM_NOAA_2(:,(4:11)),2);% AMJJASON means
    
    
    load('C:\Users\Machine\matlab_storage_of_output_files\ERA-Interim_PCs_Z500_lim0-360_-20_-90_1979-2011_AMJJASON_varimax.mat');  % AMJJASON
    title_str2='SAM';

    
    f1=fig('units','inches','width',11.5,'height',6.5,'font','Helvetica','fontsize',18,'border','on');
    box on
    
    plot(MA_PCs_save(:,1),MA_PCs_save(:,2),'-*k','linewidth',3)
    hold on
    yrs=[1979:2015]';
    %plot(yrs,SAM_x,'*-r','linewidth',3)
    
  %  plot([1979:2021]',SAM_NOAA_eW,'*-r','linewidth',3)
     plot([1979:2011]',SAM_NOAA_eW(1:33),'*-r','linewidth',3)
    legend('SAM this study', 'SAM NOAA','Location','northwest')
    
    
                     set(gca,...
          'linewidth',3,...
          'FontWeight','bold' ); 
      
  xlabel('Year','FontWeight','bold','FontSize',18 );   
  ylabel('SAM','FontWeight','bold','FontSize',18 ); 
  
  x1=SAM_NOAA_eW(1:33);
  y1=MA_PCs_save(:,2);
  
[r,p, nv]=corrcoef_df_c2(detrend(x1),detrend(y1));
rp=[r(2),p(2)]
p_level(p(2))
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change to correct - symbol
ax=gca;
set(ax,'ticklabelinterpreter','none')  %or 'tex' but not 'latex'
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
      
set(gca,'XLim',[1978 2012])    
    
    save_fig=1;
if save_fig==1
    site='RI';
    filename=['SAM_indices_comparison'];
    filedir ='C:\Users\Machine\matlab_storage_of_output_files\figures\';
    savefilename_c=strcat(filedir,filename);

    % save as png 
    orient landscape
    cd('G:\My Drive\ISO_CFA\matlab')
    export_fig('-png','-painters', '-depsc','-opengl', '-r190',savefilename_c);  % ,'-nocrop'         
    export_fig('-pdf','-painters', '-depsc','-opengl', '-r190',savefilename_c);  % ,'-nocrop'
    cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab')  
end 
    
    
    
    
    