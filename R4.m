

    load('SAM_NOAA_2.mat')
    title_str1='SAM NOAA';
    %SAM_x=nanmean(SAM_NOAA_2,2);% annual means
    SAM_x=nanmean(SAM_NOAA_2(:,(4:11)),2);% AMJJASON means
    
    
    load('C:\Users\Machine\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2011_AMJJASON_c52_annual_mean_varimax_c2.mat');  % AMJJASON
    title_str2='SAM';

    
    f1=fig('units','inches','width',11.5,'height',6.5,'font','Helvetica','fontsize',18,'border','on');
    box on
    
    plot(MA_PCs_save(:,1),MA_PCs_save(:,2),'-*k','linewidth',3)
    hold on
    yrs=[1979:2015]';
    plot(yrs,SAM_x,'*-r','linewidth',3)
    
    legend('SAM this study', 'SAM NOAA','Location','northwest')
    
    
                     set(gca,...
          'linewidth',3,...
          'FontWeight','bold' ); 
      
  xlabel('Year','FontWeight','bold','FontSize',18 );   
  ylabel('SAM','FontWeight','bold','FontSize',18 );  
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change to correct - symbol
ax=gca;
set(ax,'ticklabelinterpreter','none')  %or 'tex' but not 'latex'
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
      
    
    
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
    
    
    
    
    