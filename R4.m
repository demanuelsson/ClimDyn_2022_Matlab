

    load('SAM_NOAA_2.mat')
    title_str1='SAM NOAA';
    %SAM_x=nanmean(SAM_NOAA_2,2);% annual means
    SAM_x=nanmean(SAM_NOAA_2(:,(4:11)),2);% AMJJASON means
    
    
    load('C:\Users\Machine\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2011_AMJJASON_c52_annual_mean_varimax_c2.mat');  % AMJJASON
    title_str2='SAM';

    
    f1=fig('units','inches','width',11.5,'height',7.5,'font','Helvetica','fontsize',18,'border','on');
    box on
    
    figure
    plot(MA_PCs_save(:,1),MA_PCs_save(:,2),'-*')
    hold on
    yrs=[1979:2015]';
    plot(yrs,SAM_x,'.-r')
    
    legend('SAM this study', 'SAM NOAA')
    
    