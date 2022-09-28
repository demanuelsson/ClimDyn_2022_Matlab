%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RICE_d18O_anomaly_DE21.m
% d18O time series
% Daniel Emanuelsson 2021
% MATLAB 2020b
%%%%%%%%%%%%%%%%%%%%%%%%%%
% uses imput from composite_SAM_PSA1_in_phase_ERA_I_z500_2mT_HadISST_SIC_c15.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% close all


        iso_alt_nr=2;
        
        if iso_alt_nr==1
             load('C:\Users\Machine\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c23.mat');
        elseif iso_alt_nr==2
           load('C:\Users\Machine\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_AMJJASON_c24.mat'); % AMJJASON 
%             load('C:\Users\benman\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_clim_c24.mat'); % clim removed, no difference             
        end

f1=fig('units','inches','width',18,'height',5,'font','Helvetica','fontsize',18,'border','on'); 

    
      hold on
      box on
      
      
%plot(MA_save(:,1),anomaly(nanmoving_average(MA_save(:,3),10)),'-k','LineWidth',2.5)
plot(MA_save((1980:end),1),MA_save((1980:end),3),'-k','LineWidth',2.5)

plot(MA_save((1980:end),1),MA_save((1980:end),3),'.k','MarkerSize',22)

wr_c=nanmean(MA_save((1980:end),3));
lrc3=hline(wr_c,':');
  
set(lrc3,'Color',[0.1,0.1,0.1],'LineWidth',1.5);
  

wr_c1=nanmean(MA_save((1980:end),3))+nanstd(MA_save((1980:end),3)); % std lines
lrc4=hline(wr_c1,':');
set(lrc4,'Color',[0.9,0.0,0.0],'LineWidth',1.5);  


wr_c1=nanmean(MA_save((1980:end),3))-nanstd(MA_save((1980:end),3)); % std lines
lrc4=hline(wr_c1,':');
set(lrc4,'Color',[0.9,0.0,0.0],'LineWidth',1.5); 


yr_s=1979;
yr_e=2011;

set(gca,'XLim',[1978 2012])
set(gca,'YLim',[-27 -19])
%set(gca,'XLim',[1900 2011])

t_c=[1979:2011]';

% from composite_SAM... file
load('C:\Users\Machine\matlab_storage_of_output_files\SAM_PSA_index_c1.mat')

plot(t_c(ind_c1,1),wr_c,'*k','MarkerSize',14,'linewidth',1.2 ) % years

plot(t_c(ind_c2,1),wr_c,'ok','MarkerSize',14,'linewidth',1.2 ) % years

 
              set(gca,...
          'linewidth',2,...
          'FontWeight','bold' );

     xlabel('Year','FontWeight','bold','FontSize',20 );
     ylabel(['{\delta}^1^8O (',char(8240),')'],'FontWeight','bold','FontSize',20 );
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x y width height] for the figure, increase margins
    pos1 = get(gca, 'Position');
    pos1(1) = 0.2;
    pos1(3) = 0.6;
    set(gca, 'Position', pos1)
    
%%%%%%%%%%%%
 hs1=['{\delta}^1^8O']; %(',char(8240),')'];   

l2=legend(['RI ', hs1]);
 set(l2,'FontSize',14, 'FontWeight', 'bold', 'location', 'NorthWest');    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% letter
letter='a';

TextBox = uicontrol('style','text');
letter_position=   [200 380 60 60];
   
set(TextBox,'String',letter,'position',letter_position,'FontSize',32 ); % ,'FontWeight','bold'
                 set(TextBox,'foregroundcolor', [0 0 0], ...
            'backgroundcolor', [1 1 1]);
%%%%%%%%%%%%%%%%        

% Change to correct minus symbol
ax=gca;
set(ax,'ticklabelinterpreter','none')  %or 'tex' but not 'latex'
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
xticklabels(ax, strrep(xticklabels(ax),'-','–'));        
        
%%%%%%%%%%%%         
save_fig=1;

if save_fig==1
iso_str='d18O';

filename=['RI_anomaly_',iso_str,'_',num2str(yr_s),'-',num2str(yr_e),'_fig'];


filedir ='C:\Users\Machine\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);
cd('G:\My Drive\ISO_CFA\matlab')
 % save as png 
 orient landscape
 export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c);            
 export_fig('-pdf','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c);          
cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab')  
 
end            