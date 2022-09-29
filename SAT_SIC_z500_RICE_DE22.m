%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAT_SIC_z500_RICE_DE21.m
% d18O corr w. SAT, SIC and Z500
% Daniel Emanuelsson, 2021
% MATLAB 2020b 
% Github version 2 
% Moving correlation, scatter plots and moving neff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
% * IPO_vline.m DE
% * fig.m fileexchange
% * export_fig.m fileexchange
%%%%%%%%%%%%%%%%%%%%%%
clear
close all
%%%%%
addpath C:\Users\Machine\matlab_lib\Data
addpath C:\Users\Machine\matlab_lib\Data\ERA_interim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        iso_alt_nr=2;       
        if iso_alt_nr==1
             load('C:\Users\Machine\matlab_storage_of_output_files\RICE_combined_Deep_1213B_c23.mat');
             iso_season_str='annual';
        elseif iso_alt_nr==2
            
              load('C:\Users\Machine\matlab_storage_of_output_files\RI_combined_Deep_1213B_annual_means_no_summer_AMJJASON_c24.mat'); % AMJJASON             
            iso_season_str='ExtendedWinter';  
        end

    stacked_record_annual_Ma=MA_save((1901:end),:);
    iso_str='d18O';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig   %%% moving correlation with RICE d18O
sliding_w=11;%%%%%%%%%%%%%%%%%%%%%%%%%%%
add_yr=2; % (0) 2009 (2) 2011
    

h=fig('units','inches','width',14,'height',4.5,'font','Helvetica','fontsize',16,'border','on');
   box on
   hold on   
   set(gca,'YLim',[-1 1])
 %%%%%%%%%%%%%%            
addpath C:\PHD\matlab_lib\Library
%   edit IPO_vline
%   IPO_vline(15)
   hold on 
   
load('C:\Users\Machine\matlab_storage_of_output_files\Ma_ADP_index_HadISST_SIC.mat'); % ADP index   ADP_index.m  
load('C:\Users\Machine\matlab_storage_of_output_files\SST_Nino-4_HadISST_1979_2012.mat');

      %%%%%%%%%%%%%%%%%%%%%%% RICE
le=length(stacked_record_annual_Ma(80:110+add_yr));

Ma_corr2=NaN(le,3); % use NaN instead of zero because NPGO index is shorter

% iso
Ma_corr2(:,1)=detrend(stacked_record_annual_Ma((80:110+add_yr),3));
Ma_corr2(:,2)=detrend(Ma_save_nino(1:31+add_yr));

Ma_corr2(:,3)=detrend(MA_ADP_save((8:38+add_yr),2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable 1: Nino-4 SST
% not using function to get better overview
% 

windowSize=sliding_w;

[N,M] = size(Ma_corr2);

correlationTS = nan(N, 1);
correlationTS_p= nan(N, 1);
correlationTS_rlo=nan(N, 1);
correlationTS_rup=nan(N, 1);


% correlationTS = nan(N, M-1);
% correlationTS_p= nan(N, M-1);
% correlationTS_rlo=nan(N, M-1);
% correlationTS_rup=nan(N, M-1);


indexColumn=1;
 %for t = windowSize+1:N % (Original)
for t = floor(windowSize/2)+1:(N-floor(windowSize/2)) % Modified
        % 6:105:
     
%      dataMatrix=Ma_corr4;
 %C = corrcoef(dataMatrix(t-windowSize:t, :)); %  (Original)
 % C = corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function 
 
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :)); % Change to centered moving correlation window
  %[C, p, rlo, rup ]= corrcoef(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :), 'alpha', alpha); 
  
  X=Ma_corr2(t-floor(windowSize/2):t+floor(windowSize/2),1);
  Y=Ma_corr2(t-floor(windowSize/2):t+floor(windowSize/2),2); 
  
  
  [C, p, rlo, rup ]= corrcoef_df(X,Y); % Change to centered moving correlation window 
 
 
%     idx = setdiff(1:M, [indexColumn]);
%     correlationTS(t, :) = C(indexColumn, idx);
%     correlationTS_p(t, :)= p(indexColumn, idx);
%     correlationTS_rlo(t, :)= rlo(indexColumn, idx);
%     correlationTS_rup(t, :)= rup(indexColumn, idx);
    
    
%    idx = setdiff(1:M, [indexColumn]);
    correlationTS(t, :) = C(1, 2);
    correlationTS_p(t, :)= p(1, 2);
    correlationTS_rlo(t, :)= rlo(1, 2);
    correlationTS_rup(t, :)= rup(1, 2);
     
    
end


clear X Y t C p rlo rup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable 2: ADP

correlationTS_c = nan(N, 1);
correlationTS_p_c= nan(N, 1);
correlationTS_rlo_c=nan(N, 1);
correlationTS_rup_c=nan(N, 1);


for t = floor(windowSize/2)+1:(N-floor(windowSize/2)) % Modified
        % 6:105:
     
%      dataMatrix=Ma_corr4;
 %C = corrcoef(dataMatrix(t-windowSize:t, :)); %  (Original)
 % C = corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function 
 
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :)); % Change to centered moving correlation window
  %[C, p, rlo, rup ]= corrcoef(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :), 'alpha', alpha); 
  
  X=Ma_corr2(t-floor(windowSize/2):t+floor(windowSize/2),1);
  Y=Ma_corr2(t-floor(windowSize/2):t+floor(windowSize/2),3); 
  
  
  [C, p, rlo, rup ]= corrcoef_df(X,Y); % Change to centered moving correlation window 
 
 
%     idx = setdiff(1:M, [indexColumn]);
%     correlationTS(t, :) = C(indexColumn, idx);
%     correlationTS_p(t, :)= p(indexColumn, idx);
%     correlationTS_rlo(t, :)= rlo(indexColumn, idx);
%     correlationTS_rup(t, :)= rup(indexColumn, idx);
    
    
%    idx = setdiff(1:M, [indexColumn]);
    correlationTS_c(t, :) = C(1, 2);
    correlationTS_p_c(t, :)= p(1, 2);
    correlationTS_rlo_c(t, :)= rlo(1, 2);
    correlationTS_rup_c(t, :)= rup(1, 2);
     
    
end


clear X Y t C p rlo rup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  [correlationTS2, correlationTS_p2, correlationTS_rlo2, correlationTS_rup2]=movingCorrelation_c(Ma_corr2,sliding_w,1, 0.05);    


% edit corrcoef
% edit corrcoef_df
% edit movingCorrelation_c


color_code=[0.3,0.3,0.3];

h70=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rup(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code);
h7=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS(:,1),'-','LineWidth',4,'MarkerSize',14,'Color',color_code);
h71=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rlo(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code);




h80=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rup_c(:,1)*(-1),'--','LineWidth',2,'MarkerSize',14,'Color',[0.6,0.0,0.0]); 
h8=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_c(:,1)*(-1),'-','LineWidth',4,'MarkerSize',14,'Color',[0.6,0.0,0.0]);
h81=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rlo_c(:,1)*(-1),'--','LineWidth',2,'MarkerSize',14,'Color',[0.6,0.0,0.0]); 

set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
    
   %%%%% Legend %%%%%%%%%%%%
    hs61='PC1 SAT';
    hs89='PC2 SAT';
    hs62='PC3 SAT';
    hs7='Nino-4 SST'; 
    hs8='ADP (−1)'; 
        
    hs7_c=['r ({\delta}^1^8O, ',hs7,')'];
    hs8_c=['r ({\delta}^1^8O, ',hs8,')'];

  h_leg=legend([ h7 h8 ],  hs7_c, hs8_c);
 %  h_leg=legend([  h8 ],   hs8_c);
  set(h_leg, 'location', 'SouthWest','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change to correct - symbol
ax=gca;
set(ax,'ticklabelinterpreter','none')  %or 'tex' but not 'latex'
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
  
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[80 350 60 60];
  letter='d';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
 %%%%%%%%%%%%%%%%%%%%%%%%%
        set(gca,'XLim',[1977 2012])
        %hl=hline([-0.5,0,0.5]);
        hl=hline([0]);
        set(hl,'LineWidth',3,'Color',[0.0,0.0,0.0]);
        
%         xlabel('Year','FontWeight','bold','FontSize',18 );
         ylabel('Correlation','FontWeight','bold','FontSize',20 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
% save fig.
end_str='_c2';
filename=['ADP_RICE_w_',num2str(sliding_w),'_', iso_str,'_',iso_season_str, end_str]; % poor chooise of name

filedir ='C:\Users\Machine\matlab_storage_of_output_files\figures\';

savefilename_c=strcat(filedir,filename);
% save as png 
orient landscape
% increase quality of picture by increase to -r300
cd('G:\My Drive\ISO_CFA\matlab')
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
export_fig('-pdf','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c); % PNG-nocrop'
cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab') 
%% scatter plots
close all


%alt=5; % (1) ADP (2) nino-4 (3) SAM (4) PSA1 (5) PSA2

for alt=1:4

 h=fig('units','inches','width',6,'height',6,'font','Helvetica','fontsize',16,'border','on');  
   box on
   hold on
   
  x1=stacked_record_annual_Ma((80:110+add_yr),3);
  
  if alt==1
    y1=MA_ADP_save((8:38+add_yr),2);
  elseif alt==2
    y1=Ma_save_nino(1:31+add_yr);
  elseif alt>2
      

    eof_alt_nr=1;
    if eof_alt_nr==1
        load('C:\Users\Machine\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2011_AMJJASON_c52_annual_mean_varimax_c2.mat') %Z500
    end
        
      if alt==3
        pc_nr=2;
        factr=1;
      elseif alt==4
        pc_nr=3;
        factr=1;
      elseif alt==5
        pc_nr=4;
        factr=-1;
      end
      y1=MA_PCs_save((1:33), pc_nr)*factr;
      
      if alt==2
          
      end
    
  end
   
%plot(x1,y1,'*','MarkerSize',10,'Color',[0.0,0.0,0.0])
plot(y1,x1,'*','MarkerSize',10,'Color',[0.0,0.0,0.0])

 if alt==1 || alt==2 %|| alt==3
  %  plot(x1(end:end),y1(end:end),'*','MarkerSize',10,'Color',[1,0.0,0.0])
    plot(y1(end:end),x1(end:end),'*','MarkerSize',10,'Color',[1,0.0,0.0])
 elseif alt==4 || alt==3 % year 2010 and 2011 as red *
  %  plot(x1(end-1:end),y1(end-1:end),'*','MarkerSize',10,'Color',[1,0.0,0.0]) 
    plot(y1(end-1:end),x1(end-1:end),'*','MarkerSize',10,'Color',[1,0.0,0.0]) 
 end


           set(gca,...
          'linewidth',3,...
          'FontWeight','bold' );
      
ylabel(['{\delta}^1^8O (',char(8240),')'])

x3=x1;
y3=y1;

if alt==1 || alt==2  
x1=x1(1:end-1);
y1=y1(1:end-1);

elseif alt==3 || alt==4 %|| alt==1 || alt==2  
    
x1=x1(1:end-2);
y1=y1(1:end-2);    

end

pol_1=polyfit(y1,x1,1);
%x2=[-29:0.5:-15];
x2=[-3:0.01:3];
polv_1=polyval(pol_1,x2);

plot(x2,polv_1,'-k','LineWidth',2)
set(gca,'YLim',[-28 -18])

%%%%%%%%%
show_slope_nr2=1; % present during polynya

if show_slope_nr2==1 && (alt==3 || alt==4 )
    
% pol_2=polyfit(y3(32:33),x3(32:33),1);    
% polv_2=polyval(pol_2,x2);
if alt==3
offset_y=3.5;
elseif alt==4
offset_y=4.25;    
end
plot(x2,offset_y+polv_1,'--r','LineWidth',0.5)    
end


%%%%%%%%%



if alt==1
    
    
  title_str='ADP SIC';
    
elseif alt==2
    
set(gca,'XLim',[-0.8 0.8])  

 title_str='Nino-4 SST (°C)';
    
elseif alt==3 
  title_str='SAM (m s.d.^−^1)';  
elseif alt==4 
  title_str='PSA1 (m s.d.^−^1)';
elseif alt==5 
  title_str='PSA2 (m s.d.^−^1)';  
 
end

xlabel(title_str)


%correlation
%  factor_nr
[r,p, nv]=corrcoef_df_c2(detrend(x1),detrend(y1));
rp=[r(2),p(2)]
p_level(p(2))
%nv
% edit corrcoef_df_c2
% help getdof

%  x=X;
%  y=Y;

 nc=size(x1);
   r1 = corrcoef(x1(2:end),x1(1:end-1),'rows','pairwise');
   r2 = corrcoef(y1(2:end),y1(1:end-1),'rows','pairwise');
%    r1 = corrcoef(x(2:end),x(1:end-1),'rows','pairwise');
%    r2 = corrcoef(y(2:end),y(1:end-1),'rows','pairwise');
%    %r1 = acf(x,1);
%    %r2 = acf(y,2);
%    
  %  nv = nv*(1-r1(2)*r2(2))/(1+r1(2)*r2(2))
  
    neff = round(nc(1)*(1-r1(2)*r2(2))/(1+r1(2)*r2(2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if alt==1
x_text=0.075;    
y_text=0.06;
letter='e';
elseif alt==2
x_text=0.075;    
y_text=0.95;

letter='f';
elseif alt==3
y_text=0.98;
x_text=0.9;
letter='b';
elseif alt==4
y_text=0.98;
x_text=0.1;
letter='c';
elseif alt==5
y_text=0.98;
x_text=0.1;
letter='d';
end

 r_str=num2str(round(rp(1)*100)/100);
 r_str_c=strrep(r_str, '-', '–');
 axestext_c(x_text,y_text,['r = ',r_str_c,', ',p_level(p(2)),', n_e_f_f = ',num2str(neff)],'FontWeight','bold','FontSize',14 );   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change to correct minus symbol
ax=gca;
set(ax,'ticklabelinterpreter','none')  %or 'tex' but not 'latex'
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[0 485 35 65];
%   letter='b';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',30 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%          
% save fig.         
filename=['scatter_plot_',num2str(alt),'_', iso_str,'_',iso_season_str,'_c1']; % poor chooise of name
filedir ='C:\Users\Machine\matlab_storage_of_output_files\figures\';

savefilename_c=strcat(filedir,filename);
% save as png 
orient landscape
cd('G:\My Drive\ISO_CFA\matlab')
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
export_fig('-pdf','-painters', '-depsc','-opengl', '-r190',savefilename_c); % PNG-nocrop'

cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab') 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Z500 PCs moving corr, Fig 5a
% 
    color_code_1=[0.1,0.1,0.1];
    color_code_2= [0.9,0.0,0.0];
    color_code_3= [0.0,0.0,0.8];


    eof_alt_nr=1;
    if eof_alt_nr==1
        load('C:\Users\Machine\matlab_storage_of_output_files\ERA-Interim_PCs_z500_lim0-360_-20_-90_1979-2011_AMJJASON_c52_annual_mean_varimax_c2.mat') %Z500
    end

MA_PCs_save_z500=MA_PCs_save;  
% MA_PCs_save_monthly_z500=MA_PCs_save;
clear MA_PCs_save MA_PCs_save_monthly  
%%%%%%%%%%%%%%%%%%%%%%% RICE
% pol_nr11=1; % polarity positive SAM, as displayed in Fig 2 2014
% pol_nr12=1; % PSA patterns as in Kidson 1988 his fig 4b,c 
% pol_nr13=-1;
pol_nr11=1; % polarity positive SAM, as displayed in Fig 2  2012
pol_nr12=1; % PSA patterns as in Kidson 1988 his fig 4b,c 
pol_nr13=-1;

fac_nr11=1; 
fac_nr12=1;
fac_nr13=1;

le=length(stacked_record_annual_Ma(80:110+add_yr)); % 1979-2011
Ma_corr=NaN(le,6);
Ma_corr4(:,1)=detrend(stacked_record_annual_Ma((80:110+add_yr),2));
Ma_corr4((1:31+add_yr),2)=detrend(MA_PCs_save_z500((1:31+add_yr),2))*pol_nr11*fac_nr11; % z500 PC1 SAM  (factor from Regression ERA-i)
Ma_corr4((1:31+add_yr),3)=detrend(MA_PCs_save_z500((1:31+add_yr),3))*pol_nr12*fac_nr12; % z500 PC2 PSA1 
Ma_corr4((1:31+add_yr),4)=detrend(MA_PCs_save_z500((1:31+add_yr),4))*pol_nr13*fac_nr13; % z500 PC3 PSA2 

%edit movingCorrelation_c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable 1: Z500 PC1
% not using function to get better overview

clear correlationTS correlationTS_p correlationTS_rlo correlationTS_rup
clear correlationTS_c correlationTS_p_c correlationTS_rlo_c correlationTS_rup_c

windowSize=sliding_w;

[N,M] = size(Ma_corr4);

correlationTS = nan(N, 1);
correlationTS_p= nan(N, 1);
correlationTS_rlo=nan(N, 1);
correlationTS_rup=nan(N, 1);


% correlationTS = nan(N, M-1);
% correlationTS_p= nan(N, M-1);
% correlationTS_rlo=nan(N, M-1);
% correlationTS_rup=nan(N, M-1);


indexColumn=1;
 %for t = windowSize+1:N % (Original)
for t = floor(windowSize/2)+1:(N-floor(windowSize/2)) % Modified
        % 6:105:
     
%      dataMatrix=Ma_corr4;
 %C = corrcoef(dataMatrix(t-windowSize:t, :)); %  (Original)
 % C = corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function 
 
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :)); % Change to centered moving correlation window
  %[C, p, rlo, rup ]= corrcoef(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :), 'alpha', alpha); 
  
  X=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),1);
  Y=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),2); 
  
  
  [C, p, rlo, rup ]= corrcoef_df(X,Y); % Change to centered moving correlation window 
 
 
%     idx = setdiff(1:M, [indexColumn]);
%     correlationTS(t, :) = C(indexColumn, idx);
%     correlationTS_p(t, :)= p(indexColumn, idx);
%     correlationTS_rlo(t, :)= rlo(indexColumn, idx);
%     correlationTS_rup(t, :)= rup(indexColumn, idx);
    
    
%    idx = setdiff(1:M, [indexColumn]);
    correlationTS(t, :) = C(1, 2);
    correlationTS_p(t, :)= p(1, 2);
    correlationTS_rlo(t, :)= rlo(1, 2);
    correlationTS_rup(t, :)= rup(1, 2);
     
    
end


clear X Y t C p rlo rup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable 2: Z500 PC2

correlationTS_c = nan(N, 1);
correlationTS_p_c= nan(N, 1);
correlationTS_rlo_c=nan(N, 1);
correlationTS_rup_c=nan(N, 1);


for t = floor(windowSize/2)+1:(N-floor(windowSize/2)) % Modified
        % 6:105:
     
%      dataMatrix=Ma_corr4;
 %C = corrcoef(dataMatrix(t-windowSize:t, :)); %  (Original)
 % C = corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function 
 
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :)); % Change to centered moving correlation window
  %[C, p, rlo, rup ]= corrcoef(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :), 'alpha', alpha); 
  
  X=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),1);
  Y=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),3); 
  
  
  [C, p, rlo, rup ]= corrcoef_df(X,Y); % Change to centered moving correlation window 
 
    
%    idx = setdiff(1:M, [indexColumn]);
    correlationTS_c(t, :) = C(1, 2);
    correlationTS_p_c(t, :)= p(1, 2);
    correlationTS_rlo_c(t, :)= rlo(1, 2);
    correlationTS_rup_c(t, :)= rup(1, 2);
        
end

clear X Y t C p rlo rup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable 3: Z500 PC3

correlationTS_c3 = nan(N, 1);
correlationTS_p_c3= nan(N, 1);
correlationTS_rlo_c3=nan(N, 1);
correlationTS_rup_c3=nan(N, 1);


for t = floor(windowSize/2)+1:(N-floor(windowSize/2)) % Modified
        % 6:105:
     
%      dataMatrix=Ma_corr4;
 %C = corrcoef(dataMatrix(t-windowSize:t, :)); %  (Original)
 % C = corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-windowSize:t, :)); % switch to UW Steig corr function 
 
 % [C, p, rlo, rup ]= corrcoef_df(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :)); % Change to centered moving correlation window
  %[C, p, rlo, rup ]= corrcoef(dataMatrix(t-floor(windowSize/2):t+floor(windowSize/2), :), 'alpha', alpha); 
  
  X=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),1);
  Y=Ma_corr4(t-floor(windowSize/2):t+floor(windowSize/2),4); 
  
  
  [C, p, rlo, rup ]= corrcoef_df(X,Y); % Change to centered moving correlation window 
 
    
%    idx = setdiff(1:M, [indexColumn]);
    correlationTS_c3(t, :) = C(1, 2);
    correlationTS_p_c3(t, :)= p(1, 2);
    correlationTS_rlo_c3(t, :)= rlo(1, 2);
    correlationTS_rup_c3(t, :)= rup(1, 2);
     
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  [correlationTS4, correlationTS4_p, correlationTS4_rlo, correlationTS4_rup]=movingCorrelation_c(Ma_corr4,sliding_w,1, 0.05);
 
 h=fig('units','inches','width',14,'height',4.5,'font','Helvetica','fontsize',16,'border','on');
   box on
   hold on 
   set(gca,'YLim',[-1 1])
 %%%%%%%%%%%%%%
shade_nr=0;
        if shade_nr==1
            
%             addpath C:\PHD\matlab_lib\Library
 %           IPO_vline(10)
            hold on
        
        end
     % PC1
       plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rlo(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_1);
    h1=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS(:,1),'-','LineWidth',4,'MarkerSize',14,'Color',color_code_1);
       plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rup(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_1);
    
    % PC2
       plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rlo_c(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_2);
    h2=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_c(:,1),'-','LineWidth',4,'MarkerSize',14,'Color',color_code_2);     
       plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rup_c(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_2);
       
    % PC3   
%        plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rlo_c3(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_3);
%     h3=plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_c3(:,1),'-','LineWidth',4,'MarkerSize',14,'Color',color_code_3);
%        plot(stacked_record_annual_Ma((80:110+add_yr),1), correlationTS_rup_c3(:,1),'--','LineWidth',2,'MarkerSize',14,'Color',color_code_3);
    
               set(gca,...
          'linewidth',3,...
          'FontWeight','bold' ); 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend
    
    hs1='SAM';
    hs2='PSA1';
%    hs3='PSA2';
    
    
    hs1_c=['r ({\delta}^1^8O, ',hs1,')'];% (',num2str(fac_nr11),'))'];
    hs2_c=['r ({\delta}^1^8O, ',hs2,')'];
 %   hs3_c=['r ({\delta}^1^8O, ',hs3,')'];
 
    
  %h_leg=legend([ h1 h2 h3],  hs1_c, hs2_c, hs3_c); 
  h_leg=legend([ h1 h2],  hs1_c, hs2_c); 
  set(h_leg, 'location', 'SouthWest','EdgeColor',[0 0 0],'FontSize',18,'FontWeight','bold' ); % ,'color','none'
  
   title_str='Z500 PCs';   
 axestext_c(0.9999999,1.0,title_str,'FontWeight','bold','FontSize',20 ); 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[80 350 60 60];
  letter='a';
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
 %%%%%%%%%%%%%%%%%%%%%%%%%
        set(gca,'XLim',[1977 2012])
      %  hl=hline([-0.5,0,0.5]);
        hl=hline([0]);
        set(hl,'LineWidth',3,'Color',[0.0,0.0,0.0]);
%        xlabel('Year','FontWeight','bold','FontSize',18 );
        ylabel('Correlation','FontWeight','bold','FontSize',20 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change to correct minus symbol
ax=gca;
set(ax,'ticklabelinterpreter','none')  %or 'tex' but not 'latex'
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
xticklabels(ax, strrep(xticklabels(ax),'-','–'));        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fig
end_str='_c2';
filename=['z500_RICE_ERA_interim_mov_corr_w_',num2str(sliding_w),'_', iso_str,'_',iso_season_str,end_str];       
filedir ='C:\Users\Machine\matlab_storage_of_output_files\figures\';


savefilename_c=strcat(filedir,filename);
% save as png 
orient landscape
cd('G:\My Drive\ISO_CFA\matlab')
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
export_fig('-pdf','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c);
cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% neff moving window
% suppl. material 
close all

show_two_nr=1; % 0,1 


PC_nr=3; % Also used for non PCs index, Nino-4 and ADP 
if PC_nr==1 % SAM 
    column_nr=2;
     letter='a';
     name_str='SAM';
elseif PC_nr==2 % PSA1
     column_nr=3;
     letter='b';
     name_str='PSA1';
elseif PC_nr==3 % Nino-4     
     letter='f';
     name_str='Nino-4';
elseif PC_nr==4 % ADP 
     letter='e';
     name_str='ADP';    
end

iso_str_c='{\delta}^1^8O';

tic
if PC_nr==1 || PC_nr==2
    wr_x=detrend(MA_PCs_save_z500((1:31+add_yr),column_nr))*pol_nr11*fac_nr11; % Z500 PC1 SAM  (factor from Regression era)
elseif PC_nr==3
    wr_x=detrend(Ma_save_nino(1:33));
elseif PC_nr==4
    wr_x=detrend(MA_ADP_save((8:40),2));
end
wr_y=detrend(MA_save((1980:2012),3));


neff_store=[];
yrs_store=[];

for i=1:1:23
    
    
 wr_mx=wr_x(i:i+10);
 wr_my=wr_y(i:i+10);
 middle_yr=1979+5+i;
 
  nc=size(wr_mx);
   r1 = corrcoef(wr_mx(2:end),wr_mx(1:end-1),'rows','pairwise');
   r2 = corrcoef(wr_my(2:end),wr_my(1:end-1),'rows','pairwise');
%    r1 = corrcoef(x(2:end),x(1:end-1),'rows','pairwise');
%    r2 = corrcoef(y(2:end),y(1:end-1),'rows','pairwise');
%    %r1 = acf(x,1);
%    %r2 = acf(y,2);
%    
  %  nv = nv*(1-r1(2)*r2(2))/(1+r1(2)*r2(2))
  
    neff = round(nc(1)*(1-r1(2)*r2(2))/(1+r1(2)*r2(2))); % Bretherton et al. 1999 their eq. 31
 
    neff_store=[  neff; neff_store]; 
    yrs_store=[  middle_yr; yrs_store];
    
    
end

 h=fig('units','inches','width',14,'height',4.5,'font','Helvetica','fontsize',16,'border','on');
   box on
   hold on 
   
  set(gca,'XLim',[1977 2012])
  if PC_nr==1 || PC_nr==2
     set(gca,'YLim',[8 18]) 
  elseif PC_nr==3 || PC_nr==4
     set(gca,'YLim',[5 13])       
  end
   
color_code_1=[0.0,0.0,0.0];
%plot(yrs_store, neff_store,'-','Color',color_code_1,'MarkerSize',20);

plot(yrs_store, neff_store,'-','Color',color_code_1,'LineWidth',2.5);

%lrc3=hline(median(neff_store),'--');
lrc3=hline(11,'--');
set(lrc3,'Color',[0.1,0.1,0.1],'LineWidth',1.5);

                 set(gca,...
          'linewidth',3,...
          'FontWeight','bold' ); 
      
        xlabel('Year','FontWeight','bold','FontSize',18 );
        ylabel('n_e_f_f','FontWeight','bold','FontSize',20 ); 
        name_str_1=name_str;
%         title(['n_e_f_f','(',iso_str_c,', ',name_str_1,')'])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % letter
  letter_pos=[80 350 60 60];
 
          TextBox = uicontrol('style','text');
          set(TextBox,'String',letter,'position',letter_pos,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize
          set(TextBox,'foregroundcolor', [0 0 0], ...
         'backgroundcolor', [1 1 1]);
 
if show_two_nr==1 && (PC_nr==1 || PC_nr==3)   
    
    if PC_nr==1
    load('C:\Users\Machine\matlab_storage_of_output_files\z500_RI_ERA_interim_mov_neff_w_11_d18O_PSA1_ExtendedWinter_c2.mat') 
    color_code_2= [0.9,0.0,0.0];
    elseif PC_nr==3
    load('C:\Users\Machine\matlab_storage_of_output_files\z500_RI_ERA_interim_mov_neff_w_11_d18O_ADP_ExtendedWinter_c2.mat') 
    color_code_2= [0.6,0.0,0.0];
    end
    
    

    plot(yrs_store, neff_store,'-','Color',color_code_2,'LineWidth',2.5);
    
    
     str_1=['n_e_f_f','(',iso_str_c,', ',name_str_1,')'];
     str_2=['n_e_f_f','(',iso_str_c,', ',name_str,')'];
     
     legend(str_1, str_2,'Location','northwest')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change to correct minus symbol
ax=gca;
set(ax,'ticklabelinterpreter','none')  %or 'tex' but not 'latex'
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
xticklabels(ax, strrep(xticklabels(ax),'-','–'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fig
end_str='_c2';
filename=['z500_RI_ERA_interim_mov_neff_w_',num2str(sliding_w),'_', iso_str,'_',name_str,'_',iso_season_str,end_str];       

filedir ='C:\Users\Machine\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,filename);
% save as png 
cd('G:\My Drive\ISO_CFA\matlab')
orient landscape
export_fig('-png','-nocrop','-painters', '-depsc','-opengl', '-r170',savefilename_c); % PNG-nocrop' 
export_fig('-pdf','-nocrop','-painters', '-depsc','-opengl', '-r190',savefilename_c);
cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab')  
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% save ts
filedir ='C:\Users\Machine\matlab_storage_of_output_files\';
savefilename_c=strcat(filedir,filename);
save(savefilename_c,'yrs_store', 'neff_store','name_str') 

