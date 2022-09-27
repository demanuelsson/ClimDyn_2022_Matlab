%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map_EOF_z500_SAT_SIC_DE22.m
% Daniel Emanuelsson
% Matlab 2021a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    addpath('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\configfiles\')
    



 eval('config_map')


end
toc        
%%   

 f1=fig('units','inches','width',12,'height',12,'font','Helvetica','fontsize',16,'border','on');
 
 
 
 
bedmap2('gl','color','k','LineWidth',1.5); 
% title('bedmap2 gl','fontname','courier')
hold on
% mapzoom('Roosevelt Island','mapwidth',3500,'inset','southwest','size',.38)
 mapzoom('Roosevelt Island','mapwidth',7000)

lat_c=-86;lon_c=-125;


bedmap2('surfc','zvals',0:1000:4000,'color','-r','linewidth',2);

bedmap2('patchshelves','LineWidth',1,'linestyle','-','facecolor',[1 1 1],'edgecolor',[0 0 0],'linewidth',1.5);

bedmap2('graticule','lats',-90:10:-60,'lons',0:30:180,'color',[0.6 0.6 0.6],'LineWidth',1.0,'linestyle',':');


 
plotm(-79.3628, -161.7009,'.','MarkerSize',26, 'color',[0.9 0.0 0]); %  Deep core from Darcy google earth file
% plotm(-79.3628, -161.7009,'o','MarkerSize',9, 'color',[0.2 0.65 0],'LineWidth',1.8)
 
 h2=plotm(-80, -165,'.','MarkerSize',20,'Color',[0.45,0.56,0.61]); % Margaret AWS

 linem([-65 -45],[-60 -60],'k','LineWidth',2)  % ABS/Weddell Sea
 linem([-75 -45],[-130 -130],'k','LineWidth',2) % ABS/Ross Sea
 linem([-70 -45],[160 160],'k','LineWidth',2) % Ross Sea/western Pacific Ocean
 linem([-70 -45],[20 20],'k','LineWidth',2)



% h4=plotm(-65.245556,-64.258056,'.r','MarkerSize',20);  % -65.245556S, -64.258056W   Faraday Station

% h4=plotm(-80, -120,'.r','MarkerSize',20); % Byrd Station



%  plotm(-79.3750,-161.875,'.m','MarkerSize',20);
% % plotm(-79.3750,-161.750,'.r','MarkerSize',20); % used
%  plotm(-79.3750,(era_lon(68)-360),'.g','MarkerSize',20);


su1='RICE ice core';
su2='AWS Maragaret';
su3='Ice core sites';
su4='Research stations';
 


fw='italic';
nr=11;


% scarlabel('Ross Ice Shelf','fontsize',nr,'fontangle','italic','fontweight','bold')

scarlabel('South Pole','fontsize',nr,'fontangle','italic','fontweight','bold')

% scarlabel('Ross Ice Shelf','fontsize',nr-1,'fontangle','italic','fontweight','bold')    


% nr_c=8;
% h33=scarlabel('Ross Ice Shelf','fontsize',nr_c,'fontangle','italic','fontweight','bold','color',[1.0 1.0 1.0]);
%             pos_c = get(h33,'position');  
%             pos_c(2)=pos_c(2)-0.008;  
%             set(h33,'pos',pos_c);

 textm(-81,-165,'Ross Ice','fontsize',nr,'fontangle','italic','FontWeight','bold','color',[0 0 0])            
 textm(-79.9,-169.5,'Shelf','fontsize',nr,'fontangle','italic','FontWeight','bold','color',[0 0 0])            
            
% textm(

%scarlabel('Ross Sea','fontsize',nr)
%scarlabel('Roosevelt Island','fontsize',nr)
% scarlabel('Siple Dome','fontsize',nr)
%scarlabel('West Antarctic Ice Sheet','fontsize',nr)
% scarlabel('Bellingshausen Sea','fontsize',nr)
% scarlabel('Wilkes Land','fontsize',nr)
% scarlabel('Amundsen Sea','fontsize',nr)
% scarlabel('Amundsen Sea','fontsize',nr)
%scarlabel('Transantarctic Mountains','fontsize',nr,'rotation',-50,'color','red','fontweight','bold','fontangle','italic')
%scarlabel('Marie Byrd Land','fontsize',nr)
%scarlabel('Victoria Land','fontsize',nr,'fontangle','italic')
%scarlabel('Byrd','fontsize',nr,'fontangle','italic')
%scarlabel('Ellsworth Land','fontsize',nr,'fontangle','italic')

%,'Box','off');

x_c1=[ -77.3 :-0.1: -78.3]'; %length(x_c1)
y_c1=[ -165 :0.1: -163.5]';
y_c2=y_c1(2:end-4); %length(y_c2)
plotm(x_c1,y_c2,'-k','LineWidth',2,'MarkerSize',10)% ,'color',[0.6 0.0 0]

 textm(-76.5,-159.3,'Roosevelt','fontsize',nr,'FontWeight','bold','color',[0.6 0.0 0])
 textm(-75.6,-163.5,'Island','fontsize',nr,'FontWeight','bold','color',[0.6 0.0 0])

% 
%  textm(-78,-160,'Roosevelt','fontsize',nr,'fontweight','bold')
%  textm(-77.2,-162,'Island','fontsize',nr,'fontweight','bold')

 textm(-75,245,'Marie Byrd Land','fontsize',nr,'fontangle','italic','rotation',-20,'fontweight','bold')


    
    textm(-78.4, -117.8,'WAIS','fontsize',nr,'fontweight','bold','fontangle','italic')
    textm(-67,-43,'Weddell Sea','fontsize',nr+2)
    textm(-71,-155,'Ross Sea','fontsize',nr+2)
    textm(-67.0, 250,'Amundsen','fontsize',nr+2)
    textm(-67.3, 246.5,'Sea','fontsize',nr+2)
    textm(-63.0, 273,'Bellingshausen','fontsize',nr+2)
    textm(-64.2, 270.2,'Sea','fontsize',nr+2)
    

textm(-74.0,266.5, 'Ellsworth Land','fontsize',nr,'fontweight','bold','fontangle','italic')    
textm( -65.4,-57, 'Antarctic Peninsula','fontsize',nr,'rotation',-30,'fontangle','italic','fontweight','bold') 
textm(-88,188,'Transantarctic Mountains','fontsize',nr,'rotation',-43,'fontangle','italic','fontweight','bold')
 

    textm(-80,32,'80{\circ}S','fontsize',nr,'FontWeight','bold')
    textm(-70,117,'70{\circ}S','fontsize',nr,'FontWeight','bold')    
    textm(-60,-59.3,'60{\circ}W','fontsize',nr,'FontWeight','bold') 

 % lon markers
 textm(-66,150,'150{\circ}E','fontsize',nr,'FontWeight','bold')
 textm(-69,-180,'180{\circ}W','fontsize',nr,'FontWeight','bold')
 textm(-68,-150,'150{\circ}W','fontsize',nr,'FontWeight','bold')
 textm(-60,-120,'120{\circ}W','fontsize',nr,'FontWeight','bold') 
 textm(-60,-90,'90{\circ}W','fontsize',nr,'FontWeight','bold') 

             
      textm(-82,66,'4000','fontsize',10,'rotation',0,'fontangle','italic') 
      textm(-72.2,100,'3000','fontsize',10,'rotation',60,'fontangle','italic')
      
         
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annotation(f1,'rectangle',...
    [0.17 0.33 0.5 0.59],...
    'LineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% letter

TextBox = uicontrol('style','text');

           set(TextBox,'String',letter,'position',letter_position,'FontWeight','bold','FontSize',34 ); % x position y position xsize ysize

            set(TextBox,'foregroundcolor', [0 0 0], ...
            'backgroundcolor', [1 1 1]);
  
     
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

    qual_str='-r250';

% save fig.  [change dir folder, or stay in the matlab dir]
% 
 

savefilename_c=strcat(filedir,filename);
% saveas(f1,[savefilename_c],'fig')

% save as png 
orient landscape

cd('G:\My Drive\ISO_CFA\matlab')
export_fig('-png','-painters', qual_str, savefilename_c); % func needs ghostscript ,'-nocrop'
export_fig('-pdf','-painters', qual_str, savefilename_c); % func needs ghostscript ,'-nocrop'
cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%