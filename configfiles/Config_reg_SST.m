% Config_corr_SIC

data_dir='C:\Users\Machine\matlab_lib\Data\HadISST_data\';
filedir ='C:\Users\Machine\matlab_storage_of_output_files\';

 name='SST';
 type_nr=1;    % (1) annual (or seasonal but one value for each year) (2) monthly
 param_nr=1;% (1-SST/ 2 SIC)

 site='RICE';
  yr_s=1979;
  yr_e=2011; 

 sea_nr=7; % AMJJASON
 season='AMJJASON';
 PC_nr=3; % (2) SAM (3) PSA1 [use] (4) PSA2

 set_yr=1;
 
 iso_nr=7;
 iso='PCs';
 
 
 

  proj_nr=2;
  proj='mercator';
  lat1=-90;
  lat2=20;
  lon1=0;
  lon2= 359.9;
  box_use=0;  % (0/1)
  RSAS_SIC_box=0;
  
  show_colorbar=1;
  color_alt=7; 
  c_limit=0.6;
  
        x_move_colorbar=-0.048;
        y_move_colorbar=0.002;
        adj_width=0.0;
        adj_height=-0.05;
        
   lock_scalebar=1;
        
  
  co_pa=[.88 .18 .88];  % contour   
  line_w_nr=2;

letter='c';
letter_pos= [175 525 50 50];
letter_size=22;%28

coast_nr=1;

x_lab=0.96;y_lab=0.63; % unit position

show_title=0;
Nino_34=1; % show nino-34 box

label_size=20;
label_vert=-300;
label_across=-64.9;
%label_across=19.99;

