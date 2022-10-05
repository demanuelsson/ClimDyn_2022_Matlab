function [coord]=star_coord_WA(site)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antarctic coordinates
% For ice-core drill sites and AWSs
% Daniel Emanuelsson
%%%%%%%%%%%%%

if strcmp(site,'RI')==1 || strcmp(site,'RICE')
    
%-79.3628, -161.7009%  Deep core from Darcy's google earth file

    coord=[-79.3628, -161.7009];  

elseif strcmp(site,'Marg')==1 || strcmp(site,'AWSMarg')==1
    
    coord=[-80, -165];  

elseif strcmp(site,'WDC')==1 || strcmp(site,'ITASE 2000-1')
    
     coord=[-79.467, -112.085];    

%     79.467 -112.085

elseif strcmp(site,'ITASE_2001_5')==1    ||    strcmp(site, 'ITASE 2001-5')
    
%     77.06S, 89.14W

    coord=[-77.06, -89.14 ]; % updated fits 100 m magnitude contour
              
elseif strcmp(site,'Byrd')==1    
    
  %  80°S, 120°W
    
     coord=[-80, -120];
     
elseif strcmp(site,'ITASE 2000-5')==1    || strcmp(site,'ITASE_2000-5')==1  || strcmp(site,'ITASE_2000_5')
        
     coord=[-77.68, -124];
     
elseif strcmp(site,'ITASE 2000-4')==1   || strcmp(site,'ITASE_2000-4')==1   
        
     coord=[-78.08, -120.08];
     

 elseif strcmp(site,'Faraday')==1    
    
   %  -65.245556S, -64.258056W
     
     coord=[-65.245556, -64.258056 ];    
     
     
 elseif strcmp(site,'Gomez')==1    
     
    coord=[-73.59, -70.36 ];  
    
 elseif strcmp(site,'Palmer')==1 
    coord=[-73.86, -65.46 ];
    
 elseif strcmp(site,'Palmer')==1 
    coord=[-73.86, -65.46 ]; 
    
 elseif strcmp(site,'James Ross I')==1   || strcmp(site,'James Ross Island')==1 
  coord=[ -64.20, -57.68 ];
  
  
  elseif strcmp(site,'Siple Station')==1 
    coord=[-75.92, -84.25 ];  
  
  elseif strcmp(site,'Dyer Plateau')==1  
   coord=[-70.67, -64.87 ];
   
  elseif strcmp(site,'Ferrigno')==1 
   
   coord=[-74.57, -86.9 ];
   
  elseif strcmp(site,'Bruce Plateau')==1 
   
   coord=[-66.03, -64.07 ];  
   
  elseif strcmp(site,'IND33')==1 
   coord=[-71.5, 10.2 ]; 
   
   elseif strcmp(site,'Rendezvous')==1 
   coord=[-74.45, -78.16 ];   
   
   elseif strcmp(site,'Jurassic')==1 
   coord=[-74.33	,-73.06 ];   
   
end