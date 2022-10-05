
function MA_PCs_save=annual_mean_DE(yr_s, yr_e,era_year_num, MA_PC_time_series,sea_nr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% annual_mean_DE.m
% annual means 
% -= Input =-
% yr_s = start year
% yr_e = end year
% era_year_num = high-resolution time series e.g. from ERA-Inerim monthly
% values
% MA_PC_time_series = Matix or vector of values to average. e.g. monthly PCs values from
% EOFs
% -= Output
% Matrix
% (:,1) years
% (:,2.......) mean values
% DE 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_annual=yr_s:yr_e;

[a b c]=size(MA_PC_time_series);

MA_PCs_store=zeros(length(e_annual)-1,b);
MA_PCs_save=zeros(length(e_annual)-1,(b+1));

date_stack_annual=zeros(length(e_annual)-1,1);

 
 for i = 1:length(e_annual)-1
     
        if sea_nr==2 % DJF
        stacked_time =era_year_num(find(era_year_num >= (e_annual(i)-1/12) & era_year_num <(e_annual(i)+2/12)));
        Dummy_temp=MA_PC_time_series(find(era_year_num >= (e_annual(i)-1/12) & era_year_num <(e_annual(i)+2/12)),:);
        
        elseif sea_nr==3 % MAM, 345 ok
        stacked_time =era_year_num(find(era_year_num >= (e_annual(i)+2/12) & era_year_num <(e_annual(i)+5/12)));
        Dummy_temp=MA_PC_time_series(find(era_year_num >= (e_annual(i)+2/12) & era_year_num <(e_annual(i)+5/12)),:);
        
        elseif sea_nr==4 % JJA
        stacked_time =era_year_num(find(era_year_num >= (e_annual(i)+5/12) & era_year_num <(e_annual(i)+8/12)));
        Dummy_temp=MA_PC_time_series(find(era_year_num >= (e_annual(i)+5/12) & era_year_num <(e_annual(i)+8/12)),:); 
        
        elseif sea_nr==5 % SON
        stacked_time =era_year_num(find(era_year_num >= (e_annual(i)+8/12) & era_year_num <(e_annual(i)+11/12)));
        Dummy_temp=MA_PC_time_series(find(era_year_num >= (e_annual(i)+8/12) & era_year_num <(e_annual(i)+11/12)),:);
        
        
        elseif sea_nr==6 % MAMJJA, 345678 ok, check by just running what's in the ()
        stacked_time =era_year_num(find(era_year_num >= (e_annual(i)+2/12) & era_year_num <(e_annual(i)+8/12)));
        Dummy_temp=MA_PC_time_series(find(era_year_num >= (e_annual(i)+2/12) & era_year_num <(e_annual(i)+8/12)),:);       
        
        
        elseif sea_nr==1
        stacked_time =era_year_num(find(era_year_num >= e_annual(i) & era_year_num <e_annual(i+1)));
        Dummy_temp=MA_PC_time_series((find(era_year_num >= e_annual(i) & era_year_num <e_annual(i+1))),:);
        
        end
        
        date_stack_annual(i)=floor(mean(stacked_time(1:end)));
        MA_PCs_store(i,:)=nanmean(Dummy_temp((1:end),:));  

        clear Dummy_temp
        
 end
 
 MA_PCs_save(:,1)=date_stack_annual;
 MA_PCs_save(:,(2:(b+1)))=MA_PCs_store;

end

