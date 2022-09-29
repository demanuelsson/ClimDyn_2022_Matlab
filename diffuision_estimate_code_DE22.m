%% SPC diffusion project code
% clear
% close all
% adopted for RI, Palmer by DE 2021 and added transfer function correction of isotopes 


% Author: Emma C. Kahle
% Last Update: 11/3/2020
% Details in corresponding publication: 
% Kahle, E. C., Holme, C., Jones, T. R., Gkinis, V., & Steig, E. J. (2018). 
% A Generalized Approach to Estimating Diffusion Length of Stable Water 
% Isotopes From Ice-Core Data. Journal of Geophysical Research: 
% Earth Surface, 123(10), 2377– 2391. https://doi.org/10.1029/2018JF004764

% This is the main code to run to create diffusion length estimates from
% ice core water isotope data. This can be utilized with other water
% isotope data sets as well.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    eval('Config_diffusion')
end
toc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up run options

% run options (1 turns a setting on, 0 turns a setting off)
run_op.q_1gauss = 1;   % use single gaussian parameterization
run_op.q_2gauss = 0;   % use double gaussian parameterization
run_op.q_plots = 0;    % make plots of each spectrum fit (1 cues plots, 0 does not plot)
run_op.remove_peak = 0;    % remove annual peak from spectra
run_op.uncertainty = 0;    % calculate uncertainty on sigma estimate
run_op.noise_add = 0;      % also calculate noise-adding technique

% save options
save_spectra = 1; % save spectra info for each isotope
save_sigma = 1;     % save diffusion lengths for each isotope

% spectra options
ARu = 100;      % number of poles for MEM analysis
method = 'MEM'; % choose 'MEM' or 'MTM'
% spec_info.win_units = 'yr';    % choose to make windows either in constant time 'yr' or constant depth 'm'
% spec_info.win_len = 100;    % in whichever unit you specify above (either years or meters)
% spec_info.step_len = 10;       % in whichever unit you specify above (either years or meters)
spec_info.win_units = 'm';    % choose to make windows either in constant time 'yr' or constant depth 'm'
spec_info.win_len = 10;    % in whichever unit you specify above (either years or meters)
spec_info.step_len = 10;       % in whichever unit you specify above (either years or meters)


% Define Holocene and Glacial regimes - can use different schemes for
% different time periods
glacial_start = 10000; % year glacial regime begins

%% Import/create spectra
% Identify data files to with water isotope data. Data must be evenly
% sampled at half-cm intervals.
tic

core_nr=2; % (1) SP (2) RI (3) Palmer (4) WDC
path_alt_nr=1; % (1) local laptop (2) P-drive

    if core_nr==1
        site='SP';
    elseif core_nr==2
        site='RI';
    elseif core_nr==3
        site='Palmer';
     elseif core_nr==4
        site='WDC';       
    end


% Set parameters for creating spectra
spec_info.save = save_spectra;
noise_toggle = 0;

% names of datafiles

if core_nr==1
    datafile17 = 'Version1_201520162017_5mm_interp.mat';
    datafile = 'SPC_ISO_halfcm_INSTAAR_linearNaN_depth_d18o_dD.txt';
    if path_alt_nr==1
    path_c1='C:\Users\Machine\matlab_lib\Library\diffusion_length_estimate\';
    elseif path_alt_nr==2
    path_c1='P:\Ice Chemistry\Palmer\Daniel\diffusion_iso\iso_data\';    
    end
elseif core_nr==2
    datafile = 'MA_RI_iso_5mm_depth_c2.mat'; % saved from RI_ice_core_iso_depth_age_c23.m line 251
    if path_alt_nr==1
    path_c1=filedir;
    elseif path_alt_nr==2
    path_c1='P:\Ice Chemistry\Palmer\Daniel\diffusion_iso\iso_data\';      
    end
elseif core_nr==3    % from Palmer_iso_005m_record_c1.m
    datafile = 'MA_Palmer_iso_5mm_depth_c1.mat';
    if path_alt_nr==1
    path_c1='C:\Users\benman\matlab_storage_of_output_files\';
    elseif path_alt_nr==2
    path_c1='P:\Ice Chemistry\Palmer\Daniel\diffusion_iso\iso_data\';     
    end
elseif core_nr==4
    datafile = 'Jones_2017_diffusion_WDC_jgrf20648-sup-0003-datasets2.mat';
    if path_alt_nr==1
    path_c1='C:\Users\benman\matlab_lib\Data\';
    elseif path_alt_nr==2
    path_c1='P:\Ice Chemistry\Palmer\Daniel\diffusion_iso\iso_data\';
    end
end



% load water isotope data
%  data17 = load(datafile17);    % import data file
 data_i = importdata([path_c1,datafile]);    % import data file

% interp d17O data to rest of data
%  d17O = interp1(data17.depthi, data17.d17Oi, data_i(:,1));

% pull all data in one matrix
%  data = [data_i(:,1),data_i(:,2),data_i(:,3),d17O];
    if core_nr==1
    data = [data_i(:,1),data_i(:,2),data_i(:,3)];  % depth, d18O, dD
    elseif core_nr==2 || core_nr==3 
    data = [data_i((2:end),1),data_i((2:end),4),data_i((2:end),3)];  % depth, d18O, dD
    elseif core_nr==4
    data = [data_i(:,1),data_i(:,3),data_i(:,2)];
    end

% Create spectra with given parameters
if path_alt_nr==1
addpath C:\Users\Machine\matlab_lib\Library\diffusion_length_estimate\
elseif path_alt_nr==2
addpath P:\Ice Chemistry\Palmer\Daniel\diffusion_iso\Kahle_functions\
end

[iso_spec] = create_spectra(method,data,spec_info,ARu,noise_toggle);
% [iso] = create_spectra(method,data,spec_info,ARu,noise_toggle);
% edit create_spectra
% function includes paths to files 

toc
%% Autofit 
tic
% Set fit parameters (for d18o)
% HOLOCENE
p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = [.1 0 .01 0 0 0];
p_ub = [100 .12 10 .5 7e-04 .015];  
parameters_d18o_hol = [p_i;p_lb;p_ub];
% GLACIAL
p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = [.1 0 .01 0 0 0];
p_ub = [100 .12 10 .5 5e-04 .015];  
parameters_d18o_glac = [p_i;p_lb;p_ub];

parameters_d18o(:,:,1) = parameters_d18o_hol;
parameters_d18o(:,:,2) = parameters_d18o_glac;

[p_fit_d18o] = auto_fit2(iso_spec.d18o, parameters_d18o, run_op, spec_info, glacial_start);
%[p_fit_d18o] = auto_fit3(iso_spec.d18o, parameters_d18o, run_op, spec_info, glacial_start);
% [p_fit_d18o] = auto_fit2(iso.d18o, parameters_d18o, run_op, spec_info, glacial_start);
% edit auto_fit2


% % Set fit parameters (for d17o)
% % HOLOCENE
% p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
% p_lb = [.01 0 .01 0 9e-05 0];
% p_ub = [100 .12 10 .1 7e-04 .015];  
% parameters_d17o_hol = [p_i;p_lb;p_ub];
% % GLACIAL
% p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
% p_lb = [.01 0 .01 0 9e-05 0];
% p_ub = [100 .12 10 .1 2e-04 .015];  
% parameters_d17o_glac = [p_i;p_lb;p_ub];

% parameters_d17o(:,:,1) = parameters_d17o_hol;
% parameters_d17o(:,:,2) = parameters_d17o_glac;
% 
% [p_fit_d17o] = auto_fit(iso_spec.d17o, parameters_d17o, run_op, spec_info, glacial_start);

% Set fit parameters (for dD) 
% HOLOCENE
p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = [.1 0 .01 0 1e-4 0];
p_ub = [100 .12 10 .5 1e-01 .015];  
parameters_dD_hol = [p_i;p_lb;p_ub];
% GLACIAL
p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = [.1 0 .01 0 1e-4 0];
p_ub = [100 .12 10 .5 1e-02 .015];  
parameters_dD_glac = [p_i;p_lb;p_ub];

parameters_dD(:,:,1) = parameters_dD_hol;
parameters_dD(:,:,2) = parameters_dD_glac;

% run fit
[p_fit_dD] = auto_fit2(iso_spec.dD, parameters_dD, run_op, spec_info, glacial_start);
% [p_fit_dD] = auto_fit2(iso.dD, parameters_dD, run_op, spec_info, glacial_start);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edit auto_fit2
%
if run_op.q_1gauss==1
    fit_str='gauss1';
elseif run_op.q_2gauss==1
    fit_str='gauss2';    
end

toc

% %% autofit, but as code
% tic 
% % loop 
% % Auto fit function
% % Use the dD and d18o spectra structures with the nonlinear least squares
% % fit routine to match a Gaussian (or summation of Gaussians) to each
% % spectra window
% 
% % Assign depths/ages of each window
% %age = iso.P(1,:); 
%   age = iso_spec.d18o.P(1,:);  % to run manually DE
% % half-cm sample size
% ss = .005;     
% 
% if run_op.q_1gauss == 1
% % p_fit = zeros(size(iso.P,2),5,3);    % initialize matrix to hold fit values
%  p_fit = zeros(size(iso_spec.d18o.P,2),5,3);      % to run manually DE
% 
% elseif run_op.q_2gauss == 1
% %p_fit = zeros(size(iso.P,2),7,3);    % initialize matrix to hold fit values
% p_fit = zeros(size(iso_spec.d18o.P,2),7,3); 
% end
% p_fit(:,1,1) = age;                 % plug in depth values
% 
% % iso_mat3d = iso.P;                  % create 3d matrix with Power and confidence levels
% % iso_mat3d(:,:,2) = iso.cu;          % to allow fit for 95% confidence of power estimate
% % iso_mat3d(:,:,3) = iso.cl;          % as well
% 
% 
% iso_mat3d = iso_spec.d18o.P;                  % create 3d matrix with Power and confidence levels
% iso_mat3d(:,:,2) = iso_spec.d18o.cu;          % to allow fit for 95% confidence of power estimate
% iso_mat3d(:,:,3) = iso_spec.d18o.cl;          % as well
% 
% 
% % calculate uncertainty of sigma if set to do so
% if run_op.uncertainty == 1
%     m_range = 1:3;
% elseif run_op.uncertainty == 0
%     m_range = 1;
% end
% 
% % calculate fit for each spectra
% for m = m_range
% %for n = 1:size(iso.P,2) % 34
% for n = 1:size(iso_spec.d18o.P,2) % 34  
%     
%     % if spectra contains nans (ie is d17o before there is data), do not
%     % attempt to fit this spectra
% %     if sum(isnan(iso.P(:,n))) > 0
%     if sum(isnan(iso_spec.d18o.P(:,n))) > 0      
%         
%         continue
%     end
%     
% % k = 2*pi*iso.f(3:end);     % import frequency data
%  k = 2*pi*iso_spec.d18o.f(3:end);     % import frequency data
% f = k/(2*pi);               % define frequency vector
% 
%     if m == 1
%         G_data = iso_mat3d(3:end,n,m);   % import power data, ignore depth and first 2 power values
%     % set up ranges for uncertainty
%     elseif m == 2      
%         % for holocene spectra
%         if age(n) <= glacial_start
%             
% %         divide_start_up = find(f>=5,1,'first');% original
% %         divide_end_up = find(f>=7,1,'first');
%         
% %         divide_start_up = find(f>=8,1,'first'); % DE
% %         divide_end_up = find(f>=10,1,'first');
%         divide_start_up = find(f>=2,1,'first'); % DE
%         divide_end_up = find(f>=4,1,'first');        
%         
%         G_data = [iso_mat3d(3:divide_start_up-1,n,2);...
%             linspace(iso_mat3d(divide_start_up,n,2),iso_mat3d(divide_end_up,n,3),divide_end_up-divide_start_up)';...
%             iso_mat3d(divide_end_up:end,n,3)];
%         % for glacial spectra
%         elseif age(n) >glacial_start
%         divide_start_up = find(f>=8,1,'first');
%         divide_end_up = find(f>=11,1,'first');
%         G_data = [iso_mat3d(3:divide_start_up-1,n,2);...
%             linspace(iso_mat3d(divide_start_up,n,2),iso_mat3d(divide_end_up,n,3),divide_end_up-divide_start_up)';...
%             iso_mat3d(divide_end_up:end,n,3)];
%         end
%     elseif m == 3
%        % holocene
%         if age(n) <= glacial_start
%         divide_start_up = find(f>=5,1,'first');
%         divide_end_up = find(f>=7,1,'first');
%         G_data = [iso_mat3d(3:divide_start_up-1,n,3);...
%             linspace(iso_mat3d(divide_start_up,n,3),iso_mat3d(divide_end_up,n,2),divide_end_up-divide_start_up)';...
%             iso_mat3d(divide_end_up:end,n,2)];
%         % glacial
%         elseif age(n) >glacial_start
%         divide_start_up = find(f>=8,1,'first');
%         divide_end_up = find(f>=11,1,'first');
%         G_data = [iso_mat3d(3:divide_start_up-1,n,3);...
%             linspace(iso_mat3d(divide_start_up,n,3),iso_mat3d(divide_end_up,n,2),divide_end_up-divide_start_up)';...
%             iso_mat3d(divide_end_up:end,n,2)];
%         end
%     end
%    
%     
% 
% 
% % % even log spacing code from Tyler Jones
% % % wavelength band to test - full range of frequencies
% % freq_high = 1/100;     %meters
% % freq_low = 1/2;   %meters
% % 
% % % number of evenly spaced logarithmic points to test withn wavelength band
% % % this can cause an error if too large
% % n_log = 16;
% % 
% % % equal log-frequency spacing, for power law fit
% %  
% % % determine b, for 10^b, for each wavelegth
% %     b_high = log10(abs(1/freq_low));
% %     b_low = log10(abs(1/freq_high));
% %  
% %     % Create a vector of 'n_log' logarithmically spaced points in the interval [10^b_high,10^b_low].
% %     logspacing = logspace(b_high,b_low,n_log);
% %     
% %     % preallocate vectors
% %     f_eq = ones(n_log-1,1);
% %     P_eq = ones(n_log-1,1);
% %     
% %     for i=1:n_log-1
% %         i_eq = find(f >= logspacing(i) & f < logspacing(i+1));    
% %         f_eq(i) = mean(f(i_eq));
% %         P_eq(i) = mean(G_data(i_eq));
% %     end
% % % end even log spacing code from Tyler Jones
% 
% % temporary assignment to debug the rest of code
% f_eq = f;
% k_eq = 2*pi*f_eq;   % logspace wavenumber
% P_eq = G_data;
%     
% % remove annual peak if specified
% if run_op.remove_peak ==1
%     % remove annual peak - NEED TO MAKE THIS WORK FOR DEPTH WORLD
%     % cut out one-year signal to find more accurate regression
%     f_conc = f_eq;
%     peak_start = nan(size(iso.P,2),1);
%     peak_end = nan(size(iso.P,2),1);
% %     peak_start = nan(size(iso_spec.d18o.P,2),1);
% %     peak_end = nan(size(iso_spec.d18o.P,2),1);    
%     
%     peak_start(1) = find (f_conc >= 3.8,1,'first'); % find beginning of annual peak
%     peak_end(1) = find (f_conc >= 4.7,1,'first'); % find end of annual peak
%     peak_start(2) = find (f_conc >= 3.6,1,'first'); 
%     peak_end(2) = find (f_conc >= 4.7,1,'first');
%     peak_start(3) = find (f_conc >= 3.6,1,'first'); 
%     peak_end(3) = find (f_conc >= 5.05,1,'first');
%     peak_start(4) = find (f_conc >= 3.6,1,'first'); 
%     peak_end(4) = find (f_conc >= 5.05,1,'first');
%     peak_start(5) = find (f_conc >= 3.6,1,'first'); 
%     peak_end(5) = find (f_conc >= 5.1,1,'first');
%     peak_start(6) = find (f_conc >= 3.6,1,'first'); 
%     peak_end(6) = find (f_conc >= 5.45,1,'first');
%     peak_start(7) = find (f_conc >= 3.6,1,'first'); 
%     peak_end(7) = find (f_conc >= 5.4,1,'first');
%     peak_start(8) = find (f_conc >= 3.6,1,'first'); 
%     peak_end(8) = find (f_conc >= 5.7,1,'first');
%     peak_start(9) = find (f_conc >= 4.7,1,'first'); 
%     peak_end(9) = find (f_conc >= 6.4,1,'first');
%     peak_start(10) = find (f_conc >= 4.7,1,'first'); 
%     peak_end(10) = find (f_conc >= 6.9,1,'first');
%     peak_start(11) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(11) = find (f_conc >= 6.7,1,'first');
%     peak_start(12) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(12) = find (f_conc >= 5.0,1,'first');
%     peak_start(13) = find (f_conc >= 5.8,1,'first'); 
%     peak_end(13) = find (f_conc >= 7.3,1,'first'); 
%     peak_start(14) = find (f_conc >= 5.9,1,'first'); 
%     peak_end(14) = find (f_conc >= 7.6,1,'first');
%     peak_start(15) = find (f_conc >= 7,1,'first'); 
%     peak_end(15) = find (f_conc >= 8.5,1,'first');
%     peak_start(16) = find (f_conc >= 7.4,1,'first'); 
%     peak_end(16) = find (f_conc >= 9.1,1,'first');
%     peak_start(17) = find (f_conc >= 8.2,1,'first'); 
%     peak_end(17) = find (f_conc >= 9.5,1,'first');
%     peak_start(18) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(18) = find (f_conc >= 5.0,1,'first');
%     peak_start(19) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(19) = find (f_conc >= 5.0,1,'first');
%     peak_start(20) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(20) = find (f_conc >= 5.0,1,'first');
%     peak_start(21) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(21) = find (f_conc >= 5.0,1,'first');
%     peak_start(22) = find (f_conc >= 8.9,1,'first'); 
%     peak_end(22) = find (f_conc >= 10.7,1,'first');
%     peak_start(23) = find (f_conc >= 7.2,1,'first'); 
%     peak_end(23) = find (f_conc >= 9.8,1,'first');
%     peak_start(24) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(24) = find (f_conc >= 5.0,1,'first');
%     peak_start(25) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(25) = find (f_conc >= 5.0,1,'first');
%     peak_start(26) = find (f_conc >= 9.2,1,'first'); 
%     peak_end(26) = find (f_conc >= 11.2,1,'first');
%     peak_start(27:end) = ones(1,size(iso.P,2)-26)*find(f_conc >= 4.9,1,'first'); 
%     peak_end(27:end) = ones(1,size(iso.P,2)-26)*find(f_conc >= 5.0,1,'first');
% %     peak_start(27:end) = ones(1,size(iso_spec.d18o.P,2)-26)*find(f_conc >= 4.9,1,'first'); 
% %     peak_end(27:end) = ones(1,size(iso_spec.d18o.P,2)-26)*find(f_conc >= 5.0,1,'first');    
%     
% %     peak_start(28) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(28) = find (f_conc >= 5.0,1,'first');
% %     peak_start(29) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(29) = find (f_conc >= 5.0,1,'first');
% %     peak_start(30) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(30) = find (f_conc >= 5.0,1,'first');
% %     peak_start(31) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(31) = find (f_conc >= 5.0,1,'first');
% %     peak_start(32) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(32) = find (f_conc >= 5.0,1,'first');
% %     peak_start(33) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(33) = find (f_conc >= 5.0,1,'first');
% %     peak_start(34) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(34) = find (f_conc >= 5.0,1,'first');
% %     peak_start(35) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(35) = find (f_conc >= 5.0,1,'first');
% %     peak_start(36) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(36) = find (f_conc >= 5.0,1,'first');
% %     peak_start(37) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(37) = find (f_conc >= 5.0,1,'first');
% %     peak_start(38) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(38) = find (f_conc >= 5.0,1,'first');
% %     peak_start(39) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(39) = find (f_conc >= 5.0,1,'first');
% %     peak_start(40) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(40) = find (f_conc >= 5.0,1,'first');
% %     peak_start(41) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(41) = find (f_conc >= 5.0,1,'first');
% %     peak_start(42) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(42) = find (f_conc >= 5.0,1,'first');
% %     peak_start(43) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(43) = find (f_conc >= 5.0,1,'first');
% %     peak_start(44) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(44) = find (f_conc >= 5.0,1,'first');
% %     peak_start(45) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(45) = find (f_conc >= 5.0,1,'first');
% %     peak_start(46) = find (f_conc >= 4.9,1,'first'); 
% %     peak_end(46) = find (f_conc >= 5.0,1,'first');
% %     f_conc(peak_start(n):peak_end(n)) = []; % remove annual peak
% %     data_conc = G_data; % rename ln(P) vector
% %     data_conc(peak_start(n):peak_end(n)) = []; % remove annual peak
%     peak_length = length(f_conc(peak_start(n):peak_end(n))); % data points in peak
%     P_eq(peak_start(n):peak_end(n)) = linspace(P_eq(peak_start(n)),P_eq(peak_end(n)),peak_length) ; % remove annual peak
% end
% 
% % Single Gaussian fit    
% if run_op.q_1gauss == 1
%     
% % define model for fit
% % G = @(p) p(1)*exp(-k_eq.^2*p(2)^2)...      % Gaussian diffusion function
% %     + p(3).^2*ss./abs(1-p(4)*exp(-1i*k_eq*ss)).^2; % Red noise function
% 
% G = @(p) p(1)*exp(-k_eq.^2*p(2)^2)...      % Gaussian diffusion function
%      + p(3).^2*ss./abs(1-p(4)*exp(-1i*k_eq*ss)).^2; % Red noise function
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%DE
% sigma = 0.1; % Noise level
% s = @(y) (G(y)).^2;
% wf = @(y) (s(y).^ 2)./(s(y).^ 2 + sigma^2 + eps); % like Johnsen 1997,   
% % s = @(p) (G(p)).^2;
% % wf = @(p) (s(p).^ 2)./(s(p).^ 2 + sigma^2 + eps); % like Johnsen 1997,  
% 
% 
% load('C:\Users\benman\matlab_storage_of_output_files\MA_RI_iso_5mm_depth_c2.mat') % RI data
% x=MA_save((4000:8000),4);
% addpath C:\Users\benman\matlab_lib\Library\gspbox-0.7.5\gspbox
% gsp_start;
% 
% N = length(MA_save((4000:8000),4)); % Number of nodes
% %N = length(iso_spec.d18o.P(3:end,1)); % Number of nodes
% G1 = gsp_sensor(N);
% G1 = gsp_compute_fourier_basis(G1);
% G1 = gsp_adj2vec(G1);
% 
% xbar = gsp_filter(G1,wf,x);
% %%%%%%%%%%%%%%%%%%
% 
% % H = @(p)
% % calculate residual
% residual = @(p) log10(G(p)) - log10(P_eq);
% 
% %Initial parameters and limits for 1 gaussian fit
% % for holocene spectra
% if age(n) <= glacial_start
% % p_i = parameters(1,1:4,1);     % [P_1, sigma_1, sigma_noise, ar1]
% % p_lb = parameters(2,1:4,1);
% % p_ub = parameters(3,1:4,1);
% 
% p_i = parameters_d18o(1,1:4,1);     % [P_1, sigma_1, sigma_noise, ar1]
% p_lb = parameters_d18o(2,1:4,1);
% p_ub = parameters_d18o(3,1:4,1);
% 
% 
% 
% % for glacial spectra
% elseif age(n) > glacial_start
% p_i = parameters(1,1:4,2);     % [P_1, sigma_1, sigma_noise, ar1]
% p_lb = parameters(2,1:4,2);
% p_ub = parameters(3,1:4,2);
% end
% 
% % Make auto fit
% options = optimset('Display', 'off');
% p_fit_temp = lsqnonlin(residual,p_i,p_lb,p_ub,options); % solves non-linear least squares problems
% p_fit(n,2:size(p_fit_temp,2)+1,m) = p_fit_temp;    % plug in current fit values to p_fit matrix
% 
% % Plot fits, if specified
% if run_op.q_plots == 1
%     
% gauss_1 = p_fit_temp(1)*exp(-k_eq.^2*p_fit_temp(2)^2);    % Diffusion Gaussian fit
% noise = p_fit_temp(3).^2*ss./abs(1-p_fit_temp(4)*exp(-1i*k_eq*ss)).^2;   % noise fit
% 
% % Plot fit results
% if isnan(iso.P(:,1)) > 0  % if iso = d17o   
%     fig(n+length(iso.P(1,:)),'units','inches','width',8.5,'height',5,'font','Helvetica','fontsize',15);
% elseif iso.P(4,1) > 1    % if iso = dD, plot in different figures
%     fig(n+2*length(iso.P(1,:)),'units','inches','width',8.5,'height',5,'font','Helvetica','fontsize',15);
% else                    % if iso = d18o
%     fig(n,'units','inches','width',8.5,'height',5,'font','Helvetica','fontsize',15);
% end
% 
% %semilogy(f_eq,P_eq,'b','LineWidth',1.2)     % plot log-binned spectrum
% semilogy(f,G_data,'b','LineWidth',1.2)     % plot original spectrum
% hold on
% if m ==1    % plot solid lines for actual data and fits
% semilogy(f_eq,G(p_fit_temp),'r','LineWidth',1.2)
% semilogy(f_eq,gauss_1,'color', [198 120 15] ./ 255,'LineWidth',1.2)
% semilogy(f_eq,noise,'k','LineWidth',1.2)
% 
% 
% else        % plot dashed lines for uncertainty levels
%     semilogy(f_eq,G(p_fit_temp),'r','LineWidth',0.7,'linestyle','--')
%     hold on
%     semilogy(f_eq,gauss_1,'color', [198 120 15] ./ 255,'LineWidth',0.7,'linestyle','--')
%     semilogy(f_eq,noise,'k','LineWidth',0.7,'linestyle','--')
% end
% xlabel('Frequency [cycles/m]')
% ylabel(['PSD [',char(8240),'^2 \cdot m]'])
% title(['Age = ' num2str(age(n)) ' yr BP'])
% legend('Data','Total Fit', 'Gaussian', 'Noise')
% 
% %%%%added by DE%%%%%%%%%
% 
% 
% save_MA_nr=1;
% 
% if save_MA_nr==1 && m ==1
% filedir ='C:\Users\benman\matlab_storage_of_output_files\';    
% 
% 
%    count_str=[];
%    site='RI';
%    end_str='_c1';
%    chem_str='iso';
%    type_str='depth';
%    
%    if round(age(n))<10
%    name_str1=[site,'_', chem_str,'_fit_and_noise','_00',num2str(round(age(n))),'_',type_str,end_str];
%    
%    elseif round(age(n))>=10 && round(age(n))<100
%    name_str1=[site,'_', chem_str,'_fit_and_noise','_0',num2str(round(age(n))),'_',type_str,end_str];
%    
%    else
%    name_str1=[site,'_', chem_str,'_fit_and_noise','_',num2str(round(age(n))),'_',type_str,end_str];
%    end
% 
%  savefilename =[filedir, name_str1 '.mat']; 
%  save(savefilename,'f_eq','gauss_1','noise'); % 
% 
% 
% end
% 
% %%%%%%%%%%%
% 
% 
% 
% if iso.P(4,1) > 1       % determine if we're scaling axes for dD
%     axis([0 100 5*10^-6 5*10^1])
% else                    % of if we're scaling axes of d18o or d17o
%     axis([0 100 7*10^-7 10^0])
% end
% set(gca,'YminorTick','off')
% 
% end
% end
% 
% 
% % % Double Gaussian Fit
% % if run_op.q_2gauss == 1
% % 
% % % Define model for fit
% % G = @(p) p(1)*exp(-k_eq.^2*p(2)^2)...      % Gaussian diffusion function
% %     + p(3).^2*ss./abs(1-p(4)*exp(-1i*k_eq*ss)).^2 ... % Red noise function
% %     + p(5)*exp(-k_eq.^2*p(6)^2);               % second gaussian function'
% %     
% % 
% % % calculate residual of model - data
% % residual = @(p) log10(G(p)) - log10(P_eq);
% 
% 
% end
% end
% 
% toc
%% save fig (for option with plot ticked)
tic
save_fig_nr=0;

if save_fig_nr==1
  

filename=['Diffusion_spectra_',site,'_Kahle_',fit_str];
%filedir ='C:\Users\benman\matlab_storage_of_output_files\figures\';
savefilename_c=strcat(filedir,'figures\',filename);

 % save as png 
orient landscape
 export_fig('-png','-painters', '-depsc','-opengl', '-r190',savefilename_c);  % ,'-nocrop'         
 export_fig('-pdf','-painters', '-depsc','-opengl', '-r190',savefilename_c);  % ,'-nocrop'  

end


% save all fit parameters
save_MA_nr=0;

if save_MA_nr==1
%filedir ='C:\Users\benman\matlab_storage_of_output_files\';    
   end_str='_c1';
   chem_str='iso';
   type_str='depth';
   name_str1=[site,'_', chem_str,'_diffusion_parameters','_',type_str,end_str];


 savefilename =[filedir, name_str1 '.mat']; 
 save(savefilename,'p_fit_d18o','p_fit_dD'); % 


end

toc
%% Run noise_adding technique if specified at top
if run_op.noise_add == 1
    noise_toggle = 1;
    
    % create spectra with noise added
    [dD_NA, d18o_NA,~,win_params] = create_spectra_WDCagewin(method,datafile,win_len,step_len,ARu,noise_toggle);
    
    % autofit
    % d18o
    p_i = [1 .06 2.5e-05 .03];     % [P_1, sigma_1, sigma_noise, ar1]
    p_lb = [.1 0 .01 0];
    p_ub = [100 .12 10 .5];  
    parameters_d18o_NA = [p_i;p_lb;p_ub];
    [p_fit_d18o_NA] = auto_fit(d18o_NA, parameters_d18o_NA, 1, 0, q_plots, remove_peak, uncertainty);

    % dD
    p_i = [10 .06 5 .1];     % [P_1, sigma_1, sigma_noise, ar1]
    p_lb = [.1 0 .01 0];
    p_ub = [10^4 .12 100 .5];  
    parameters_dD_NA = [p_i;p_lb;p_ub];
    [p_fit_dD_NA] = auto_fit_compare_Tyler(dD_NA, parameters_dD_NA, 1, 0, q_plots, remove_peak, uncertainty);
    % edit auto_fit_compare_Tyler
    
end


%% Plot results    [run thining correction below then run figs]
% Diff length
% dD fits
auto_age_dD = iso_spec.dD.P(1,:);
auto_sigmaD = p_fit_dD(:,3,1);
auto_sigma_upD = p_fit_dD(:,3,2);
auto_sigma_dnD = p_fit_dD(:,3,3);

% auto_age_NA = dD_NA.P(1,:);
% auto_sigmaD_NA = p_fit_dD_NA(:,3,1);
% auto_sigma_upD_NA = p_fit_dD_NA(:,3,2);
% auto_sigma_dnD_NA = p_fit_dD_NA(:,3,3);

% d18o fits
auto_age_d18o = iso_spec.d18o.P(1,:);
auto_sigma18o = p_fit_d18o(:,3,1);
auto_sigma_up18o = p_fit_d18o(:,3,2);
auto_sigma_dn18o = p_fit_d18o(:,3,3);

% auto_age_NA = d18o_NA.P(1,:);
% auto_sigma18o_NA = p_fit_d18o_NA(:,3,1);
% auto_sigma_up18o_NA = p_fit_d18o_NA(:,3,2);
% auto_sigma_dn18o_NA = p_fit_d18o_NA(:,3,3);

% % d17o fits
% d17o_start = find(p_fit_d17o(:,3,1)>0, 1, 'first'); % find start of d17o data
% auto_age_d17o = iso_spec.d17o.P(1,d17o_start:end);
% auto_sigma17o = p_fit_d17o(d17o_start:end,3,1);
% auto_sigma_up17o = p_fit_d17o(d17o_start:end,3,2);
% auto_sigma_dn17o = p_fit_d17o(d17o_start:end,3,3);
% 
% % adjust uncertainty ranges to account for any upper bounds that actually
% % came out as lower
% for ii = 1:size(auto_sigma17o,1)
%     if auto_sigma_up17o(ii)>auto_sigma17o(ii)   % if upper bound is actually lower
%         auto_sigma_up17o(ii) = auto_sigma17o(ii) - (auto_sigma17o(ii-1)-auto_sigma_up17o(ii-1));   % set it to uncertainty of previous window
%     end
% end

% auto_age_NA = d17o_NA.P(1,:);
% auto_sigma17o_NA = p_fit_d17o_NA(:,3,1);
% auto_sigma_up17o_NA = p_fit_d17o_NA(:,3,2);
% auto_sigma_dn17o_NA = p_fit_d17o_NA(:,3,3);


% Plot d18o and dD for NA and DG - 1 panel figure
fig('units','inches','width',10,'height',6,'font','Helvetica','fontsize',16);
if run_op.uncertainty == 1
% ax1 = shadedErrorBar(auto_age_d18o,auto_sigma18o,[auto_sigma18o-auto_sigma_dn18o auto_sigma_up18o-auto_sigma18o]','lineProps','-b','patchSaturation',0.5); 
ax1 = shadedErrorBar(auto_age_d18o,auto_sigma18o,[auto_sigma18o-auto_sigma_dn18o auto_sigma_up18o-auto_sigma18o]','lineProps','-b','patchSaturation',0.5); 


% hold on
% ax2 = shadedErrorBar(auto_age_NA/1000,auto_sigma18o_NA,[auto_sigma18o_NA-auto_sigma_dn18o_NA auto_sigma_up18o_NA-auto_sigma18o_NA]','-ok',0.5);   
elseif run_op.uncertainty == 0
ax1 = plot(auto_age_d18o,auto_sigma18o,'b','LineWidth',1);

 hold on
% plot(auto_age_d18o,auto_sigma18o./lamb_vq,'b','LineWidth',2);% correction
% lamb from below

auto_sigma18o_de=auto_sigma18o./thinning_func_vq;

plot(auto_age_d18o,auto_sigma18o_de,'b','LineWidth',2);
% ax2 = plot(auto_age_NA/1000,auto_sigma18o_NA);
end
% if run_op.uncertainty ==1
% elseif run_op.uncertainty ==0
% end

% % d17o
% hold on
% if run_op.uncertainty == 1
% ax2 = shadedErrorBar(auto_age_d17o/1000,auto_sigma17o,[auto_sigma17o-auto_sigma_dn17o auto_sigma_up17o-auto_sigma17o]','-c',0.5);   
% % hold on
% % ax2 = shadedErrorBar(auto_age_NA/1000,auto_sigma17o_NA,[auto_sigmad17o_NA-auto_sigma_dn17o_NA auto_sigma_up17o_NA-auto_sigma17o_NA]','-ok',0.5);   
% elseif run_op.uncertainty == 0
% ax2 = plot(auto_age_d17o/1000,auto_sigma17o,'-c','LineWidth',1);
% % ax2 = plot(auto_age_NA/1000,auto_sigma17o_NA);
% end


% dD
if run_op.uncertainty == 1
ax3 = shadedErrorBar(auto_age_dD,auto_sigmaD,[auto_sigmaD-auto_sigma_dnD auto_sigma_upD-auto_sigmaD]','lineProps','-r','patchSaturation',0.5);   
% hold on
% ax2 = shadedErrorBar(auto_age_NA/1000,auto_sigmaD_NA,[auto_sigmaD_NA-auto_sigma_dnD_NA auto_sigma_upD_NA-auto_sigmaD_NA]','-ok',0.5);   
elseif run_op.uncertainty == 0
ax3 = plot(auto_age_dD,auto_sigmaD,'-r','LineWidth',1);
% ax2 = plot(auto_age_NA/1000,auto_sigmaD_NA);

% hold on
% plot(auto_age_dD,auto_sigmaD./lamb_vq,'r','LineWidth',2);
auto_sigmaD_de=auto_sigmaD./thinning_func_vq;
plot(auto_age_dD,auto_sigmaD_de,'r','LineWidth',2);


end
%xlabel('Age [ka]')
%xlabel('Age [a]')




xlabel('Depth [m]')
ylabel('Diffusion Length [m]')
title([site,' ',fit_str])

%legend('\delta^{18}O','\delta^{17}O','\delta D')
legend('\delta^{18}O','destrained','\delta D','destrained')

%    set(gca,'XLim',[0 63]) % 1800 date
   set(gca,'XLim',[0 300])  
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_fig=0;
if save_fig==1
    filename=['Diffusion_length_',site,'_spectral_Kahle_',fit_str];
    filedir ='C:\Users\benman\matlab_storage_of_output_files\figures\';
    savefilename_c=strcat(filedir,filename);

    % save as png 
    orient landscape
    % increase quality of picture by increase to -r500
    export_fig('-png','-painters', '-depsc','-opengl', '-r190',savefilename_c);  % ,'-nocrop'         
    export_fig('-pdf','-painters', '-depsc','-opengl', '-r190',savefilename_c);  % ,'-nocrop' 
end
%% diff length sigma(yr) in time [ClimDyn fig]
% 

auto_sigma18o_de_a=auto_sigma18o_de./lamb_vq;
auto_sigmaD_de_a=auto_sigmaD_de./lamb_vq;



%%%%%%%%% time

 [excel_layer_depth,reporttxt] = xlsread('G:\My Drive\ISO_CFA\Age scale 2017\rice_2700-year_timescale_and_accumulation_updatedapr2018.csv'); % 2018  May, Winstrup


% auto_age_d18o
depth_wr=auto_age_d18o;
age_vq = interp1(excel_layer_depth(:,1),excel_layer_depth(:,2),depth_wr); % depth, age, depth 5mm

fig('units','inches','width',10,'height',6,'font','Helvetica','fontsize',16);

%plot(auto_age_d18o,auto_sigma18o_de_a,'b','LineWidth',2); % [yrs]
plot(age_vq,auto_sigma18o_de_a,'b','LineWidth',2); % [yrs]

set(gca,...
          'linewidth',2,...
          'FontWeight','bold' );

hold on
%plot(auto_age_d18o,auto_sigmaD_de_a,'r','LineWidth',2);
plot(age_vq,auto_sigmaD_de_a,'r','LineWidth',2);
ylabel('Diffusion Length [yr]')
xlabel('Year [A.D.]')
legend('\delta^{18}O','\deltaD')

save_fig=1;
if save_fig==1
    filename=['Diffusion_length_yr_',site,'_spectral_Kahle_',fit_str];
%    filedir ='C:\Users\benman\matlab_storage_of_output_files\figures\';
    savefilename_c=strcat(filedir,'figures\',filename);

    % save as png 
    orient landscape
    cd('G:\My Drive\ISO_CFA\matlab')
    export_fig('-png','-painters', '-depsc','-opengl', '-r190',savefilename_c);  % ,'-nocrop'         
    export_fig('-pdf','-painters', '-depsc','-opengl', '-r190',savefilename_c);  % ,'-nocrop'
    cd('G:\My Drive\ClimDyn_oct2022_R4\ClimDyn_R4_2022_Matlab\ClimDyn_R4_2022_Matlab') 
end


% save_nr=1; % (1/0)

% if save_sigma==1
%    filename=['RI_diffusion_length_spectral_Kahle_',fit_str];   
%     filedir ='C:\Users\benman\matlab_storage_of_output_files\';
%     save([filedir,filename],'iso_spec','auto_sigma18o','auto_sigmaD','auto_age_d18o');
% end

%% Accumulation, layer thickness lambda  
% correcting layer thinning
% from here and below added by DE
% RI data
tic
if path_alt_nr==1 
[excel_layer_depth,reporttxt] =...
    xlsread('G:\My Drive\ISO_CFA\Age scale 2017\rice_2700-year_timescale_and_accumulation_updatedapr2018.csv');
elseif path_alt_nr==2
[excel_layer_depth,reporttxt] =...
    xlsread('P:\Ice Chemistry\Palmer\Daniel\diffusion_iso\iso_data\rice_2700-year_timescale_and_accumulation_updatedapr2018.csv');    
end
lamb_c1=excel_layer_depth((1:2012),6);% layer width, not corrected for thinning, from Winstrup 2017 [m], pg. 760
% (:,7) density
thinning_func=excel_layer_depth((1:2012),8); %thinning annual
lamb_corrected=lamb_c1./thinning_func;% devide with thinning function explained in Winstrup [m]

density=excel_layer_depth((1:2012),7)*1000; % density [g/cm3]      *1000 to get kg/m3

depth=excel_layer_depth((1:2012),1);
depth_sigma=auto_age_d18o';% from autofit 
 
lamb_vq = interp1( depth , lamb_corrected , depth_sigma); % Lambda, wavelength, layer thickness (check)
thinning_func_vq = interp1( depth , thinning_func , depth_sigma);
density_vq = interp1( depth , density , depth_sigma) ;

% lamb_vq=lamb_vq./density_vq; % how to calculate ice eqivalent

toc
  
%% guassian filter, loop through each ice core interval
 close all
 tic 
nr=length(auto_age_d18o);
 
%  del_x_store_store=NaN(nr,2);
 del_x_store=[];
 
%  del_x_measured=MA_save_sub_an((st_ind:end),column_nr);
  % loop
%   direction_type_nr=1; % (1) forward (2) back-ward diffusion, not right eq for this  
%     gaussFilter = exp(-(2 * pi * auto_sigma18o/lamb_vq )^ 2);
  
 gaussFilter_store=[]; 
 ind_c_all=[];
 
 
 filter_alt=7; % transfer function alt (1)

 
 if filter_alt==6 || filter_alt==7
 files = subdir('C:\Users\benman\matlab_storage_of_output_files\RI_iso_fit_and_noise*.mat'); % saved from auto_fit
 end
 

 
    for ib=1:nr-14
           if filter_alt==1
%            gaussFilter = exp(-2*( pi * auto_sigma18o(ib)/lamb_vq(ib))^ 2);                %[Holme eq. 10       (1) Jones eq 2 for converstion from m to yrs], 
%         gaussFilter = exp(-(2 * pi * auto_sigma18o(ib)*f)^ 2); % (1) Jones eq 2, Holme eq. 10,
%               
%          
           elseif filter_alt==12
           gaussFilter = exp(-2*(pi * auto_sigma18o(ib)/(0.05* thinning_func_vq(ib)) )^ 2);
              % Gkinis et al. 2014
           elseif filter_alt==2
           gaussFilter = (1/(auto_sigma_dn18o(ib)*sqrt(2*pi)))*exp(-(depth_sigma(ib)^2/(2*(auto_sigma18o(ib)^2)))); % (2) Kahle eq 6, just zero
%           Kuttel eq.4
           elseif filter_alt==3
            gaussFilter =sqrt(6/(pi*auto_sigma_dn18o(ib)^2))*exp(-depth_sigma(ib).^2/(auto_sigma_dn18o(ib)^2));
            
           elseif filter_alt==4 
             gaussFilter = fspecial('gaussian',[10 1],auto_sigma_dn18o(ib)); % vector
             % yfilt = conv (MA_save(ind_c,4), gaussFilter.^0.75,'same');
           elseif filter_alt==5 % same as alt 1, but vector   
            sz = 30;    % length of gaussFilter vector
            x = linspace(-sz / 2, sz / 2, sz); 
            gaussFilter = exp(-2*x.^2*(pi * auto_sigma18o(ib)/0.05 )^ 2); % same as alt 1, but vector 
%           gaussFilter = exp(-x.^2/(2* auto_sigma18o(ib)^ 2)); % as online example
           % https://stackoverflow.com/questions/6992213/gaussian-filter-on-a-vector-in-matlab
%            gaussFilter = gaussFilter / sum (gaussFilter); % normalize 

            elseif filter_alt==6 % Johnsen 1977, Gkinis 2011  only wiener part
                
               filename=files(ib).name;
               load(filename); 
               gaussFilter=(abs(gauss_1).^2)./((abs(gauss_1).^2)+(abs(noise)).^2); % Wiener filter
         %      gaussFilter=real(ifft2(gauss_1))./((abs(gauss_1).^2)+(abs(noise)).^2); % Wiener filter
         
         elseif filter_alt==7 

             
            f=iso_spec.d18o.f(3:end);
            k = 2*pi*iso_spec.d18o.f(3:end);  
%             f=0.05;
    if core_nr==1
            gaussFilter = exp(-2*f.^2*(pi * auto_sigma18o(ib))^ 2);  % Holme eq (10)
    else
%             gaussFilter = exp(-2*f.^2*(pi * auto_sigma18o(ib)/lamb_vq(ib))^ 2);  % Holme eq (10)

            auto_sigma18o_de_ib=auto_sigma18o(ib)/thinning_func_vq(ib);
            
%             gaussFilter = exp(-0.25*f.^2* pi^2 *auto_sigma18o_de_ib^ 2); % looks better
          
            gaussFilter = exp(-(1/2).*k.^2 .*auto_sigma18o_de_ib^ 2); 
     %       gaussFilter = exp(-0.25*f.^2*(pi * auto_sigma18o(ib))^ 2);  %
            % right direction but mainly increase noise  
    end
%                filename=files(ib).name;
%                load(filename); 
%                optimFilter=(abs(gauss_1).^2)./((abs(gauss_1).^2)+(abs(noise)).^2); % Wiener filter               
%                gaussFilter=gaussFilter1.*optimFilter;
               
            
           end
           
           
           
%  % f_c=((1./auto_sigma_dn18o(ib).*sqrt(2*pi))).*exp((-0.5.*(depth_sigma(ib).^2))./(auto_sigma_dn18o(ib).^2));   % or loop   Kuttel eq. 4 (for smoothing of singal) Jones eq 1
%   f_c=((1./(diff_l_x_add_c.*sqrt(2*pi)))).*exp((-0.5.*(x_c(ib).^2))./(diff_l_x_add_c.^2));   % or loop   Kuttel eq. 4 (for smoothing of singal) Jones eq 1
% 
 
%   gaussFilter_store=[gaussFilter_store; gaussFilter];
    gaussFilter_store=[gaussFilter_store, gaussFilter]; % alt6
    end
 toc
  
 %%   normalize
 
 gaussFilter_store_norm=[];
 
 for ib=1:nr-14
 % gaussFilter = gaussFilter_store(:,ib) ./ sum (gaussFilter_store,2); % normalize, check wr_c2=sum (gaussFilter_store,2)
  gaussFilter = gaussFilter_store(:,ib) / sum (gaussFilter_store(:,ib)); % normalize, similar to other functions 
  
  gaussFilter_store_norm=[gaussFilter_store_norm, gaussFilter];
  
 end
 
 gaussFilter_store=gaussFilter_store_norm;
 
 %% make correction of iso data, convolution
 %
%     close all
    % load iso
%     clear  del_x_store
  if path_alt_nr==1  
    load('C:\Users\benman\matlab_storage_of_output_files\MA_RI_iso_5mm_depth_c2.mat') % RI data
    
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  % smooth the data with a gaussian filter first   
      hsize=16;sigma=3; %
% % hsize=30;sigma=5; %
     myfilter=fspecial('gaussian',hsize ,sigma);
% % 
% % %     x_filter_gaussian = imfilter(x, myfilter, 'replicate');
%       x_gaussian_smoothed= imfilter(data(:,2), myfilter, 'replicate');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
  elseif path_alt_nr==2 
    load('P:\Ice Chemistry\Palmer\Daniel\diffusion_iso\iso_data\MA_RI_iso_5mm_depth_c2.mat') % RI data  
  end
   % for ib2=1:4
     for ib2=1:nr-14  % NaNs at the end
 %  for ib2=1:292  
 %  for ib2=1:10 
       
        
           
%              gaussFilter_c1 = gaussFilter_store(ib2) / sum (gaussFilter_store(1:293)); % normalize   NaNs at the end
                if filter_alt<=3
                gaussFilter_c1 = gaussFilter_store(ib2);
                elseif filter_alt>=4
                gaussFilter_c1 = gaussFilter_store(:,ib2); 
                
%                 gaussFilter_c1=[repmat(gaussFilter_c1(1)+0.5,10,1);gaussFilter_c1];
%                 gaussFilter_c1=gaussFilter_c1/sum(gaussFilter_c1);
                end
%              gaussFilter_c1 = gaussFilter_store_norm(ib2);
              % index for correction
              % find iso interval to correct using depth
              if ib2==1
              ind_c=find(data(:,1)<=depth_sigma(ib2) & data(:,1)>0);
              else
              ind_c=find(data(:,1)<=depth_sigma(ib2) & data(:,1)>depth_sigma(ib2-1));    
              end
              
             %  apply filter
             %  yfilt = filter (gaussFilter_c1,1, y); 
             %  yfilt = filter (gaussFilter_c1,1, MA_save(ind_c,4)); % shifted
               
%              yfilt = conv (MA_save(ind_c,4), gaussFilter_c1, 'same');
%               yfilt = cconv (MA_save(ind_c,4), 1/gaussFilter_c1);
%                yfilt = deconv (MA_save(ind_c,4), gaussFilter_c1);
             %  yfilt = deconv (nanmoving_average(MA_save(ind_c,4),10), gaussFilter_c1);
     %           yfilt = deconv (MA_save(ind_c,4), gaussFilter_c1);
            %    yfilt = deconv (y_filter_gaussian(ind_c), gaussFilter_c1); 
      %%%%%%%%%%%%%%%%%%%%%      
       % K=wiener2(MA_save(ind_c,4),[25 1]);     % should be done in freq domain 
%         K=wiener2(f,[25 1]);  %use k_eq instead? and fft, check example code
       %  yfilt = deconv (data(ind_c,2), gaussFilter_c1);     % filter not applied in freq domain
             x_gaussian_smoothed= imfilter(data(ind_c,2), myfilter, 'replicate');
       
         % yfilt = deconv (nanmoving_average(data(ind_c,2),5), gaussFilter_c1);     %   
%            yfilt = deconv (x_gaussian_smoothed, gaussFilter_c1); % makes it shorter
           yfilt = conv (x_gaussian_smoothed, gaussFilter_c1,'same');
           
      %%%%%%%%%%%%%%%%% 
%         H = deconv (data(ind_c,2), gaussFilter_c1);
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               signal_en=length(yfilt);
               signal_st=6;
                
               del_x_store=[del_x_store; yfilt];
               ind_c_all=[ind_c_all;ind_c];
          
    end
    
%
% 
%         % apply filter
%         yfilt = filter (gaussFilter,1, y); 
  
toc
% Alt1
% plot as is
%   figure
%     plot(MA_save(ind_c,1),MA_save(ind_c,4))
%  hold on
%     plot(MA_save(ind_c,1)+0.3, yfilt,'.-r') 
    
% Alt2
% remove mean
%  figure % one section
%     plot(MA_save(ind_c,1),MA_save(ind_c,4))
%  hold on
%     plot(MA_save(ind_c,1)+0.3, yfilt,'.-r')
 figure % one section

 ind_f=[signal_st:signal_en]';
 
 plot(data(ind_c(ind_f),1)+0.005,yfilt(ind_f),'.-r') 
        hold on 
%          plot(data(ind_c,2),'b')
  plot(data(ind_c,1),x_gaussian_smoothed,'b')
 
  set(gca,'YLim',[-30 -15])  

         
%          
%          figure
%           plot(nanmoving_average(MA_save(ind_c,4),5),'g')
%         hold on
%           plot( nanmoving_average(yfilt,5),'.-g') 
        
%legend('d18O','d18O_0')
legend('filtered','signal smoothed')    

%% plot corrected data

figure
plot(data(ind_c_all,1),del_x_store,'-k') 
% jump in y-axis
% values missing at the start, pad record
legend('d18O_0','location','southeast')
size(del_x_store);
xlabel('depth (m)')

%% plot transfer functions
% does not contain sigma cfa
[a, b]=size(gaussFilter_store);

figure
hold on

for i=1:b
    
plot(log(iso_spec.d18o.f(3:end)),gaussFilter_store(1:end,i),'.-')


end
xlim([0 4])
xlabel('log (freq)')