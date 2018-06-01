% This code is for selecting ground motions by the use of MCMCGMS methods
% The target means and covariances corresponding to a 
% pre-defined scenario earthquake are obtained by using Generalized Conditional Spectrum.
% Single-component or two-component motions can be selected for different objectives.

% This version is created by Yuan Shi and Peter J. Stafford
% Civil and Environmental Engineering
% Imperial College

% Last update: 13/09/2017

% Further details are provided in the following document:

%  Shi Y and Stafford PJ (2018). "Markov-chain Monte Carlo Ground-Motion Selection Algorithms for 
%  Conditional Intensity Measure Targets" Earthquake Engineering & Structural Dynamics (in revision)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Please note:
%   This code is used for presenting an example of case 1 shown in the above article for users. 
%   (i.e. records are selected to match to a target distribution constructed by spectral acceleration 
%   at different periods and 5-95% significant duration)
%   This code is intended to assist readers in coding record selection by the use of MCMCGMS methods for 
%   their own purposes. It is easy to adjust this example code to consider other intensity measures cases.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ouput variables:
% Index_GM  : cells to save records ID of selected records

% Input variables:
% Only the NGA W1 database is used in this code as an example. This databased workspace is named as 'rec_selection_meta_data.mat'.
% The database can be alternatively applied as a user's purposes with the information of spectral acceleration spectra, 5-95% signifiant duration and causal parameters.

% To illustrate a example (only Spectral acceleration and 5-95% significant duration considered)
% s_main                   : Short name for Spectral Acceleration
% ds_main                  : Short name for 5-95% Significant Duration 

% The Campbell and Bozorgnia (2008) Ground Motion Model (GMM) is used in this code for spectral acceleration. 
% The Bommer, Stafford and Alarcon (2009) Ground Motion Model is used in this code for 5-95% significant duration.
% These GMMs can be changed if desired and additional information may be required by GMMs.

% Variable definitions
% variables used for GMM
% Input                    : A structure with parameters that specify the rupture scenario

% .Mmagnitude              : Earthquake magnitude
% .Rrupture                : Earthquake distance
% .deltainput              : Average dip of the rupture place (degree)
% .lambdainput             : Rake angle (degree)
% .Vs30input               : Shear wave velocity averaged over top 30 m (m/s)
% .Zvsinput                : Depth to the 2.5 km/s shear-wave velocity horizon (km)
% .arbinput                : =1 for arbitrary component sigma 
%                            =0 for average component sigma
% .Ztor                    : Depth to the top of coseismic rupture (km)
% .Epsilon                 : Epsilon value for a target scenario


% Input.nGM                : The target number of selected records

% Specific inputs for MCMCGMS 1
% selection_method_MCMCGMS_1 : Indicator for MCMCGMS method 1; if =1 means MCMCGMS method 1 is applied otherwise not applied

% Input_MCMCGMS_1          : A structure with parameters that specify input of MCMCGMS method 1

% .niter                   : Number of iteration 
% .nRep                    : Number of replication used in energy check
% .pre_select_input        : Indicator of PDF preselection and if =0 means no PDF preselection needed (default=1; i.e. PDF preselection needed)
% .gamma                   : The value of gamma used in noise covariance
% .W_SA                    : Normalized weight factor used for SA in MCMCGMS_1
% .W_Ds_595                : Normalized weight factor used for Ds_595 in MCMCGMS_1   

% specific inputs for MCMCGMS 2
% selection_method_MCMCGMS_2 : Indicator for MCMCGMS method 2; if =1 means MCMCGMS method 2 is applied otherwise not applied

% Input_MCMCGMS_2          : A structure with parameters that specify input of MCMCGMS method 2

% .nloop                   : number of loop for optimization
% .pre_select_input        : indicator of PDF preselection and if =0 means no PDF preselection needed (default=1; i.e. PDF preselection needed)
% .mean                    : Weight factor for mean
% .std                     : Weight factor for standard deviation
% .ske                     : Weight factor for skewness
% .ske_dev_check           : Indicator for whether skewness is considered in the optimization; if =1 yes other no (default=1; i.e. skewness accounted for in the optimization)
% .W_SA                    : Normalized weight factor used for SA in MCMCGMS_2
% .W_Ds_595                : Normalized weight factor used for Ds_595 in MCMCGMS_2

% Please note the weight factor of skewness is suggested as 0.03 times weight factor of mean. This can be changed if desired.


% Please note if all weight factors for intensity measures are set as 1.0 (i.e., W_SA = W_Ds_595 = 1) 
% this means no weighting factor are explicitly applied, which is implemented in the article. 
% if normalized weight factors are used for intensity measures, sum of weight factors (i.e., W_SA + W_Ds_595 = 1) should be equal to one.   

% Binallow                 : A structure with parameters that build a bin to pre select records based on casual parameters (i.e. magnitude and distance)
% .M                       : Earthquake magnitude
% .D                       : Earthquake distance

% Scalefactor_allow        : The maximum allowable scale factor
% T1_input                 : The period of conditional spectral acceleration
% PerTgt_input             : Period range of interest 
% plot_indicator           : Indicator for plot distribution match and if = 0 means no plots (default=1; i.e. plots)
% need_distribution        : Indicator to obtain distribution of selected records and if = 0 means no distribution needed (default=1; i.e. distribution needed)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
% User inputs begin
% Short names for intensity measures
s_main              = 'SA';                    
ds_main             = 'Ds595';
%The target number of selected records     
Input.nGM           = 50;
% User inputs to specify the target earthquake rupture scenario (right lateral strike slip is assumed here)
Input.Mmagnitude    = 6.5;
Input.Rrupture      = 10;
Input.deltainput    = 90;      
Input.lambdainput   = 180;    
Input.Vs30input     = 400;     
Input.Zvsinput      = 2;          
Input.arbinput      = 0;      
Input.Ztor          = 0;           
Input.Epsilon       = 2;         
            
% User inputs for MCMCGMS method 1
selection_method_MCMCGMS_1       = 1;
Input_MCMCGMS_1.niter            = 50;         
Input_MCMCGMS_1.nRep             = 10;           
Input_MCMCGMS_1.pre_select_input = 1;
Input_MCMCGMS_1.gamma            = 1;
Input_MCMCGMS_1.W_SA             = 1.0;
Input_MCMCGMS_1.W_Ds_595         = 1.0;
% User inputs for MCMCGMS method 2
selection_method_MCMCGMS_2       = 1;
Input_MCMCGMS_2.nloop            = 2;                 
Input_MCMCGMS_2.pre_select_input = 1;  
Input_MCMCGMS_2.mean             = 1.0;            
Input_MCMCGMS_2.std              = 2.0;
Input_MCMCGMS_2.ske              = 0.03*Input_MCMCGMS_2.mean; 
Input_MCMCGMS_2.ske_dev_check    = 1;
Input_MCMCGMS_2.W_SA             = 1.0;                 
Input_MCMCGMS_2.W_Ds_595         = 1.0;

% User inputs for preselection for casual parameters and scale factor
% This preselection is not considered in this code (by changed the value this preselection can be considered if desired),
% since the study is focused on PDF preselection. 
Binallow.M                       = [0 inf];
Binallow.D                       = [0 inf];
Scalefactor_allow                = inf;                

% Other parameters to select records and evaluate selections 
T1_input            = 1;               
PerTgt_input        = [0.1,0.12,0.133,0.16,0.19,0.22,0.25,0.3,0.35,0.4,0.48,0.55,0.65,0.75,0.9,0.95,1.2,1.4,1.6,1.9,2.2,2.6,3,3.5,4.2,4.8,5.5,6.5,7.5,9];
plot_indicator      = 1; 
need_distribution   = 1;
% User inputs end here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save selection method indicator
selection_method_indicator=[selection_method_MCMCGMS_1,selection_method_MCMCGMS_2];

% put T1 in the period vector
if ~any(PerTgt_input == T1_input)
    PerTgt_input = [PerTgt_input(PerTgt_input<T1_input) T1_input PerTgt_input(PerTgt_input>T1_input)];
end

% Save location of conditional period in the period vector
Index_T1_input = find(PerTgt_input == T1_input);

% Screen database 
if Input.arbinput == 1
	[data_after_preselection,GM] = database_input_arbitrary_SA_DS595(PerTgt_input,Binallow);
else 
	[data_after_preselection,GM] = database_input_average_SA_DS595(PerTgt_input,Binallow);
end

% Target distribution
[Tgt_distro_output,Database_afterScaling] = get_target_distribution_upload(Input, s_main, ds_main, data_after_preselection, PerTgt_input, T1_input, Scalefactor_allow, GM);

fprintf('Target distribution obtained \n '); 

% MCMCGMS method 1
if selection_method_MCMCGMS_1 == 1
   [MCMCGMS_1_output] = MCMCGMS_method1_for_case_study_SA_PLUS_DS595_upload(Input_MCMCGMS_1, Tgt_distro_output, Database_afterScaling,Index_T1_input, Input);

% Selected records ID  
   Index_GM{1} = MCMCGMS_1_output.Index_optim;
end

fprintf('MCMCGMS method 1 done \n ');

% MCMCGMS method 2
if selection_method_MCMCGMS_2 == 1
   [MCMCGMS_2_output] = MCMCGMS_method2_for_case_study_SA_PLUS_DS595_upload(Input_MCMCGMS_2, Tgt_distro_output, Database_afterScaling, Index_T1_input, Input);

% Selected records ID 
   Index_GM{2} = MCMCGMS_2_output.Index_GM_after_OP;
end 

fprintf('MCMCGMS method 2 done \n ');  

% Selected records ID 
GMID_final = cell(length(selection_method_indicator),1);
% Selected records
GM_final = cell(length(selection_method_indicator),1);
% Mean of selected records
Mean_GM = cell(length(selection_method_indicator),1);
% Dispersion of selected records
Dis_GM = cell(length(selection_method_indicator),1);

for i = 1:length(selection_method_indicator)
	if selection_method_indicator(i) == 1
% Save GMID 
      GMID_final{i} = Database_afterScaling.GMID (Index_GM{i}); 
% Save records selected 
      GM_final{i} = Database_afterScaling.GM_plot(Index_GM{i},:); % format is log
% The mean and dispersion of records selected set
      Mean_GM{i} = mean(GM_final{i});
      Dis_GM{i} = std(GM_final{i});
  end
end  

if selection_method_MCMCGMS_1 == 1
% Plot the distribution match MCMCGMS method 1
              figure(1)
              loglog(Tgt_distro_output.TgtIM_plot,exp(Tgt_distro_output.TargetMean_plot),'LineWidth',2)
              hold on
              loglog(Tgt_distro_output.TgtIM_plot,exp(Mean_GM{1}),'--','LineWidth',2)
              axis([min(Tgt_distro_output.TgtIM_plot) max(Tgt_distro_output.TgtIM_plot) min(min(exp(Mean_GM{1})), min(exp(Tgt_distro_output.TargetMean_plot))), max(max(exp(Mean_GM{1})), max(exp(Tgt_distro_output.TargetMean_plot)))])
              set(gca,'XTickLabel',{'10^{-1}s' '1s' '10^{1}s' 'Ds5-95%'})
              xlabel('T (s) and Ds5-95%')
              ylabel('Median IM')
              title('Response spectrum mean match by MCMCGMS method 1 ((g) for SA and (s) for Ds5-95%)');
              set(gca,'fontsize',9,'fontname','Times New Roman')
              hleg=legend('exp(Target mean lnIM)','exp(Mean of selected lnIM)','Location','SW');
              set(hleg,'fontsize',9,'fontweight','normal','fontname','Times New Roman');
              legend boxoff
% figure for dispersion
              figure(2)
              semilogx(Tgt_distro_output.TgtIM_plot,sqrt(diag(Tgt_distro_output.TargetCov_plot))','LineWidth',2)
              hold on
              semilogx(Tgt_distro_output.TgtIM_plot,Dis_GM{1},'--','LineWidth',2)
              axis([min(Tgt_distro_output.TgtIM_plot) max(Tgt_distro_output.TgtIM_plot) 0 1])
              set(gca,'XTickLabel',{'10^{-1}s' '1s' '10^{1}s' 'Ds5-95%'})
              xlabel('T (s) and Ds5-95%')
              ylabel('Standard deviation of lnIM')
              title('Response spectrum standard deviation match by MCMCGMS method 1 ((g) for SA and (s) for Ds5-95%)');
              set(gca,'fontsize',9,'fontname','Times New Roman')
              hleg=legend('Target standard deviation of lnIM','Standard deviation of selected lnIM');
              set(hleg,'fontsize',9,'fontweight','normal','fontname','Times New Roman');
              legend boxoff               
end
if selection_method_MCMCGMS_2 == 1
% Plot the distribution match MCMCGMS method 2
              figure(3)
              loglog(Tgt_distro_output.TgtIM_plot,exp(Tgt_distro_output.TargetMean_plot),'LineWidth',2)            
              hold on
              loglog(Tgt_distro_output.TgtIM_plot,exp(Mean_GM{2}),'--','LineWidth',2)
              axis([min(Tgt_distro_output.TgtIM_plot) max(Tgt_distro_output.TgtIM_plot) min(min(exp(Mean_GM{2})), min(exp(Tgt_distro_output.TargetMean_plot))), max(max(exp(Mean_GM{2})), max(exp(Tgt_distro_output.TargetMean_plot)))])
              set(gca,'XTickLabel',{'10^{-1}s' '1s' '10^{1}s' 'Ds5-95%'})
              xlabel('T (s) and Ds5-95%')
              ylabel('Median IM')
              title('Response spectrum mean match by MCMCGMS method 2 ((g) for SA and (s) for Ds5-95%)');
              set(gca,'fontsize',9,'fontname','Times New Roman')
              hleg=legend('exp(Target mean lnIM)','exp(Mean of selected lnIM)','Location','SW');
              set(hleg,'fontsize',9,'fontweight','normal','fontname','Times New Roman');
              legend boxoff             
% figure for dispersion
              figure(4)
              semilogx(Tgt_distro_output.TgtIM_plot,sqrt(diag(Tgt_distro_output.TargetCov_plot))','LineWidth',2)
              hold on
              semilogx(Tgt_distro_output.TgtIM_plot,Dis_GM{2},'--','LineWidth',2)
              axis([min(Tgt_distro_output.TgtIM_plot) max(Tgt_distro_output.TgtIM_plot) 0 1])
              set(gca,'XTickLabel',{'10^{-1}s' '1s' '10^{1}s' 'Ds5-95%'})
              xlabel('T (s) and Ds5-95%')
              ylabel('Standard deviation of lnIM')
              title('Response spectrum standard deviation match by MCMCGMS method 2 ((g) for SA and (s) for Ds5-95%)');
              set(gca,'fontsize',9,'fontname','Times New Roman')
              hleg=legend('Target standard deviation of lnIM','Standard deviation of selected lnIM');
              set(hleg,'fontsize',9,'fontweight','normal','fontname','Times New Roman');
              legend boxoff              
end               
