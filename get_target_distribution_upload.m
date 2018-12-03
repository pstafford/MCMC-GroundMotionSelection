function [Tgt_distro_output,Database_afterScaling]=get_target_distribution_upload(Input, s_input, ds_input, data_after_preslection, PerTgt_input_from_main, T1_from_main, max_Scalefactor,GM)
                                                                                  
% This function is used to get target distribution
% This code is only used for a test case (only SA and 5-95DS are considered)
% This code can be easily generalized to account for different intensity measures
% The Campbell and Bozorgnia (2008) Ground Motion Model is used in this code for SA. 
% The Bommer, Stafford and Alarcón (2009) Ground Motion Model (GMM) is used in this code for Ds595.
% The Baker and Jayaram (2008) Correlation function is used to calculation correlation of SA
% The Bradley (2011) Correlation function is used to calculation correlation between SA and Ds595

% Tgt_distro_output               : A structure with parameters that specify output of target distribution
% .TargetMean_plot                : Target mean for plot
% .TargetCov_plot                 : Target Covariance for plot 
% .TgtIM_plot                     : Target intensity vector for plot
% .TargetCov_PD                   : Target Covariance which is positive definite 
% .Max_change1_Corr               : Maximum changes on correlation matrix due to nearest positive definite matrix algorithm applied

% Database_afterScaling           : A structure with parameters that specify output of database after scaling
% .GM_plot                        : Database after scaling for plot
% .GMID                           : Record ID corresponding to database after preselection on limitation of scaling 
% .GMMR                           : Magnitude and distance corresponding to database after preselection on limitation of scaling
% .SF                             : Scale factor corresponding to database after preselection on limitation of scaling


% Inputs
Mmagnitude = Input.Mmagnitude;      % Magnitude
Rrupture = Input.Rrupture;          % Distance (km)
Rjb = Input.Rrupture;               % Assumed same for example

% Inputs for ground motion prediction equation (right lateral strike slip)
deltainput = Input.deltainput;      % Average dip of the rupture place (degree)
lambdainput = Input.lambdainput;    % Rake angle (degree)
Vs30input = Input.Vs30input;        % Shear wave velocity averaged over top 30 m (m/s)
Zvsinput = Input.Zvsinput;          % Depth to the 2.5 km/s shear-wave velocity horizon (km)
arbinput = Input.arbinput;          % 1 for arbitrary component sigma, 0 for average component sigma
Ztor = Input.Ztor;                  % Depth to the top of coseismic rupture (km)
Epsilon = Input.Epsilon;            % Epsilon is assumed 2

s = s_input;                        % Short name of SA
ds = ds_input;                      % Short name of Ds595
T1 = T1_from_main;                  % conditional period
PerTgt = PerTgt_input_from_main;    % interesting period range
% Input end

% Target IM vector
TgtIM = [PerTgt,100]; % 100 is set here for plotting Ds5-95; 
% Identify conditional IM
Index_T1 = find(PerTgt == T1);

% Matrix to save intensity measures' name
IMiNames = cell(1,length(TgtIM));
% Build intensity names vector
for i = 1:length(TgtIM)
    if i <= length(PerTgt)
	    IMiNames(i) = {s};
    else
      IMiNames(i) = {ds};
    end
end
% Matrix to save mean and standard deviation
MeanIM = zeros(1,length(TgtIM));
SdIM = zeros(1,length(TgtIM));

% Calculation of mean and standard deviation on intensity vector
for i = 1:length(TgtIM)
  if strcmp(IMiNames(i),'SA') == 1
	    [MeanIM(i),SdIM(i)] = CampbellBozorgnia_2008_SA (Mmagnitude, PerTgt(i), Rrupture, Rjb, Ztor, deltainput, lambdainput, Vs30input, Zvsinput, arbinput);
  else
        [MeanIM(i),SdIM(i)] = BSA_09_Ds(Mmagnitude,Rrupture,Vs30input,Ztor,arbinput);
  end
end

% Covariance
% Matrix to save rho and unconditional covariance
rho = zeros(length(TgtIM));
Uncondition_Cov = zeros(length(TgtIM));

% Unconditional covariance
for i = 1:length(TgtIM)
  for j = 1:length(TgtIM)
      
      if strcmp(IMiNames(i),'SA')==1
           if strcmp(IMiNames(j),'Ds595')==1
             rho(i,j)=Bradley2011_SA_Ds595_correlation(TgtIM(i));
          else
             rho(i,j)=Baker_Jayaram_correlation(TgtIM(i),TgtIM(j));
           end
      else
           if strcmp(IMiNames(j),'Ds595')==1
             rho(i,j)=1;
          else
             rho(i,j)=Bradley2011_SA_Ds595_correlation(TgtIM(j)); 
          end
       end
       Uncondition_Cov(i,j)=rho(i,j)*SdIM(i)*SdIM(j);
   end
end

% Matrix to save correlation and Covariance between IMi and IMj
rho_T1 = zeros(length(TgtIM),1);
Cov_T1_other = zeros(length(TgtIM),1);

% Covariance between T1 and other intensity measures
for i = 1:length(TgtIM)
 if strcmp(IMiNames(i),'SA') == 1
   rho_T1(i) = Baker_Jayaram_correlation(PerTgt(i),T1);
 else
   rho_T1(i) = Bradley2011_SA_Ds595_correlation(T1);
 end
   Cov_T1_other(i,:) = rho_T1(i)*SdIM(i)*SdIM(TgtIM == T1);
end

Cov_T1_other_save = Cov_T1_other*(1/SdIM(TgtIM == T1)^2)*(Cov_T1_other)'; 

% Conditional correlation matrix
% Parameters used in covariance modify 
TgtIM_PD_modify = TgtIM;
rho_PD_modify = rho;
rho_T1_PD_modify = rho_T1;
SdIM_PD_modify = SdIM;

% Delete conditional IM content
TgtIM_PD_modify(Index_T1) = [];
rho_PD_modify(Index_T1,:) = [];
rho_PD_modify(:,Index_T1) = [];
rho_T1_PD_modify(Index_T1) = [];
SdIM_PD_modify(Index_T1) = [];

% Matrix to save 
corr_codition = zeros(length(TgtIM_PD_modify));
for i = 1:length(TgtIM_PD_modify)
 for k = 1:length(TgtIM_PD_modify)
     corr_codition(i,k) = (rho_PD_modify(i,k) - rho_T1_PD_modify(i)*rho_T1_PD_modify(k))/(sqrt(1 - rho_T1_PD_modify(i)^2)*sqrt(1 - rho_T1_PD_modify(k)^2));
 end
end

[~,p] = chol(corr_codition);

if p ~= 0 
 b = ones(length(rho_PD_modify(1,:)),1); % Initialize parameters for function call
 tau = 1.0e-2;                           % Initialize parameters for function call
 tol = 1.0e-2;                          % Initialize parameters for function call
% Identify the nearest positive definite correlation matrix
 corr_codition_modify = nearestCorrelationMatrix(corr_codition, b, tau, tol);
% Matrix used to save covariance obtained with and without PD modifying correlation matrix 
 TargetCov = zeros(length(TgtIM_PD_modify));
 TargetCov_unchanged = zeros(length(TgtIM_PD_modify));

 for i = 1:length(TgtIM_PD_modify)
   for k = 1:length(TgtIM_PD_modify)
     TargetCov(i,k) = (sqrt(1 - rho_T1_PD_modify(i)^2)*sqrt(1 - rho_T1_PD_modify(k)^2))*(SdIM_PD_modify(i)*SdIM_PD_modify(k))*corr_codition_modify(i,k);
     TargetCov_unchanged(i,k) = (sqrt(1 - rho_T1_PD_modify(i)^2)*sqrt(1 - rho_T1_PD_modify(k)^2))*(SdIM_PD_modify(i)*SdIM_PD_modify(k))*corr_codition(i,k);
   end
 end

  Tgt_distro_output.Max_change1_Corr = max(abs(TargetCov-TargetCov_unchanged));
  TargetCov_PDmodify = TargetCov;
% set covariance zero at IMj   
  TargetCov_PDmodify = [TargetCov_PDmodify(1:Index_T1-1,:);zeros(1,length(TargetCov));TargetCov_PDmodify(Index_T1:end,:)];
  TargetCov_PDmodify = [TargetCov_PDmodify(:,1:Index_T1-1),zeros(length(TargetCov)+1,1),TargetCov_PDmodify(:,Index_T1:end)];
  
% Target Covariance obtained based Eq 4 shown in the article which is equal to TargetCov_unchanged used here. 
  TargetCov_Eq_4 = Uncondition_Cov - Cov_T1_other_save;

else
  TargetCov_Eq_4 = Uncondition_Cov - Cov_T1_other_save;
  TargetCov_PDmodify = TargetCov_Eq_4; 
  Tgt_distro_output.Max_change1_Corr = 0;
end  

% Matrix to save target mean
TargetMean = zeros(1,length(TgtIM));
% Conditional mean 
for i = 1:length(TgtIM)
    TargetMean(i) = log(MeanIM(i))+rho_T1(i)*Epsilon*SdIM(i);
end
% Save the location of conditional intensity in the intensity measure vector 
lnSa1 = TargetMean(PerTgt == T1);
% Target distribution obtained

% Matrix to save database after scaling
database_SA_afterScaling = zeros([length(data_after_preslection.SA(:,1)) length(PerTgt)]);

% Matrix to save scale factor
Scalefactor = zeros(1,length(data_after_preslection.SA(:,1)));

% Database after scaling 
for i = 1:length(data_after_preslection.SA(:,1))
    Scalefactor(i) = exp(lnSa1)/data_after_preslection.SA(i,Index_T1);
    database_SA_afterScaling(i,:) = log(data_after_preslection.SA(i,:)*Scalefactor(i));% log database (SA)
end

% Remove records due to maximum scale factor limitation (ID MR and scale factor matrix modified)
remove_ID_data = find(Scalefactor >= max_Scalefactor);
Scalefactor(remove_ID_data) = [];
database_afterScaling = [database_SA_afterScaling,log(data_after_preslection.Ds)]; % log database (DS) and scale factor for duration is one
database_afterScaling(remove_ID_data,:) = [];
GM.GMID_after_preselection(remove_ID_data) = [];
GM.MR_after_preselection(remove_ID_data,:) = [];

% Save scale factor, record ID and magnitude-distance matrix for output
Database_afterScaling.SF = Scalefactor;
Database_afterScaling.GMID = GM.GMID_after_preselection;
Database_afterScaling.GMMR = GM.MR_after_preselection;
% Save matrixes for plotting or selecting GMs before deleting values corresponding to conditional intensity measure
Database_afterScaling.GM_plot = database_afterScaling;
Tgt_distro_output.TgtIM_plot = TgtIM;
Tgt_distro_output.TargetMean_plot = TargetMean;
Tgt_distro_output.TargetCov_plot = TargetCov_Eq_4;
Tgt_distro_output.TargetCov_PD = TargetCov_PDmodify;
