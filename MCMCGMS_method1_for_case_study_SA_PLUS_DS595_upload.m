function [MCMCGMS_1_output]=MCMCGMS_method1_for_case_study_SA_PLUS_DS595_upload(Input_MCMCGMS_1,Tgt_distro_output,Database_afterScaling, Index_T1, Input)

% This function is used to perform MCMCGMS method 1

% Output
% MCMCGMS_1_output                : A structure with parameters that specify output of MCMCGMS method 1
% .Index_optim                    : Record ID of selected records

% Input

%Index_T1 : the location of conditional period

nGM = Input.nGM;                                       % Number of ground motions needed to be selected
niter = Input_MCMCGMS_1.niter;                         % Number of iteration
nRep = Input_MCMCGMS_1.nRep;                           % Number of replication used in energy check
TargetMean = Tgt_distro_output.TargetMean_plot;        % Target mean
TargetCov = Tgt_distro_output.TargetCov_PD;            % Target Covariance
TgtIM = Tgt_distro_output.TgtIM_plot;                  % Target intensity measure vector
database_afterScaling = Database_afterScaling.GM_plot; % Database after scaling

%%======================================================
% Modify target distribution begin
%%======================================================

% Modify database
database_afterScaling(:,Index_T1) = [];

% Modify target distribution
TgtIM(Index_T1) = [];
TargetMean(Index_T1) = [];
TargetCov(Index_T1,:) = [];
TargetCov(:,Index_T1) = [];

numIMs = length(TargetMean);
% Weight factor used for each IMi
if Input_MCMCGMS_1.W_SA == Input_MCMCGMS_1.W_Ds_595 && Input_MCMCGMS_1.W_Ds_595 == 1.0
   fprintf('No weight factors are explicitly applied to intensity measures \n '); 
   weightsIM = ones(1,numIMs);
else
   fprintf('Normalized weight factors are applied to intensity measures \n ');   
   weightsIM = [ repmat(Input_MCMCGMS_1.W_SA/(numIMs-1),1,numIMs-1), Input_MCMCGMS_1.W_Ds_595 ];
end
   
%%====================================================
% Modify end
%%====================================================

% Maximum PDF used in the energy check
epsilon_maxPDF = mvnpdf(TargetMean,TargetMean,TargetCov);

% Dimensional used in the energy check
Dimention = length(TgtIM);

% PDF pre-selection
if Input_MCMCGMS_1.pre_select_input == 1
% Save databased before PDF pre-selection  
  database_before_reduce = database_afterScaling;
% Calculation of database pdf value
  pdf_value = mvnpdf(database_afterScaling,TargetMean,TargetCov);
% Identify records whose pdf values larger than median pdf value and the record with maximum pdf value 
  database_afterScaling = database_afterScaling(pdf_value>median(pdf_value),:);
  [~,first_ID] = sort(pdf_value(pdf_value>median(pdf_value)),'descend');
else
  pdf_value_all = mvnpdf(database_afterScaling,TargetMean,TargetCov);
  [~,first_ID] = sort(pdf_value_all,'descend');
end

% Test whether enough records can be selected
if length(database_afterScaling) == nGM
  if Input_MCMCGMS_1.pre_select_input == 1
    X = ['Only ', num2str(nGM), 'records are available in the database after PDF preselection applied'];
    disp(X);
    [~,~,MCMCGMS_1_output.Index_optim] = intersect(database_afterScaling(linspace(1,nGM,nGM),:),database_before_reduce,'rows');
    return
  else  
    X = ['Only ', num2str(nGM), 'records are available in the database'];
    disp(X);
    MCMCGMS_1_output.Index_optim = linspace(1,nGM,nGM);
    return
  end  
end   

if length(database_afterScaling) < nGM 
  error('Error: no enough records, please reduce bounds casual parameters/ scale factor / preselection ');
end 

% Matrix used to save current state and corresponding records ID
% It is difficult to preallocate array
IM = [];
GMSelected = [];

% Initial selected ground motion
IM_SA_Random = first_ID(1);

%%====================================================
% Random walk algorithm with bin size
%%====================================================

% Nosie matrix build
% Mean
MeanZ = zeros(length(TgtIM),1); 

% Diagonal matrix generated to save noise Cov matrix
Diagnal_Cov_one = eye(length(TargetCov));

% Conditional standard deviation
SDIM_SaT1 = sqrt(diag(TargetCov))';

for i = 1:length(TargetCov)
Diagnal_Cov_one(i,i) = SDIM_SaT1(i);
end
% Cov matrix
CovZ = Diagnal_Cov_one./Input_MCMCGMS_1.gamma; % noise Cov matrix

% Matrix to save the result of each energy check replication 
phi_nm_nRep = zeros(1,nRep);
% Matrix to save mean energy check of each selection replication
phi_nm = zeros(1,niter);

% Cell to save
% Selected records ID for each selection replication
GM_IDforEachloop = cell(1,niter);
% Selected records for each selection replication
FinalGm = cell(1,niter);

disp('Please wait...MCMCGMS method 1');
% Selection replication begin
for n = 1:niter
   if n == 1
      IM(1,:) = database_afterScaling(IM_SA_Random,:); % The first state ground motion
      GMSelected(1) = IM_SA_Random;                    % The ground motion ID for the first ground motion
   else
      IM(1,:) = IM(i,:);
      GMSelected(1) = GMSelected(i);
   end
    i=1;
    flag = 0;
    i_forGMID = 1; % Indicator for number of records have been selected

% Matrix to save GMID and the difference between candidate state and records in the database
    GM_final_Selection = zeros(nGM,1);
    err = zeros(1,length(database_afterScaling));

% Loop for one suite of ground motion selected begin    
   while sum(GM_final_Selection ~= 0) < nGM
      
      Z(1,:) = mvnrnd(MeanZ,CovZ);
% Precandidate calculate
      IM_Candi = IM(i,:)+Z(1,:);

% Select a record which is the closest to the generation
      err = sum( weightsIM .* (IM_Candi - database_afterScaling).^2, 2 ); 

      [~,indexGM] = min(err);
      IMselect = database_afterScaling(indexGM,:);
      DeltaIM = IMselect - IM(i,:);

% Bin check
       if sum(DeltaIM < 1.96*SDIM_SaT1) == length(DeltaIM) && sum(DeltaIM > -1.96*SDIM_SaT1) == length(DeltaIM)
% Acceptance ratio recheck
         Rratio = mvnpdf(IMselect,TargetMean,TargetCov)/mvnpdf(IM(i,:),TargetMean,TargetCov);
         AccepRatio = [1,Rratio];
         Alph = rand;
          if Alph <= min(AccepRatio)
             IM(i+1,:) = IMselect;
             GMSelected(i+1) = indexGM;
% Identify a selected record 
           if sum (GM_final_Selection == GMSelected(i+1)) == 0
            GM_final_Selection(i_forGMID) = GMSelected(i+1);
            i_forGMID = i_forGMID+1;
           end
          else
            IM(i+1,:) = IM(i,:);
            GMSelected(i+1) = GMSelected(i);
          end
       else
       	  IM(i+1,:) = IM(i,:);
       	  GMSelected(i+1) = GMSelected(i);
       end      
    i = i+1;
% paucity of records in the database (e.g. select 50 records from 100 records in the database)
        if i == 1000000
          flag = 1;
         disp('paucity of records in the database, MCMCGMS method 1 will take more time (To reduce time consuming, maybe a advance noise COV applied or reduce target number of records selected). Please use MCMCGMS method 2 instead')
         break
       end
   end
   if(flag == 1)
       break
   end      
% One suite of ground motions selected
  GM_IDforEachloop{n} = GM_final_Selection;                  % Save the records selected ID
  FinalGm{n} = database_afterScaling(GM_final_Selection,:);  % Save the records (without the consideration of IMj)

% The energy check begin
% The number of observations
  n_phi = length(FinalGm{n}(:,1));

% Energy check replication begin  
  for j = 1:nRep
% Monte Carlo applied to generate simulations based on the target multivariate distribution (without the consideration of IMj)
  Simulate_Target = (mvnrnd(TargetMean,TargetCov,10000));  
% The number of simulations
  m = length(Simulate_Target(:,1));
% The energy calculate 
   phi_nm_nRep(j) = W_rlog(FinalGm{n},n_phi,Simulate_Target,m,epsilon_maxPDF,Dimention,weightsIM);
  end

% Save mean value of energy check results
  phi_nm(n) = mean(phi_nm_nRep);

  disp([num2str(round(n/niter*100)) '% done']);

end

% Selected records set with minimum energy check result
[~, Index_phi] = min(phi_nm);

% Selected records ID set with or without PDF preselection for output
if Input_MCMCGMS_1.pre_select_input == 1
  [~,~,MCMCGMS_1_output.Index_optim] = intersect(database_afterScaling(GM_IDforEachloop{Index_phi},:),database_before_reduce,'rows');
else
  MCMCGMS_1_output.Index_optim = GM_IDforEachloop{Index_phi};
end






