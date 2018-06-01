function [MCMCGMS_2_output]=MCMCGMS_method2_for_case_study_SA_PLUS_DS595_upload(Input_MCMCGMS_2,Tgt_distro_output,Database_afterScaling,Index_T1, Input)

% This function is used to perform MCMCGMS method 2

% Output
% MCMCGMS_2_output                : A structure with parameters that specify output of MCMCGMS method 2

% .Index_GM_after_OP              : Record ID of selected records after the optimization applied

% Input
% Index_T1  : the location of conditional period
nGM = Input.nGM;                                       % Number of ground motions needed to be selected
nloop = Input_MCMCGMS_2.nloop;                         % Number of loop for optimization
TargetMean = Tgt_distro_output.TargetMean_plot;        % Target mean
TargetCov = Tgt_distro_output.TargetCov_PD;            % Target Covariance
database_afterScaling = Database_afterScaling.GM_plot; % Database after scaling

%%======================================================
% Modify distribution begin
%%======================================================

% Modify database
database_afterScaling(:,Index_T1) = [];
% Modify target distribution
TargetMean(Index_T1) = [];
TargetCov(Index_T1,:) = [];
TargetCov(:,Index_T1) = [];

%%====================================================
% Modify end
%%====================================================

% PDF pre-selection
if Input_MCMCGMS_2.pre_select_input==1
% Save databased before PDF pre-selection   
  database_before_reduce=database_afterScaling;
% Calculation of database pdf value  
  pdf_value=mvnpdf(database_afterScaling,TargetMean,TargetCov);
% Identify records whose pdf values larger than median pdf value 
  database_afterScaling=database_afterScaling(pdf_value>median(pdf_value),:);
end

% Test whether enough records can be selected
if length(database_afterScaling) < nGM 
  error('Error: no enough records, please reduce bounds casual parameters/ scale factor / preselection ');
end 

% Selection algorithm 
% (this selection algorithm can be replaced by simply selecting a set of records which have the first nGM highest PDF values)
% Selected records ID
Index_GM_Final = cell(1,length(database_afterScaling(:,1)));
% Selected records
GM_Final_Selection = cell(1,length(database_afterScaling(:,1)));

iter=1; 
while iter <= length(database_afterScaling(:,1))
% Randomly select the ground motion record
    IM_SA_Random = randi(size(database_afterScaling,1));
    GM_first = database_afterScaling(IM_SA_Random,:);
% Calculate the ratio between the random record and other records in the database
% Matrix to save ratio
    Rratio = zeros(1,length(database_afterScaling(:,1)));
% Matrix to save record ID   
    Rratio_Index_GM=zeros(1,length(database_afterScaling(:,1)));

% Calculate ratio of each record in the database to randomly selected record    
   for i=1:size(database_afterScaling,1)
      Rratio(i) = mvnpdf(database_afterScaling(i,:),TargetMean,TargetCov)/mvnpdf(GM_first,TargetMean,TargetCov);
      Rratio_Index_GM(i) = i;
   end
% Sort ratio
   [DescendRatio,index_Ratio] = sort(Rratio,'descend');
% Find the target prospective GM
   Index_GM = find(DescendRatio>1 & DescendRatio~=inf);
   Index_trans = index_Ratio(Index_GM);
   Index_GM_Prospective = Rratio_Index_GM(Index_trans);

% Check the number of ground motion selected
   if length(Index_GM) >= nGM
% Select first nGM 
     Index_GM_Final{iter} = Index_GM_Prospective(1:nGM);  
     GM_Final_Selection{iter} = database_afterScaling(Index_GM_Final{iter},:);
     break  
   end
   iter=iter+1;
end

% If nGM can not be selected by using loop shown above (e.g. 50 records seleced from the database of 50 records)
% First nGM with highest PDF values are selected
if iter == length(database_afterScaling(:,1))+1 
   pdf_value_1 = mvnpdf(database_afterScaling,TargetMean,TargetCov);
   [~,pdf_value_1_order] = sort(pdf_value_1,'descend');
   Index_GM_Final{iter-1} = pdf_value_1_order(1:nGM)';  
   GM_Final_Selection{iter-1} = database_afterScaling(Index_GM_Final{iter-1},:);
end   

% Remove empty matrix
Index_emptycell = cellfun(@isempty,GM_Final_Selection);
% Selected records
GM_Final_Selection(Index_emptycell) = [];

Index_emptycell_GMindex = cellfun(@isempty,Index_GM_Final);
Index_GM_Final(Index_emptycell_GMindex) = [];

% Greedy optimization 
% The input of greedy optimization
% Selected records
GM_optim_OP = GM_Final_Selection{1};
% Selected records ID
Index_optim = Index_GM_Final{1}';
% Weight used for mean
  weight(1) = Input_MCMCGMS_2.mean;
% Weight used for standard deviation
  weight(2) = Input_MCMCGMS_2.std;
% Weight used for standard deviation
  weight(3) = Input_MCMCGMS_2.ske;
% Weight factor used for SA
  weight_Sa = Input_MCMCGMS_2.W_SA;
% Weight factor used for DS
  weight_Ds = Input_MCMCGMS_2.W_Ds_595;
% Target conditional standard deviation
  Std_Target = sqrt(diag(TargetCov));

numIMs = length(TargetMean);
% Weight factor used for each IMi
if weight_Sa == weight_Ds && weight_Ds == 1.0
   fprintf('No weight factors are explicitly applied to intensity measures \n '); 
   weightsIM = ones(1,numIMs);
else
   fprintf('Normalized weight factors are applied to intensity measures \n ');   
   weightsIM = [ repmat(weight_Sa/(numIMs-1),1,numIMs-1), weight_Ds ];
end

disp('Please wait...Optimization used in MCMCGMS method 2');

% loop begin 
for n = 1:nloop

  for i = 1:nGM
      
% The initial dev is set as a extremely big value  	
    mindevInitial = 1000000;   
% Remove one prospective ground motion selected
    GM_optim_OP(i,:) = [];
    Index_optim(i,:) = [];

    for j = 1:length(database_afterScaling(:,1))
        GmforgreedyOptest = [GM_optim_OP;database_afterScaling(j,:)];
     
        if Input_MCMCGMS_2.ske_dev_check == 1
            [mu, sig, skew] = PJSmoments(GmforgreedyOptest);
            devSkew = weightsIM .* ( skew ).^2;
        else 
            [mu, sig] = PJSmoments(GmforgreedyOptest);
            devSkew = zeros(1,length(TargetMean));
        end
        
% Mean and standard deviation
        devMeanGmtest = weightsIM .* ( mu - TargetMean ).^2;
        devSigmaGmtest = weightsIM .* ( sig - Std_Target' ).^2;
% Dev check
      mindevtest = weight(1)*sum(devMeanGmtest)+weight(2)*sum(devSigmaGmtest)+weight(3)*sum(devSkew);

% Check whether the distribution match of selected record set is improved
     if  ((mindevtest < mindevInitial) &&  ~any(Index_optim == j))
          mindevInitial = mindevtest;
          Indeximpro = j;
     end

   end
% Add new a record and its record ID   
   GM_optim_OP = [GM_optim_OP(1:i-1,:);database_afterScaling(Indeximpro,:);GM_optim_OP(i:end,:)];
   Index_optim = [Index_optim(1:i-1);Indeximpro;Index_optim(i:end)];
 
    frac_progress = ((n-1)*nGM + i)/(nloop*nGM);
    rfrac_progress = round(frac_progress*100);
    if rem(rfrac_progress,10) == 0
        disp([num2str(rfrac_progress),'% done']);
    end
 end
end


% Selected records ID set (after optimization applied) with or without PDF preselection for output
if Input_MCMCGMS_2.pre_select_input == 1
  [~,~,MCMCGMS_2_output.Index_GM_after_OP] = intersect(database_afterScaling(Index_optim,:),database_before_reduce,'rows');
else
  MCMCGMS_2_output.Index_GM_after_OP = Index_optim;
end

