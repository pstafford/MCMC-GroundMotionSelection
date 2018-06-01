function [data_after_preslection,GM] = database_input_average_SA_DS595(Input_periods,Binallow)
% This function is used to preselect records based on bin size of casual parameters (i.e. magnitude and distance).
% Output
% data_after_preslection                    : A structure with parameters that specify output of database after preselection
%.SA                                        : Spectral acceleration
%.Ds                                        : 5-95% significant duration

% GM                                        : A structure with parameters that specify output of records ID and casual parameters after preselection
%.GMID_after_preselection                   : Record ID
%.MR_after_preselection                     : Casual parameters (i.e. magnitude and distance)

% Periods of interesting spectral accelerations
PerTgt = Input_periods;  
% Load the database
load rec_selection_meta_data
% Matrix to save Target periods
Tgperiod = zeros(1,length(PerTgt));
% Ground motion records are scaled to have the same conditional SA (find relative period)
for i = 1:length(PerTgt)
    [~,Tgperiod(i)] = min (abs(Periods - PerTgt(i)));
end
% Read prespective GMID 
GMIDread = xlsread('ID.xlsx');
GM.GMID = GMIDread(2:2:end,:);
% Read magnitude and distance 
Magnitude_Distance_read = xlsread('Magnitude_Distance.xlsx');
GM.MR = Magnitude_Distance_read(2:2:end,:);
% Geometric mean of DS
data_Dsread = xlsread('Ds595data.xlsx');
data_Ds1 = data_Dsread(1:2:end-1,:);
data_Ds2 = data_Dsread(2:2:end,:);
data_Ds_before_preselection = sqrt(data_Ds1.*data_Ds2);

% Matrix to save Sa after arranging and corresponding file name
Sa_1arrange = zeros(length(GM.GMID),length(Sa_1(1,:)));
Sa_2arrange = zeros(length(GM.GMID),length(Sa_1(1,:)));

% Rearrange database based on GMID
for i = 1:length(GM.GMID)
   Sa_1arrange(i,:) = Sa_1(GM.GMID(i),:);
   Sa_2arrange(i,:) = Sa_2(GM.GMID(i),:);
end
% Geometric mean of SA
Geo_SA = sqrt(Sa_1arrange.*Sa_2arrange);
database_SA_before_preselection = Geo_SA(:,Tgperiod);

% Pre selection based on bin size
recValidMag =  GM.MR(:,1) > Binallow.M(1)  & GM.MR(:,1) < Binallow.M(2);
recValidDist = GM.MR(:,2)> Binallow.D(1)    & GM.MR(:,2) < Binallow.D(2);

% After pre selection
allowedIndex_data = find(recValidMag & recValidDist); 
data_after_preslection.SA = database_SA_before_preselection(allowedIndex_data,:);
data_after_preslection.Ds = data_Ds_before_preselection(allowedIndex_data);
GM.GMID_after_preselection = GM.GMID(allowedIndex_data);
GM.MR_after_preselection = GM.MR(allowedIndex_data,:);
