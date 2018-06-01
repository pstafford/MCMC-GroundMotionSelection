function [rho]=Bradley2011_SA_Ds595_correlation(inputs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Purpose: To provide the empirical correlation of SA (T=0.01-10s) and Ds595
%based on Bradley (2011)

%Reference: Bradley, B. A., 2011. Correlation of Significant Duration with
%Amplitude and Cumulative Intensity Measures and Its Use in Ground Motion 
%Selection, Journal of Earthquake Engineering,  15, 809-832.

% Input variables:
% inputs - a structure with the following elements
%       .Ti or .Tj - the vibration period of the SA ordinate

% Output variables:
% rho - the predicted correlation coeff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start of code

T=inputs;

if (T<0.01)|(T>10)
    fprintf('Error: period range outside of applicability \n');
    return;
end

%coefficients from Table 3 of Bradley (2011)
a=[-0.41 -0.41 -0.38 -0.35 -0.02 0.23 0.02];
b=[ 0.01  0.04  0.08  0.26  1.40  6.0 10.0];

for n=2:length(a)
   if T<=b(n)
       rho=a(n-1)+log(T/b(n-1))/log(b(n)/b(n-1))*(a(n)-a(n-1));
       break;
   end
end

%end of script
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
