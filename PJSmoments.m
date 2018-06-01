function [ mu, sig_unbiased, skew ] = PJSmoments( x, dim )
%PJSMOMENTS function to compute first three moments from a sample
%   This function essentially just computes the mean, standard deviation
%   and skewness using the same equations as Matlab's built in functions.
%   However, this function avoids the repetition of internal computations
%   and so is significantly faster when computing all three moments
%   simultaneously.
%
% Written by Peter J Stafford
% Last modified: June 1, 2018

if nargin < 2
    dim = 1;
end

N = length(x);
mu = mean(x,dim);
dxmu = x-mu;
dev2 = dxmu.*dxmu;
sig2 = sum(dev2, dim)/N;
sig_unbiased = sqrt( sig2 * N / (N-1) );
if nargout > 2
    skew = sum( dev2 .* dxmu, dim )/N ./ sig2.^(3/2);
end

end

