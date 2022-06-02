function out = mean_sigma_var(val, sigma)

if nargin < 2 || isempty(sigma), sigma = 1; end

out = mean(val) + sigma * var(val); 
