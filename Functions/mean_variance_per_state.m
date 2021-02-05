function out = mean_variance_per_state(state_uncert, id, sigma, obj)
%% Mean variance multi-flightstate cost function

if nargin < 3 || isempty(sigma), sigma = 6; end
if nargin < 4 || isempty(obj), obj = "Cd"; end

for i = size(id, 1):-1:1
    
    state(:,i) = [state_uncert(id(i,:)).(obj)]; 
end

out = mean(state) + sigma * var(state);

end