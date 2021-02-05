function table = tangent_cone_table(M1, g, max_theta)

if nargin < 3
    
    [~, max_tbm] = theta_beta_mach_curves(M1);
    max_theta = max_tbm(2);
end

theta = (0.001:0.001:max_theta)';
tau = asin(sin(theta) .* (((g + 1)/2) + (1./((M1 * sin(theta)).^2))).^0.5);

i = numel(tau);
M = zeros(size(tau));

% Exact Taylor-Maccoll solution?
while i > 0
    
    [~,M(i),~] = solvecone(tau(i), M1, g);
    if M(i) == M1, break; end
    i = i - 1;
end

table = [tau, M];
table(1:i, :) = [];