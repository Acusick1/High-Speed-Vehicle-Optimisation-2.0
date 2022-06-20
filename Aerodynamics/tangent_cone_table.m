function table = tangent_cone_table(M1, gamma, max_theta)
%TANGENT_CONE_TABLE Generate tangent cone table of 3D shock angles and
%downstream Mach number
%   Inputs:
%   M1 - Freestream Mach number
%   gamma - Gas specific heat ratio
%   max_theta - Upper boundary of table (default = max oblique shock) 

if nargin < 3
    
    [~, max_theta_beta_mach] = theta_beta_mach_curves(M1);
    max_theta = max_theta_beta_mach(2);
end

theta = (0.001:0.001:max_theta)';
tau = asin(sin(theta) .* (((gamma + 1)/2) + ...
           (1./((M1 * sin(theta)).^2))).^0.5);

% Exact Taylor-Maccoll solution?
for i = numel(tau):-1:1
    
    [~, M2(i), ~] = solvecone(tau(i), M1, gamma);
end

table = [tau, M2];