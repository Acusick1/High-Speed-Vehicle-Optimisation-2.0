function [theta_beta_mach, max_theta_beta_mach] = theta_beta_mach_curves(Minf, beta, gamma)
%THETA_BETA_MACH_CURVES get deflection angle tables from shockwave angle
%and freestream Mach number combinations
%   Inputs:
%   Minf - vector of increasing freestream Mach numbers
%   beta - vector of increasing shockwave angles (radians)
%   gamma - Specific heat constant of gas
%   
%   Outputs:
%   theta_beta_mach - Full theta beta mach table
%   max_theta_beta_mach - Table of maximum deflection angles which allow
%       weak shockwave at each freestream Mach number

if nargin < 1, Minf = 1:0.01:20; end
if nargin < 2, beta = 0:0.0001:pi/2; end
if nargin < 3, gamma = 1.4; end

% Ensure Mach is row vector and beta is column vector
Minf = Minf(:)';
beta = beta(:);

% Get theta for all combinations of Mach and shockwave angle
theta = atan(2 * cot(beta) .* ((Minf.^2 .* (sin(beta).^2)-1) ./ ...
    (Minf.^2 .* (gamma + cos(2 * beta)) + 2)));

% Get maximum deflection angle for each Mach number, along with location
[max_theta, id] = max(theta, [], 1);
max_beta = beta(id);

theta_mach = [Minf; theta]; 
theta_beta_mach = [[NaN; beta], theta_mach];

max_theta_beta_mach = [Minf', max_theta', max_beta];

%% Plot proof
% figure
% hold on
% for i = 1:100:size(theta, 2)
%     
%     plot(rad2deg(theta(:,i)), rad2deg(beta))
% end
% axis([0 45 0 90])
% hold off

%% Save to file
% save('thetaBetaCurves','thetaBetaM','max_thetaBetaM')