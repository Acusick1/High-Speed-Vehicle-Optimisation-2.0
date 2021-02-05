function [theta_beta_mach, max_theta_beta_mach] = theta_beta_mach_curves(mach, beta, g)

if nargin < 1, mach = 1:0.01:20; end
if nargin < 2, beta = (0:0.0001:pi/2)'; end
if nargin < 3, g = 1.4; end

theta = atan(2 * cot(beta) .* ((mach.^2 .* (sin(beta).^2)-1) ./ ...
    (mach.^2 .* (g + cos(2 * beta)) + 2)));

[maxTheta, id] = max(theta, [], 1);
maxBeta = beta(id);
beta = [NaN; beta];
theta_mach = [mach; theta]; 
theta_beta_mach = [beta, theta_mach];

max_theta_beta_mach = [mach' maxTheta' maxBeta];

% figure
% hold on
% for i = 1:100:size(theta, 2)
%     
%     plot(rad2deg(theta(:,i)), rad2deg(beta))
% end
% axis([0 45 0 90])
% hold off

% save('thetaBetaCurves','thetaBetaM','maxThetaBetaM')