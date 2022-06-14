function [match_del, match_mach] = matching_point(gamma, Pinf_P0)
%% TODO: match_del not implemented anywhere, should be used as a switch with Aerodynamics?
%MATCHING_POINT calculates point to switch from Newtonian to Prandtl-Meyer
%   Find point (inclination angle and Mach number) at which Prandtl-Meyer 
%   and Newtonian calculated pressures are equal

fun = @(Mq, gamma)(((2./(2+(gamma-1).*(Mq.^2))).^(gamma/(gamma-1)))...
    .*(1-(((gamma^2).*(Mq.^4).*((2./(2+(gamma-1).*(Mq.^2))).^(gamma/(gamma-1))))...
    ./(4*((Mq.^2)-1).*(1-((2./(2+(gamma-1).*(Mq.^2))).^(gamma/(gamma-1))))))));

wrapper_fun = @(Mq) fun(Mq, gamma);

M = bisection(wrapper_fun, 1.1 , 1.75, Pinf_P0);

if ~isnan(M)
    Q = ((2./(2+(gamma-1).*(M.^2))).^(gamma/(gamma-1)));
    P = Q.*(1-(((gamma^2).*(M.^4).*Q)./(4*((M.^2)-1).*(1-Q))));
    
    match_mach = M;
    match_del = asin(((Q-P)/(1-P))^0.5);
else
    match_mach = inf;
    match_del = inf;
end

%% OLD METHOD from super-hypersonic arbitrary-body program 1970
% % Initial guesses
% Mq = [1.35,1.75];
% Q = ((2./(2+(gamma-1).*(Mq.^2))).^(gamma/(gamma-1)));
% Pc = Q.*(1-(((gamma^2).*(Mq.^4).*Q)./(4*((Mq.^2)-1).*(1-Q))));
% 
% [~,closest] = min(abs(Pc - Pinf_P0));
% 
% it = zeros(1000,1);
% i = 2;
% 
% it(i) = Pc(closest);
% 
% % Iterate until convergence
% while it(i) ~= it(i-1)
%     
%     i = i+1;
%     
%     newMq = Mq(1)+(Pinf_P0-Pc(1)).*((Mq(2)-Mq(1))./(Pc(2)-Pc(1)));
%     
%     Q = (2/(2+(gamma-1)*(newMq.^2)))^(gamma/(gamma-1));
%     newPc = Q*(1-(((gamma^2)*(newMq^4)*Q)/(4*((newMq^2)-1)*(1-Q))));
%     
%     if newPc < Pinf_P0
%         
%         Mq(1) = newMq;
%         Pc(1) = newPc;
%     else
%         Mq(2) = newMq;
%         Pc(2) = newPc;
%     end
%     
%     [~,closest] = min(abs(Pc - Pinf_P0));
%     
%     it(i) = Pc(closest);
%     
% end
% 
% P = Pc(closest);
% 
% matchdel = asin(((Q-P)/(1-P))^0.5);
% matchMach = Mq(closest);