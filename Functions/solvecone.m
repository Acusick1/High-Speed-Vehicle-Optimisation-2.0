function [thetac,Mc,sol]=solvecone(thetas,Minf,gamma)
% Solves the right circular cone at zero angle of attack in supersonic flow
% thetas - shock angle [degrees]
% Minf - upstream Mach number
% gamma - ratio of specific heats

% Minf = 4.63;
% thetas = 20;
% gamma = 1.4;

% Convert to radians
thetasr=thetas;

if (thetasr<=asin(1/Minf))
    thetac=0.0; Mc=Minf; sol =[];
    return;
end

% Calculate initial flow deflection (delta) just after shock...
delta=thetabetam(thetasr,Minf,gamma);
Mn1=Minf*sin(thetasr);
Mn2=sqrt((Mn1.^2+(2/(gamma-1)))./(2*gamma/(gamma-1)*Mn1.^2-1));
M2=Mn2./(sin(thetasr-delta));

% All values are non-dimensionalized!
% Calculate the non-dimensional velocity just after the oblique shock
V0=(1+2./((gamma-1)*M2.^2)).^(-0.5);
% Calculate velocity components in spherical coordinates
Vr0=V0.*cos(thetasr-delta);
Vtheta0=-V0.*sin(thetasr-delta);
% Calculate initial values for derivatives
dVr0=Vtheta0;
% Assemble initial value for ODE solver
y0=[Vr0;dVr0];

% Set up ODE solver
% Quit integration when vtheta=0
% See: coneevent.m
% Set relative error tolerance small enough to handle low M
options=odeset('Events',@coneevent,'RelTol',1e-5);
% Solve by marching solution away from shock until either 0 degrees or flow
% flow tangency reached as indicated by y(2)==0.
% See cone.m
[sol]=ode15s(@cone,[thetasr 1e-10],y0,options,gamma);
% Check if we have a solution, as ode15s may not converge for some values.
[n,m]=size(sol.ye);
thetac=0.0;
Mc=Minf;
% If ODE solver worked, calculate the angle and Mach number at the cone.
if (n>0 && m>0 && abs(sol.ye(2))<1e-10)
    thetac=sol.xe*180.0/pi;
    Vc2=sol.ye(1)^2+sol.ye(2)^2;
    Mc=((1.0/Vc2-1)*(gamma-1)/2)^-0.5;
end
end

function [dy]=cone(theta,y,gamma)
% y is a vector containing vr, vr’
% Governing equations are continuity, irrotationality, & Euler’s equation.
dy=zeros(2,1);

dy(1)=y(2);
dy(2)=(y(2)^2*y(1)-(gamma-1)/2*(1-y(1)^2-y(2)^2)*(2*y(1)+y(2)*cot(theta)))...
    /((gamma-1)/2*(1-y(1)^2-y(2)^2)-y(2)^2);
end

function [value,isterminal,direction]=coneevent(~,y,~)
% Check cone solution for point where vtheta=0
% theta - current angle
% y - current solution vector
% gamma - ratio of specific heats

value=zeros(2,1);
isterminal=zeros(2,1);
direction=zeros(2,1);

%Quit if Vr goes negative (which shouldn’t happen!)
value(1)=1.0;
if (y(1)<0.0)
    value(1)=0.0;
end
isterminal(1)=1;
direction(1)=0;

%Quit if Vtheta goes positive (which occurs at the wall)
value(2)=1.0;
if (y(2)>0.0)
    value(2)=0.0;
end
isterminal(2)=1;
direction(2)=0;
end

function [theta]=thetabetam(beta,M,gamma)
% Return theta for beta-theta-M relationship for oblique shocks
% beta - shock angle in radians
% M - upstream Mach number
% gamma - ratio of specific heat
%Cut off at Mach wave angle
if (beta<=asin(1/M))
    theta=0;
    return
end

theta=atan(2*cot(beta).*((M*sin(beta)).^2-1)./(M^2*(gamma+cos(2*beta))+2));
end