function u = nln_hjb(states)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% NONLINEAR HJB CONTROL LAW
%
% Brent Wallace
%
% 3/18/2021
%
% *************************************************************************
%                               CALL SYNTAX:
% *************************************************************************
%
% u = nln_hjb(states)
%
% *************************************************************************
%                             	INPUTS: 
% *************************************************************************
%                      
%   states              	: 4-dimensonal state vector:
%                               states = [ x x'  theta  theta' ]^T
%
%
% *************************************************************************
%                             	OUTPUTS: 
% *************************************************************************
% 
%   u                   : Control signal (N).
%
% *************************************************************************
% *************************************************************************
% *************************************************************************


% *************************************************************************
%
% SYSTEM PARAMETERS
% 

M = 0.5;                % Cart mass (kg).
m = 0.2;                % Pendulum mass (kg).
l = 0.3;                % Pendulum length (m).
g = 9.8;                % Gravitational field constant (m/s^2).
b = 0.1;                % Cart velocity damping constant (N/(m/s)).

J = m * l^2;            % Pendulum moment of inertia (kg-m^2).
% J = 0.006;            % Pendulum moment of inertia (kg-m^2).

% *************************************************************************
%
% CONTROLLER PARAMETERS
% 

% c = 0.1, epsilon = 0.1 works out to ~ 45deg initial tilt. Settling time ~ 1
%   sec
% c = 0.1, epsilon = 1 works out to ~ 35deg initial tilt. Settling time ~ 5
%   sec.
% c = 0.01, epsilon = 1 works out to ~ 30deg initial tilt. Oscillatory.

c = 0.1;                  % System output transformation constant.
epsilon = 0.45;            % Controller gain.

% *************************************************************************
%
% UNPACK STATE VECTOR
%
%   states = [ x x'  theta  theta' ]^T
% 

x = states(1);                      % Cart x - position.
xdot = states(2);                   % Cart x - velocity.
theta = states(3);                  % Pendulum angle theta.
thetadot = states(4);               % pendulum angular velocity theta'.


% *************************************************************************
%
% CONTROL LAW
% 

z1 = x;
z2 = xdot;
znew = c*z1 + z2;

% n1 = theta - pi;
n1 = theta;
n2 = theta + thetadot + ( (m*l) / J ) * xdot * cos(theta);

ynum1 = -n1 + n2 + ((m*l*c)/J)*z1*cos(n1) + ((m*g*l)/J)*sin(n1);
ynum2 = ((m*l*c)/J)*z1*sin(n1) * (-n1 + n2 + ((m*l*c)/J)*z1*cos(n1));

yden1 = -((m*l)/J)*cos(n1) - (((m*l)/J)*sin(n1))*(-n1 + n2);
yden2 = -((2*(m^2)*(l^2)*c)/(J^2))*z1*cos(n1)*sin(n1);

ye = -2*(ynum1+ynum2)/(yden1+yden2);

u = -(1/epsilon)*(znew-ye);

