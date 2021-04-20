function xddot = cip_dyn(xu)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CART INVERTED PENDULUM DYNAMICAL EQUATIONS
%
% Brent Wallace
%
% 3/18/2021
%
% State space equations for cart inverted pendulum model (with
% translational damping, no rotational damping, pendulum mass concentrated
% at end). See Ogata, K. "Modern Conrol Engineering" 3rd ed. pp. 106.
%
% *************************************************************************
%                               CALL SYNTAX:
% *************************************************************************
%
% xddot = cip_dyn(xu)
%
% *************************************************************************
%                             	INPUTS: 
% *************************************************************************
%                      
%   xu                   : 5-dimensonal state vector, the first 4
%                           dimensions of which consist of the state:
%                               x = [ x x'  theta  theta' ]^T
%                           the last entry is the input force u.
%
%
% *************************************************************************
%                             	OUTPUTS: 
% *************************************************************************
% 
%   xdot                : 2-dimensonal vector of second order derivatives:
%                           xddot = [ x''  theta'' ]^T
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


% *************************************************************************
%
% UNPACK STATE VECTOR, CONTROL
%
%   xu = [ x x'  theta  theta'  u ]^T
% 

% x = xu(1);                  % Cart x - position. (NOT NEEDED)
xdt = xu(2);                % Cart x - velocity.
th = xu(3);                 % Pendulum angle theta.
thdt = xu(4);               % pendulum angular velocity theta'.
u = xu(5);                  % Control signal u (input force).


% *************************************************************************
%
% DYNAMICAL EQUATIONS
% 

% "Effective" mass seen in the xddot equation.
meff = M + m * (1 - cos(th)^2);      

% x''
xddt = ...
    ( u + m * l * sin(th) * thdt^2 - m * g * sin(th) * cos(th) ...
    - b * xdt) / meff;

% thata''
thddt = g / l * sin(th) - 1 / l * xddt * cos(th);


% *************************************************************************
%
% PACK OUTPUT VECTOR
%     
%   xddot = [ x''  theta'' ]^T
%

xddot = [   xddt
            thddt   ];

