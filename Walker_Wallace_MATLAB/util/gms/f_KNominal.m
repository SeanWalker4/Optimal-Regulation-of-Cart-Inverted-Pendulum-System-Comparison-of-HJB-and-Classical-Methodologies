function [Ko, F, L] = f_KNominal(P_d)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% FORM NOMINAL STABOLIZING CONTROLLER
%
% Karan Puttannaiah
% Brent Wallace  
%
% 12/14/2020
%
% Forms a nominal stabolizing MBC controller Ko for the Youla et al.
% parameterization.
%
% ***** CALL SYNTAX:
%
% [Ko, F, L] = f_KNominal(P_d)
%
% ***** INPUTS:
%
%   P_d             : Design plant (state space object).
%
% ***** OUTPUTS:
%
%   Ko              : Nominal stabolizing controller (state space object).
%   F               : MBC full state feedback gain matrix.
%   L               : MBC observer gain matrix.
%
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

% ***********************
%
% PLANT STATE SPACE MATRICES AND DIMENSIONS
%
[Ap, Bp, Cp, Dp] = ssdata(P_d);
n_x = size(Ap,1); 
n_e = size(Cp,1); 
n_u = size(Bp,2);

% ***********************
%
% GET MBC FULL STATE FEEDBACK GAIN MATRIX
%
F = lqr(Ap, Bp, 1e0*eye(n_x), 1.5e1*eye(n_u));

% ***********************
%
% GET MBC OBSERVER GAIN MATRIX
%
L = lqr(Ap',Cp',1e0*eye(n_x), 1.5e1*eye(n_e));
L=L';

% ***********************
%
% RETURN STATE SPACE OF Ko
%
Ko = ss(Ap-Bp*F-L*Cp+L*Dp*F, -L, -F, 0);