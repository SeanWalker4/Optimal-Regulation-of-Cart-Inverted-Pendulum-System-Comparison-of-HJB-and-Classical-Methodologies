function q = f_Basis(N, p, z, basis_type)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% FORM N-DIMENSIAL BASIS FOR Q PARAMETER
%
% Karan Puttannaiah
% Brent Wallace  
%
% 12/14/2020
%
% This program forms an N-dimensional (approximation) basis for the
% infinite dimensional space RH-Infty.
%
%
% ***** CALL SYNTAX:
%
% q = f_Basis(N, p, z, basis_type)
%
% ***** INPUTS:
%
%   N                   : Basis dimension (integer).
%   p                   : Basis pole location (double).
%   z                   : Basis zero location (double).
%   basis_type          : Type of basis used (integer). Has the following
%                           options:
%                           1 = fixed pole low pass
%                           2 = fixed pole all pass
%                           3 = variable pole low pass
%                           4 = variable pole all pass
%                           5 = pole-zero
%                           6 = Laguerre
%
% ***** OUTPUTS:
%
%   q                   : N-dimensional cell array of zero-pole-gain (zpk)
%                           objects containing basis terms.
%
% *************************************************************************
% *************************************************************************
% *************************************************************************


q{1} = zpk([],[],1);


switch basis_type

    % ***********************
    %
    % FIXED POLE LOW PASS
    %
    case 1
    for k=2:N
        q{k} = zpk([],-p,p)^(k-1);
    end
    
    % ***********************
    %
    % FIXED POLE ALL PASS
    %
    case 2
    for k=2:N
        q{k} = zpk(p,-p,-1)^(k-1);
    end

    % ***********************
    %
    % VARIABLE POLE LOW PASS
    %
    case 3
    for k=2:N
        q{k} = zpk([],-p*(k-1),p*(k-1));
    end

    % ***********************
    %
    % VARIABLE POLE ALL PASS
    %
    case 4
    for k=2:N
        q{k} = zpk(p*(k-1),-p*(k-1),-1);
    end

    % ***********************
    %
    % POLE-ZERO
    %
    case 5
    for k=2:N 
        q{k} = zpk(z,-p,-1)^(k-1);
    end

    % ***********************
    %
    % LAUGUERRE
    %
    case 6
    for k=2:N
        q{k} = zpk([],-p,sqrt(2*p))*zpk(p,-p,1)^(k-1); 
    end

end
