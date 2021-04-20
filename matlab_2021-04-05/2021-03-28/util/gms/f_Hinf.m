function [value, sg, varargout] = f_Hinf ...
            (M, x, T11, T12, T21, Q, vec, varargin)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% COMPUTE H-INFINITY NORM AND SUBGRADIENTS AT OF PARAMETERIZED TRANSFER
% FUNCTION MATRICES GIVEN YOULA ET AL. PARAMETER Q
%
% Karan Puttannaiah
% Brent Wallace  
%
% 12/14/2020
%
% This program provides a subgradient for the function 
%
%               phi: RH-Infty -> R,     phi(M) = ||M||_{H-Infty}
%
% at the transfer function matrix M, where M is given by:
%
%               M = sum_{k=1}^{n} M_{k} * x_{k}
%
% Where M_{k} in RH-Infty, x_{k} in R, k = 1,...,n. This program provides
% the valid subgradient
%
%               g_{phi}: RH-Infty -> R^n,
%
%                                _                     _
%               g_{phi}(M) =    |  phi^{sg}(M_{1})      |
%                               |  phi^{sg}(M_{2})      |
%                               |       ...             |
%                               |  phi^{sg}(M_{n})      |
%                               |_                     _|
%
% Where phi^{sg}(M_{k}) = Re{ u_0^H * M_{k}(j*w_0) * v_0 }, where u_0 and
% v_0 are the left and right singular values of M(j*w_0) corresponding to
% the maximum singular vale sigma_max at the peak frequency w_0 at which
% the H-Infty norm of M is achieved. See Karan Puttannaiah's Ph.D. thesis,
% pp. 76 for further documentation.
%
% ***** CALL SYNTAX:
%
% [value, sg, varargout] = f_Hinf(M, x, T11, T12, T21, Q, vec, varargin)
%
% ***** INPUTS:
%
%   M                   : Vectorized Q parameterization transfer
%                           function matrices. See f_Vectorize.m for
%                           further documentation.
%   x                   : N*n_u*n_e dimensional query point vector.
%   T11                 : (1,1) transfer function matrix of Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   T12                 : (1,2) transfer function matrix of Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   T21                 : (2,1) transfer function matrix of Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   Q                   : Approximate Youla et al. Q parameter (state space
%                           object).
%   vec                 : A vector of indices which point to the relevant
%                           output channels (i.e., rows) of the Q
%                           parameterization transfer function matrix T(Q).
%   varargin            : If this program is being evaluated for the
%                           weighted mixed sensitivity portion of the
%                           optimization, this input argument is left
%                           blank. If this program is being evaluated for a
%                           constraint functional in the optimizer, then
%                           this input argument should be populated with
%                           the constraint functional's constraint value.
%
% ***** OUTPUTS:
%
%   xk                  : Suboptimal solution (within tolerance). N*n_u*n_e
%                           dimensional vector.
%   sg                  : Subgradient at matrix M (N*n_u*n_e dimensional
%                           vector).
%   varargout           : If this program is being evaluated for the
%                           weighted mixed sensitivity portion of the
%                           optimization, this output argument is left
%                           blank. If this program is being evaluated
%                           for a constraint functional in the optimizer,
%                           then this argument returns the contraint
%                           functional's constriant value.
%
% *************************************************************************
% *************************************************************************
% *************************************************************************


% ***********************
%
% INITIALIZATION
%

n = length(x);      % Dimension of problem.

% ***********************
%
% CHECK IF CONSTRAINT FUNCTIONAL IS BEING EVALUATED
%
% If so, set the output varargout to the constraint value of the constraint
% functional.
%
if nargin == 8
 conval = varargin{1};
 varargout{1} = conval;
end


% ***********************
%
% FORM Twz(Q) = T11 + T12 * Q * T21
%
if isempty(T11)
    Twz = ss([]); 
else
    Twz = parallel(T11,series(series(T21,Q),T12));
end

% Twz = minreal(Twz);

% Select output channels of interest.
Twz = Twz(vec,:);

% ***********************
%
% EVALUATE H-INFTY NORM OF Twz
%

% [ninf, fpeak] = norm(Twz, inf, 1e-8);

% Evaluate H-Infty norm of Twz.
[ninf, fpeak] = hinfnorm(Twz, 1e-8);
value = ninf;

% ***********************
%
% PRESCALE Twz TO GET ACCURATE SVD AT PEAK FREQUENCY
%

% If fpeak is near-zero, make it a small number.
if fpeak < 1e-5
    fpeak = 1e-5; 
end

% Get a window [wmin, wmax] which contains fpeak to prescale Twz at.
if fpeak<1e-5
    wmin=1e-8; wmax=1e-2;
elseif fpeak<1e-2
    wmin=1e-4; wmax=1e0;
elseif fpeak<1e0
    wmin=1e-1; wmax=1e1;
elseif fpeak<10
    wmin=1e-1; wmax=1e2;
elseif fpeak<1e2
    wmin=1e0; wmax=1e4;
elseif fpeak<1e5
    wmin=1e3; wmax=1e7;    
elseif fpeak>=1e5
    wmin=1e4; wmax=1e10;
else
    wmin=max([1, fpeak-10]); wmax=fpeak+10;
end

% Prescale Twz at frequency window [wmin, wmax].
% See       https://www.mathworks.com/help/control/ref/prescale.html   
TwzScaled = prescale(Twz,{wmin,wmax});

% Evaluate (prescaled) w -> z TFM at frequency fpeak for subsequent SVD.
Hjwo = freqresp(TwzScaled,fpeak);
% Hjwo = evalfr(TwzScaled,fpeak);

% ***********************
%
% SVD OF PRESCALED Twz AT PEAK FREQUENCY, MAX RIGHT/LEFT SINGULAR VECTORS
%

[U,S,V] = svd(Hjwo);        % SVD at peak frequency.
if isempty(U)
    uo = []; 
    vo = []; 
else
    uo 		= U(:,1); 		% Maximum Left Singular Vector.
    vo 		= V(:,1); 		% Maximum Right Singular Vector .
end

% ***********************
%
% EVALUATE SUBGRADIENT
%
subgradient = zeros(n,1);
for i = 1:n
    Hjwo = freqresp(M{i},fpeak);
%     magHjwo = abs(Hjwo);
    subgradient(i) = real(uo'*Hjwo*vo);
end
sg = subgradient;
