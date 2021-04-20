function QN = f_FormQN(x, qk, n_u, n_e, N)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% FORM APPROXIMATE Q PARAMETER FROM FINITE DIMENSIONAL BASIS AND COORDINATE
% VECTOR x
%
% Karan Puttannaiah
% Brent Wallace  
%
% 12/14/2020
%
% Forms an approximate Q parameter Q_N from an N-dimensional basis and
% point x in R^{N*n_u*n_e}, where n_u is the number of control signals and
% n_e is the number of measured signals. Q_N is formed as follows:
%
%               Q_N = sum_{k=1}^{N} X_k * q_k
%
% Where q_k is the k-th term of the N-dimensional approximation basis, and
% X_k in R^{n_u by n_e} is the n_u by n_e coordinate of the k-th basis
% term. X_k is formed from x in R^{N*n_u*n_e} as follows:
%
%                        _                                             _
%                       | x_{k}^{11}  x_{k}^{12}     ... x_{k}^{1n_e}   |
%                       | x_{k}^{21}  x_{k}^{22}     ... x_{k}^{2n_e}   |
%               X_k =   |    ...        ...          ...    ...         |
%                       | x_{k}^{n_u1}  x_{k}^{n_u2} ... x_{k}^{n_un_e} | 
%                       |_                                             _|
%
% where,
%       x = [   x_{1}^{11} 
%               x_{1}^{12} 
%               ... 
%               x_{1}^{n_un_e} 
%               x_{2}^{11}
%               ...
%               x_{N}^{n_un_e}  ]
%
% See Karan Puttannaiah's Ph.D. thesis, pp 73 for further documentation.
%
%
% ***** CALL SYNTAX:
%
% QN = f_FormQN(x, qk, n_u, n_e, N)
%
% ***** INPUTS:
%
%   x                   : N*n_u*n_e dimensional coordinate vector.
%   qk                  : 1 by N cell array of zpk objects holding basis
%                           transfer function terms. 
%   n_u                 : Number of control signals (integer).
%   n_e                 : Number of regulated signals (integer).
%   N                   : Basis dimension (integer).
%
% ***** OUTPUTS:
%
%   Q                   : Approximate Q parameter (n_u by n_e state space
%                           object).
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

% ***********************
%
% INITIALIZATION
%

% Reshape coordinate vector x from N*n_u*n_e by 1 to n_u*n_e by N.
xtemp = reshape(x,n_u*n_e,N);

% Initialize approximate Q parameter.
QN = zeros(n_u, n_e);

% ***********************
%
% FORM APPROXIMATE Q PARAMETER
%

for i = 1:N

    % Form X_i in R^{n_u by n_e} by selecting the i-th column of xtemp (a
    % n_u*n_e dimensional vector) and reshaping it to a n_u by n_e matrix.
    X{i} = reshape(xtemp(:,i),n_u,n_e);
    
    % % Find temp = QN + X{i} * qk{i}
    % % Straight forward way.
    % temp = QN + X{i} * qk{i};

    % Alternative way. minreal in later works better when this is used,
    % i.e., order of QN found will be as expected.
    QNss = ss(QN);
    xqss = ss(series(qk{i},X{i}));

    % Form state space of additive system QNss + xqss.
    A = blkdiag(QNss.a,xqss.a);
    B = [   QNss.b
            xqss.b  ];
    C = [   QNss.c     xqss.c   ];
    D = QNss.d + xqss.d;
    temp = zpk(ss(A,B,C,D));
    
    % Take minimum realization to reduce system order.
    QN = minreal(temp,1e-6); 
    
end
QN = ss(QN);
