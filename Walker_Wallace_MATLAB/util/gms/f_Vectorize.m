 function [M Mobj Mcon] = f_Vectorize ...
        (T11, T12, T21, q, n_u, n_e, ProblemData)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% VECTORIZE OPTIMIZATION PROBLEM
%
% Karan Puttannaiah
% Brent Wallace  
%
% 12/14/2020
%
% This program takes the matrix-valued Q parameterization
%
%           T_wz = T_11 + T_12 * [ sum_{k=1}^{N} X_k * q_k ] * T_21
%
% And transforms it into an equivalent (vectorized) Q parameterization
%
%           T_wz = M_o + sum_{l=1}^{nu*ne*N} M_l x_l
%
% where,
%
%   i = 1,...,n_u
%   j = 1,...,n_e
%   l = (k-1)*nu*ne+(j-1)*nu+i;
%   M_{l} = M_{k}^{ij}
%   M_{k}^{ij} = T_{12}*B^{ij}*T_{21}*q_{k}
%
%   (see Karan Puttannaiah's Ph.D. thesis for further documentation).
%
%
% ***** CALL SYNTAX:
%
% [M Mobj Mcon] = f_Vectorize(T11, T12, T21, qk, n_u, n_e, ProblemData)
%
% ***** INPUTS:
%
%   T11                 : Q Parameterization (1,1) element (state space
%                           object).
%   T12                 : Q Parameterization (1,2) element. (state space
%                           object).
%   T21                 : Q Parameterization (2,1) element. (state space
%                           object).
%   q                   : N-dimensional cell array of zero-pole-gain (zpk)
%                           objects containing basis terms.
%   n_u                 : Number of control signals (integer).
%   n_e                 : Number of measured signals (integer).
%   ProblemData         : Data structure containing problem data (see
%                           f_GenData.m for documentation).
%
% ***** OUTPUTS:
%
%   M                   : 1 by N*n_e*n_u  cell array containing vectorized
%                           transfer function matrices M_{l}. The first nz
%                           rows of each M_{l} corresponds to the number of
%                           regulated signals associated this loop breaking
%                           point. The remaining rows correspond to the
%                           ConNum-number of constraint functionals.
%   Mobj                : 1 by N*n_e*n_u cell array containing vectorized
%                           transfer function matrices corresponding to the
%                           mixed-sensitivity portion of the optimization.
%                           i.e., the l-th entry contains the first nz
%                           rows of the transfer function matrix M_{l}.
%   Mcon                : ConNum by N*n_e*n_u cell array containing the
%                           transfer function matrices corresponding to the
%                           constraint functional portion of the
%                           optimization, where ConNum is the number of
%                           constraint functionals associated with this
%                           loop breaking point. The i-th row corresponds
%                           to the i-th constraint functional, the l-th
%                           column to the l-th vectorized term.
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

% ***********************
%
% INITIALIZATION
%

% Basis dimension.
N = size(q,2);

% Number of states in T21.
size_t21 = size(T21.a,1);

% Number of states in T12.
size_t12 = size(T12.a,1);

% ***********************
%
% VECTORIZATION
%
Mobj = {};
Mcon = {};
Bij = zeros(n_u,n_e);
for k = 1:N
    for j = 1:n_e
        for i = 1:n_u
            l = (k-1)*n_u*n_e+(j-1)*n_u+i;
            Bij = zeros(n_u,n_e);
            Bij(i,j) = 1;
            if isempty(T12)
                M{l} = ss([]); 
            else
                % State space of M_{l} = T_{12} * B^{ij} * T_{21} * q_{k}
                a = [   T21.a               zeros(size_t21,size_t12);
                        T12.b*Bij*T21.c     T12.a                       ];
                b = [   T21.b
                        T12.b*Bij*T21.d     ];
                c = [   T12.d*Bij*T21.c     T12.c   ];
                d =     T12.d*Bij*T21.d;
                M{l} = ss(a,b,c,d) * q{k};
            end
        end
    end
end

% ***********************
%
% SELECT ROWS OF M_{l} CORRESPONDING TO WEIGHTED MIXED SENSITIVITY SIGNALS
%
for k = 1:N*n_e*n_u
    Mobj{k} = M{k}(ProblemData.ObjVec,:);
end

% ***********************
%
% SELECT ROWS OF M_{l} CORRESPONDING TO EACH RESPECTIVE CONSTRAINT
% FUNCTIONAL
%
for i = 1:ProblemData.ConNum
    for k = 1:N*n_e*n_u
        Mcon{i,k} = M{k}(ProblemData.ConVec{i},:);
    end
end
