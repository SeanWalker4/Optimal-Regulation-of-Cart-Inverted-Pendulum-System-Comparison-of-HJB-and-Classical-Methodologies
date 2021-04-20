function K=f_FormK_ZamesParam(P_ss,Q,F,L);

% Form the controller from Q
% Uses Zames parameterization
% Author: Karan Puttannaiah
% Date: 09/19/2014

K=Q*inv(eye(size(P_ss*Q))-P_ss*Q);