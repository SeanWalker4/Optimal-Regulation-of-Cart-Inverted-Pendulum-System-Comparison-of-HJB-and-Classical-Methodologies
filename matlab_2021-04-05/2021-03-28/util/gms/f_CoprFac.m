function [T11rz, T12rz, T21rz, T11dz, T12dz, ...
            T21dz, T11niz, T12niz, T21niz] = f_CoprFac(...
            P_d, F, L, weights, settings)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% COPRIME FACTORIZATIONS AND Q PARAMETERIZATIONS
%
% Karan Puttannaiah
% Brent Wallace  
%
% 12/14/2020
%
% Forms the coprime factorizations and Q parameterizations of all weighting
% functions involved in the optimization.
%
% ***** CALL SYNTAX:
%
% [T11rz, T12rz, T21rz, T11dz, T12dz, ...
%             T21dz, T11niz, T12niz, T21niz] = f_CoprFac(...
%             P_d, F, L, weights, settings)
%
% ***** INPUTS:
%
%   P_d                 : Design plant (state space object).
%   F                   : MBC full state feedback gain matrix.
%   L                   : MBC observer gain matrix.
%   weights             : Data structure containing weightings (constructed
%                           in gms_main.m).
%   settings            : Data structure with the following fields:
%       loop_type       : Type of loop structure (string). See gms_main.m
%                           for further documentation.
%       youla1_zames0   : 1 = use Youla parameterization. 0 = use Zames
%                           parameterization.  
%       n_xi            : Number of states being fed back in inner loop
%                           (integer). Zero if inner-outer loop structure
%                           not used.
%       aug_integ       : 1 = augment plant at output with integrators to
%                           form design plant (if inner-outer, then only
%                           the output channels in the outer loop get
%                           augmented). 0 = don't augment.
%
% ***** OUTPUTS:
%
% Closed-loop map Q parameterizations are of the form:
%
%       T(Q) = T11 + T12 * Q * T21
%
%   T11rz               : r -> z (1,1) Q parameterization term.
%   T12rz               : r -> z (1,2) Q parameterization term.
%   T21rz               : r -> z (2,1) Q parameterization term.
%   T11dz               : d -> z (1,1) Q parameterization term.
%   T12dz               : d -> z (1,2) Q parameterization term.
%   T21dz               : d -> z (2,1) Q parameterization term.
%   T11niz              : ni -> z (1,1) Q parameterization term.
%   T12niz              : ni -> z (1,2) Q parameterization term.
%   T21niz              : ni -> z (2,1) Q parameterization term.
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

% *************************************************************************
% *************************************************************************
%
% INITIALIZATION
%
% *************************************************************************
% *************************************************************************

% ***********************
%
% LAPLACE TRANSFORM VARIABLE
%
s = tf('s');

% ***********************
%
% UNPACK SETTINGS
%

% Loop type.
loop_type = settings.loop_type;  

% Youla or Zames.
youla1_zames0 = settings.youla1_zames0;   

% Number of states fed back to inner loop.
n_xi = settings.n_xi;  

% If integral augmentation is to be performed.
aug_integ = settings.aug_integ;  


% ***********************
%
% PLANT STATE SPACE, DIMENSIONS
%
[Ap, Bp, Cp, Dp] = ssdata(P_d);

% Number of measured signals, number of control signals.
[n_e, n_u] = size(P_d);

% Number of measured signals in outer loop (relevant for inner-outer only).
n_yp = n_e - n_xi;          

% ***********************
%
% UNPACK WEIGHTING FUNCTIONS
%
W1 = weights.W1; 
W2 = weights.W2; 
W3 = weights.W3; 
Wd1 = weights.Wd1; 
Wd2 = weights.Wd2; 
Wd3 = weights.Wd3; 
Wni1 = weights.Wni1; 
Wni2 = weights.Wni2; 
Wni3 = weights.Wni3; 

% ***********************
%
% UNPACK CONSTRAINT FUNCTIONALS
%
W1c = weights.W1c; 
W2c = weights.W2c; 
W3c = weights.W3c; 
Wd1c = weights.Wd1c; 
Wd2c = weights.Wd2c; 
Wd3c = weights.Wd3c;
Wni1c = weights.Wni1c; 
Wni2c = weights.Wni2c; 
Wni3c = weights.Wni3c; 


% *************************************************************************
% *************************************************************************
%
% COPRIME FACTORIZATIONS OF PLANT AND CONTROLLER
%
% *************************************************************************
% *************************************************************************


% *************************************************************************
%
% RIGHT COPRIME FACTORIZATION -- PLANT
%

if youla1_zames0
    
    % ***********************
    %
    % YOULA PARAMETERIZATION
    %
    
    NumP.a = Ap - Bp*F; 
    NumP.b = Bp; 
    NumP.c = Cp - Dp*F; 
    NumP.d = Dp; 
    NumP = ss(NumP.a,NumP.b,NumP.c,NumP.d);
    
    DenP.a = Ap - Bp*F; 
    DenP.b = Bp; 
    DenP.c = -F; 
    DenP.d = eye(n_u); 
    DenP = ss(DenP.a,DenP.b,DenP.c,DenP.d);
    
else
    
    % ***********************
    %
    % ZAMES PARAMETERIZATION
    %
    
    NumP = P_d;
    DenP = ss(eye(n_u));
    
end


% *************************************************************************
%
% RIGHT COPRIME FACTORIZATION -- CONTROLLER
%

if youla1_zames0
    
    % ***********************
    %
    % YOULA PARAMETERIZATION
    %
    
    NumK.a = Ap - Bp*F; 
    NumK.b = -L; 
    NumK.c = -F; 
    NumK.d = zeros(n_u,n_e); 
    NumK = ss(NumK.a,NumK.b,NumK.c,NumK.d);
    
    DenK.a = Ap-Bp*F; 
    DenK.b = L; 
    DenK.c = Cp-Dp*F; 
    DenK.d = eye(n_e); 
    DenK = ss(DenK.a,DenK.b,DenK.c,DenK.d);
    
else
    
    % ***********************
    %
    % ZAMES PARAMETERIZATION
    %
    
    NumK = ss(zeros(n_u,n_e));
    DenK = ss(eye(n_e));
    
end


% *************************************************************************
%
% LEFT COPRIME FACTORIZATION -- PLANT
%

if youla1_zames0
    
    % ***********************
    %
    % YOULA PARAMETERIZATION
    %
    
    NumPt.a = Ap - L*Cp; 
    NumPt.b = Bp - L*Dp; 
    NumPt.c = Cp; 
    NumPt.d = Dp; 
    NumPt = ss(NumPt.a,NumPt.b,NumPt.c,NumPt.d);
    
    DenPt.a = Ap - L*Cp; 
    DenPt.b = -L; 
    DenPt.c = Cp; 
    DenPt.d = eye(n_e); 
    DenPt = ss(DenPt.a,DenPt.b,DenPt.c,DenPt.d);
    
else
    
    % ***********************
    %
    % ZAMES PARAMETERIZATION
    %
    
    NumPt = P_d;
    DenPt = ss(eye(n_e));
    
end


% *************************************************************************
%
% LEFT COPRIME FACTORIZATION -- CONTROLLER
%

if youla1_zames0
    
    % ***********************
    %
    % YOULA PARAMETERIZATION
    %
    
    NumKt.a = Ap - L*Cp; 
    NumKt.b = -L; 
    NumKt.c = -F;
    NumKt.d = zeros(n_u,n_e); 
    NumKt=ss(NumKt.a,NumKt.b,NumKt.c,NumKt.d);

    DenKt.a = Ap - L*Cp; 
    DenKt.b = -(Bp - L*Dp); 
    DenKt.c = -F; 
    DenKt.d = eye(n_u); 
    DenKt = ss(DenKt.a,DenKt.b,DenKt.c,DenKt.d);
    
else
    
    % ***********************
    %
    % ZAMES PARAMETERIZATION
    %
    
    NumKt = ss(zeros(n_u,n_e));
    DenKt = ss(eye(n_u));
    
end


% *************************************************************************
% *************************************************************************
%
% Q PARAMETERIZATION OF CLOSED LOOP MAPS -- BEFORE WIEGHTINGS AND
% INNER-OUTER I/O CHANNEL SELECTIONS
%
% *************************************************************************
% *************************************************************************

% S_e
SOut11 = DenK*DenPt; 
SOut12 = -NumP; 
SOut21 = DenPt;

% K*S_e
KSOut11 = NumK*DenPt; 
KSOut12 = DenP; 
KSOut21 = DenPt;

% T_e
TOut11 = NumP*NumKt; 
TOut12 = NumP; 
TOut21 = DenPt;

% S_c
SIn11 = DenP*DenKt; 
SIn12 = -DenP; 
SIn21 = NumPt;

% P*S_c
PSIn11 = NumP*DenKt; 
PSIn12 = -NumP; 
PSIn21 = NumPt;

% T_c
TIn11 = eye(n_u) - DenP*DenKt; 
TIn12 = DenP; 
TIn21 = NumPt;

% S_ni
Sni11 = SOut11; 
Sni12 = SOut12; 
Sni21 = SOut21;

% K*S_ni
KSni11 = KSOut11; 
KSni12 = KSOut12; 
KSni21 = KSOut21;

% T_ni
Tni11 = TOut11; 
Tni12 = TOut12; 
Tni21 = TOut21;


% *************************************************************************
% *************************************************************************
%
% I/O CHANNEL SELECTIONS (INNER-OUTER LOOP ONLY)
%
% *************************************************************************
% *************************************************************************

if strcmp(loop_type, 'in_out')
        
    % S_e
    SOut11 = SOut11(1:n_yp , 1:n_yp); 
    SOut12 = SOut12(1:n_yp , :); 
    SOut21 = SOut21(: , 1:n_yp);

    % K*S_e
    KSOut11 = KSOut11(: , 1:n_yp); 
    KSOut12 = KSOut12(: , :);           % No change 
    KSOut21 = KSOut21(:, 1:n_yp);

    % T_e
    TOut11 = TOut11(1:n_yp , 1:n_yp); 
    TOut12 = TOut12(1:n_yp , :); 
    TOut21 = TOut21(: , 1:n_yp);

    % S_c -- NO CHANGE

    % P*S_c
    PSIn11 = PSIn11(1:n_yp , :); 
    PSIn12 = PSIn12(1:n_yp , :); 
    PSIn21 = PSIn21(: , :);             % No change 

    % T_c -- NO CHANGE

    % S_ni
    Sni11 = SOut11(n_yp+1:end , n_yp+1:end); 
    Sni12 = SOut12(n_yp+1:end , :); 
    Sni21 = SOut21(: , n_yp+1:end);

    % K*S_ni
    KSni11 = KSni11(: , n_yp+1:end); 
    KSni12 = KSOut12(: , :);            % No change
    KSni21 = KSOut21(: , n_yp+1:end);

    % T_ni
    Tni11 = TOut11(1:n_yp , n_yp+1:end); 
    Tni12 = TOut12(1:n_yp , :); 
    Tni21 = TOut21(: , n_yp+1:end);

end


% *************************************************************************
% *************************************************************************
%
% INTRODUCE WEIGHTING FUNCTIONS TO Q PARAMETERIZATIONS
%
% *************************************************************************
% *************************************************************************

% Weighting function W1 (S_e).
T11rz1 = ss([]); 
T12rz1 = ss([]); 
T21rz1 = ss([]);
if weights.exist_W1
    
    T11rz1 = series(SOut11,W1); 
    T12rz1 = series(SOut12,W1);
    T21rz1 = SOut21;    
    
end

% Weighting function W2 (K*S_e).
T11rz2 = ss([]); 
T12rz2 = ss([]); 
T21rz2 = ss([]); 
if weights.exist_W2
    
    T11rz2 = series(KSOut11,W2); 
    T12rz2 = series(KSOut12,W2);
    T21rz2 = KSOut21;  
    
    if aug_integ
       
        T11rz2 = series(T11rz2,ss(1/(s+1e-1)));
        T12rz2 = series(T11rz2,ss(1/(s+1e-1)));
        
    end

end

% Weighting function W3 (T_e).
T11rz3 = ss([]); 
T12rz3 = ss([]); 
T21rz3 = ss([]); 
if weights.exist_W3
    
    T11rz3 = series(TOut11,W3); 
    T12rz3 = series(TOut12,W3);
    T21rz3 = TOut21;     
    
end

% Weighting function Wd1 (S_c).
T11dz1 = ss([]); 
T12dz1 = ss([]); 
T21dz1 = ss([]);   
if weights.exist_Wd1
    
    T11dz1 = series(SIn11,Wd1); 
    T12dz1 = series(SIn12,Wd1);
    T21dz1 = SIn21;   
    
end

% Weighting function Wd2 (P*S_c).
T11dz2 = ss([]); 
T12dz2 = ss([]); 
T21dz2 = ss([]);   
if weights.exist_Wd2
    
    T11dz2 = series(PSIn11,Wd2); 
    T12dz2 = series(PSIn12,Wd2);
    T21dz2 = PSIn21;     
    
    if aug_integ
       
        T11dz2 = series(T11dz2,ss(s+1e-6));
        T12dz2 = series(T12dz2,ss(s+1e-6));
        
    end
    
end

% Weighting function Wd3 (T_c).
T11dz3 = ss([]); 
T12dz3 = ss([]); 
T21dz3 = ss([]);   
if weights.exist_Wd3
    
    T11dz3 = series(TIn11,Wd3); 
    T12dz3 = series(TIn12,Wd3);
    T21dz3 = TIn21;    

end

% Weighting function Wni1 (S_ni).
T11niz1 = ss([]); 
T12niz1 = ss([]); 
T21niz1 = ss([]);   
if weights.exist_Wni1
    
    T11niz1 = series(Sni11,Wni1); 
    T12niz1 = series(Sni12,Wni1);
    T21niz1 = Sni21;    
    
end

% Weighting function Wni2 (K*S_ni).
T11niz2 = ss([]); 
T12niz2 = ss([]); 
T21niz2 = ss([]);   
if weights.exist_Wni2
    
    T11niz2 = series(KSni11,Wni2); 
    T12niz2 = series(KSni12,Wni2);
    T21niz2 = KSni21;      
    
end

% Weighting function Wni3 (T_ni).
T11niz3 = ss([]); 
T12niz3 = ss([]); 
T21niz3 = ss([]);   
if weights.exist_Wni3
    
    T11niz3 = series(Tni11,Wni3); 
    T12niz3 = series(Tni12,Wni3);
    T21niz3 = Tni21;     
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Added for Integ Augment fix
% s=tf('s');
% % T11rz1=series(W1*SOut11,ss(s+1e-6)); T12rz1=series(W1*SOut12,ss(s+1e-6)); T21rz1=SOut21;
% T11rz2=series(W2*KSOut11,ss(1/(s+1e-1))); T12rz2=series(W2*KSOut12,ss(1/(s+1e-1))); T21rz2=KSOut21;
% % % The eps in the integrator (1/(s+eps)) is picked to relatively high value
% % % ~0.1 because of a bad effect caused by bilinear transformation. When
% % % bilinear transformation is done, the weighted KS or
% % % (T11rz2+T12rz2**T21rz2) was going to a high value at low frequencies
% % % (near DC), even though the acutal weighted KS (i.e., without bilin) was
% % % low at those frequencies. 
% % T11dz1=series(Wd1*SensIn11,ss(s+1e-6)); T12dz1=series(Wd1*SensIn12,ss(s+1e-6)); T21dz1=SensIn21;
% T11dz2=series(Wd2*PSIn11,ss(s+1e-6)); T12dz2=series(Wd2*PSIn12,ss(s+1e-6)); T21dz2=PSIn21;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% *************************************************************************
% *************************************************************************
%
% CONTRAINT TRANSFER FUNCTION Q PARAMETERIZATIONS
%
% *************************************************************************
% *************************************************************************

% Constraint functional W1c (S_e).
T11rz1c = ss([]); 
T12rz1c = ss([]); 
T21rz1c = ss([]); 
for i=1:length(W1c)
    
    T11rz1c = [T11rz1c; series(SOut11,W1c{i}.tfm)]; 
    T12rz1c = [T12rz1c; series(SOut12,W1c{i}.tfm)];
    T21rz1c = [T21rz1c; SOut21];
    
end

% Constraint functional W2c (K*S_e).
T11rz2c = ss([]); 
T12rz2c = ss([]); 
T21rz2c = ss([]); 
for i=1:length(W2c)
    
    T11rz2c = [T11rz2c; series(KSOut11,W2c{i}.tfm)]; 
    T12rz2c = [T12rz2c; series(KSOut12,W2c{i}.tfm)];
    T21rz2c = [T21rz2c; KSOut21];
    
end

% Constraint functional W3c (T_e).
T11rz3c = ss([]); 
T12rz3c = ss([]); 
T21rz3c = ss([]); 
for i=1:length(W3c)
    
    T11rz3c = [T11rz3c; series(TOut11,W3c{i}.tfm)]; 
    T12rz3c = [T12rz3c; series(TOut12,W3c{i}.tfm)];
    T21rz3c = [T21rz3c; TOut21];
    
end

% Constraint functional Wd1c (S_c).
T11dz1c = ss([]); 
T12dz1c = ss([]); 
T21dz1c = ss([]); 
for i=1:length(Wd1c)
    
    T11dz1c = [T11dz1c; series(SIn11,Wd1c{i}.tfm)]; 
    T12dz1c = [T12dz1c; series(SIn12,Wd1c{i}.tfm)];
    T21dz1c = [T21dz1c; SIn21];
    
end

% Constraint functional Wd2c (P*S_c).
T11dz2c = ss([]); 
T12dz2c = ss([]); 
T21dz2c = ss([]); 
for i=1:length(Wd2c)
    
    T11dz2c = [T11dz2c; series(PSIn11,Wd2c{i}.tfm)]; 
    T12dz2c = [T12dz2c; series(PSIn12,Wd2c{i}.tfm)];
    T21dz2c = [T21dz2c; PSIn21];
    
end

% Constraint functional Wd3c (T_c).
T11dz3c = ss([]); 
T12dz3c = ss([]); 
T21dz3c = ss([]); 
for i=1:length(Wd3c)
    
    T11dz3c = [T11dz3c; series(TIn11,Wd3c{i}.tfm)]; 
    T12dz3c = [T12dz3c; series(TIn12,Wd3c{i}.tfm)];
    T21dz3c = [T21dz3c; TIn21];
    
end

% Constraint functional Wni1c (S_ni).
T11niz1c = ss([]); 
T12niz1c = ss([]); 
T21niz1c = ss([]); 
for i=1:length(Wni1c)
    
    T11niz1c = [T11niz1c; series(Sni11,Wni1c{i}.tfm)]; 
    T12niz1c = [T12niz1c; series(Sni12,Wni1c{i}.tfm)];
    T21niz1c = [T21niz1c; Sni21];
    
end

% Constraint functional Wni2c (K*S_ni).
T11niz2c = ss([]); 
T12niz2c = ss([]); 
T21niz2c = ss([]); 
for i=1:length(Wni2c)
    
    T11niz2c = [T11niz2c; series(KSni11,Wni2c{i}.tfm)]; 
    T12niz2c = [T12niz2c; series(KSni12,Wni2c{i}.tfm)];
    T21niz2c = [T21niz2c; KSni21];
    
end

% Constraint functional Wni3c (T_ni).
T11niz3c = ss([]); 
T12niz3c = ss([]); 
T21niz3c = ss([]); 
for i=1:length(Wni3c)
    
    T11niz3c = [T11niz3c; series(Tni11,Wni3c{i}.tfm)]; 
    T12niz3c = [T12niz3c; series(Tni12,Wni3c{i}.tfm)];
    T21niz3c = [T21niz3c; Tni21];
    
end


% *************************************************************************
% *************************************************************************
%
% FINALIZE OUTPUT
%
% *************************************************************************
% *************************************************************************

% Trz
T11rz = [T11rz1; T11rz2; T11rz3; T11rz1c; T11rz2c; T11rz3c]; 
T12rz = [T12rz1; T12rz2; T12rz3; T12rz1c; T12rz2c; T12rz3c]; 
T21rz = T21rz1;

% Tdz
T11dz = [T11dz1; T11dz2; T11dz3; T11dz1c; T11dz2c; T11dz3c]; 
T12dz = [T12dz1; T12dz2; T12dz3; T12dz1c; T12dz2c; T12dz3c]; 
T21dz = T21dz1;

% Tniz
T11niz = [T11niz1; T11niz2; T11niz3; T11niz1c; T11niz2c; T11niz3c]; 
T12niz = [T12niz1; T12niz2; T12niz3; T12niz1c; T12niz2c; T12niz3c]; 
T21niz = T21niz1;