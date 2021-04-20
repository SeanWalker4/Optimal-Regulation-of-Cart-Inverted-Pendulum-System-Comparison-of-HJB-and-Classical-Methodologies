function [K, gamma_opt] = gms_main(sys_struct, opt_struct)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% GENERALIZED MIXED SENSITIVITY MAIN
%
% Karan Puttannaiah
% Brent Wallace  
%
% 12/14/2020
%
% Main File for Control Design using Generalized Mixed Sensitivity(GMS)
% Computes a H-infinity based Feedback Controller based on
% multiobjective constrained convex optimization.
%
% ***** CALL SYNTAX:
%
% [K, gamma_opt] = gms_main(sys_struct, opt_struct)
%
% ***** INPUTS:
%
%   sys_struct          : Data structure with the following fields:
%       settings        : Data structure with the following fields: 
%           loop_type   : Type of loop structure (string) with the
%                           following options:
%                           'std_pk'    : Standard P-K feedback structure
%                           'in_out'    : Inner-outer feedback structure
%           aug_integ   : 1 = augment plant at output with integrators to
%                           form design plant (if inner-outer, then only
%                           the output channels in the outer loop get
%                           augmented). 0 = don't augment.
%           bilin_param : Bilinear transform (BLT) parameters. NOTE:
%                           declare only if BLT is to be performed. Data
%                           structure with the following fields:
%               p1      : Bilinear transform parameter (see above).
%               p2      : Bilinear transform parameter (see above).
%
%                           Transform has form
%
%                               s - p1
%                       s' =    ------
%                               1 - s/p2
%
%           youla1_zames0: 1 = use Youla parameterization. 0 = use Zames
%                           parameterization.   
%       P               : Data structure with the following fields:
%           P_ss        : State space representation of plant P (state
%                           space object)
%           Mi          : (INNER-OUTER ONLY) A matrix with following
%                           properties:
%                           num of col = num of states
%                           num of row = num of states being fed back in
%                           inner loop
%                           Eg: feeding back states 3 and 4:
%                           Mi =   [0 0 1 0 0 0;
%                                   0 0 0 1 0 0];
%       weights         : Closed-loop map weighting functions (state-space
%                           objects) and constraint functionals
%                           (strunctures). NOTE: declare only if the
%                           respective weighting function/constraint
%                           functional is involved in the optimization.
%        WEIGHTING FUNCTIONS:
%                           Each is a state space object.
%           W1          : Weighting function W1 (S_e).
%           W2          : Weighting function W2 (K*S_e).
%           W3          : Weighting function W3 (T_e).
%           Wd1         : Weighting function W4 (S_c).
%           Wd2         : Weighting function W5 (P*S_c).
%           Wd3         : Weighting function W6 (T_c).
%           Wni1        : Weighting function W7 (S_ni).
%           Wni2        : Weighting function W8 (K*S_ni).
%           Wni3        : Weighting function W9 (T_ni).
%        CONSTRAINT FUNCTIONALS: 
%                           Each is a struct with the following fields:
%                               tfm     : Transfer function matrix of
%                               constraint functional (state space object).
%                               Fun     : Functional type (string). Has
%                                           following options:
%                                   'f_Hinf' = H-Infty constraint.
%                                   'f_Linf' = L-Infty constraint.
%                               Val     : Functional constraint in
%                                           optimization (double).
%           W1c         : Constraint functional W1c (S_e).
%           W2c         : Constraint functional W2c (K*S_e).
%           W3c         : Constraint functional W3c (T_e).
%           Wd1c        : Constraint functional W4c (S_c).
%           Wd2c        : Constraint functional W5c (P*S_c).
%           Wd3c        : Constraint functional W6c (T_c).
%           Wni1c       : Constraint functional W7c (S_ni).
%           Wni2c       : Constraint functional W8c (K*S_ni).
%           Wni3c       : Constraint functional W9c (T_ni).    
%     
%   opt_struct          : Data structure with the following fields:
%       settings        : Optimization settings. Data structure with the
%                           following fields:
%           xmin1       : Lower bound on optimization solution (double).
%                           i.e., if nu is the number of controls, ny the
%                           number of measured signals, N the dimension of
%                           the basis, then nu*ne*N optimum xopt lies in
%                           the box [xmin1*ones(nu*ne*N),
%                           xmax1*ones(nu*ne*N)]
%           xmax1       : Upper bound on optimization solution (double).
%                           i.e., if nu is the number of controls, ny the
%                           number of measured signals, N the dimension of
%                           the basis, then nu*ne*N optimum xopt lies in
%                           the box [xmin1*ones(nu*ne*N),
%                           xmax1*ones(nu*ne*N)]
%           x01         : Initial condition for the optimization (double).
%                           i.e., if nu is the number of controls, ny the
%                           number of measured signals, N the dimension of
%                           the basis, then the optimization will
%                           initialize at x01*ones(nu*ne*N) . NOTE: must
%                           lie in the box [xmin1*ones(nu*ne*N),
%                           xmax1*ones(nu*ne*N)]
%           maxiter     : Maximum number of iterations before algorithm
%                           terminates (integer).
%           sum1_max0   : 1 = perform weighted sum optimization. 0 =
%                           perform weighted max optimization.
%           algo        : Algorithm to use (integer). Has the following
%                           options:
%                           1 = ACCPM
%                           2 = Kelly's CPM
%                           3 = SolvOpt
%       basis           : Basis parameters. Data structure with the
%                           following fields:
%           type        : Type of basis used (integer). Has the following
%                           options:
%                           1 = fixed pole low pass
%                           2 = fixed pole all pass
%                           3 = variable pole low pass
%                           4 = variable pole all pass
%                           5 = pole-zero
%                           6 = Laguerre
%           N           : Number of basis terms (integer).
%           p           : Pole location in basis (double). 
%           z           : Zero location in basis (double). 
%       obj_funct_w     : Objective function weighting parameters
%                           (nonnegative double). NOTE: declare only if the
%                           respective closed-loop map is involved in the
%                           optimization.
%           mu1         : Weighting on W1 (S_e).
%           mu2         : Weighting on W2 (K*S_e).
%           mu3         : Weighting on W3 (T_e).
%           rho1        : Weighting on W4 (S_c).
%           rho2        : Weighting on W5 (P*S_c).
%           rho3        : Weighting on W6 (T_c).
%           eta1        : Weighting on W7 (S_ni).
%           eta2        : Weighting on W8 (K*S_ni).
%           eta3        : Weighting on W9 (T_ni).
%           
%
%
% ***** OUTPUTS:
%
%   K                   : GMS H-Inf optimal controller (state space
%                           object).
%   gamma_opt           : Optimal performance measure (double).
%
% 
% EXECUTION PROCEDURE:
%
%   - Form the design plant:
%       - Define the original plant
%       - Integrator augmentation if needed
%       - Bilinear transformation values if needed
%
%   - Select weighting functions:
%       - Tradeoff param rho
%       - W for obj
%       - W for constraint
%
%   - Select optimization params:
%       - LB and UB
%       - Init point
%       - Maximum number of iterations
%
%   - Select Youla/Zames parametrization:
%       - Select Youla or Zames
%       - Initial controller
%
%   - Finite Dimensionality
%       - Basis params
%
%   - Objective function:
%       - sum/max/stacking
%
%   - Find initial controller (Ko, F, L)
%
%   - Youla parameterization
%
%   - Find Initial Q parameter using initial controller (Ko, F, L)
%
%   - Extract required data from problem setup
%
%   - Vectorize the optimization problem
%
%   - Optimization process
%       - define how subgradient is picked based on sum/max/stacking
%
%   - form Q using the optimized variables and bases
%
%   - form Controller K using the obtained Q
%
%   - Inverse bilinear transformation if needed
%
%   - Inverse of integrator augmentation if needed
%
%   - Compute OL and CL maps
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% START PROGRAM
%
% *************************************************************************
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
% *************************************************************************



% Transfer function variable
s = tf('s');

% *************************************************************************
% *************************************************************************
%
% UNPACK INPUT DATA STRUCTURES
%
% *************************************************************************
% *************************************************************************

%%
% *************************************************************************
%
% UNPACK sys_struct
%
% *************************************************************************


% *************************************************************************
%
% UNPACK SETTINGS
%

% Loop type (standard P-K structure or inner-outer structure).
loop_type = sys_struct.settings.loop_type;

% Integral augmentation at plant output control.
aug_integ = sys_struct.settings.aug_integ;

% Check if bilinear transform is to be performed.
bilinear = isfield(sys_struct.settings, 'bilin_param');

% Bilinear transform parameters (if relevant).
if bilinear
    
    p1 =  sys_struct.settings.bilin_param.p1;   % p1
    p2 =  sys_struct.settings.bilin_param.p2;   % p2
    
end

% Youla parameterization or Zames parameterization (Youla = 1, Zames = 0).
youla1_zames0 = sys_struct.settings.youla1_zames0;

% *************************************************************************
%
% UNPACK PLANT PARAMETERS
%

% Plant (state space object).
P_ss = sys_struct.P.P_ss;

% Get plant state space data.
[Ap,Bp,Cp,Dp] = ssdata(P_ss); 

% Get plant dimensions.
ni_p = size(Bp,2);          % Number of plant inputs.
no_p = size(Cp,1);          % Number of plant outputs.
nx_p = size(Ap,1);          % Number of plant states.

% States fed back in inner loop (if inner-outer structure is used).
if strcmp(loop_type, 'in_out')      % Inner-outer
   
    % Make sure matrix Mi is defined. If not, throw an error.
    if isfield(sys_struct.P, 'Mi')
        
        Mi = sys_struct.P.Mi;
        
    else
        
        error(['***** ERROR: Inner-Outer loop structure specified, ' ...
            'but variable Mi not declared. Please specify states to ' ...
            'be fed back in inner loop. See gms_main.m for ' ...
            'documentation. *****'])
        
    end
    
    % Get number of states to be fed back to inner loop.
    n_xi = size(Mi,1);
    
    
else                                % Not inner-outer
    
    n_xi = 0;
    
end


% *************************************************************************
%
% UNPACK WEIGHTING FUNCTIONS, MULTIPLY BY SCALAR WEIGHTING
%
% Apply scalar weightings to closed-loop map weighting functions for
% optimization.
% mu corresponds to weight on properties at plant output,
% rho corresponds to weight on properties at plant input.
% If eta is defined, it corresponds to weight on properties at sensor noise
% for the inner-outer loop case.
%

% Weighting function W1 (S_e).
weights.exist_W1 = isfield(sys_struct.weights, 'W1'); 
if weights.exist_W1                 % Multiply by scalar weighting.
    w = opt_struct.obj_funct_w.mu1;
    weights.W1 = w * sys_struct.weights.W1;
else                                % Else, leave field empty.
    weights.W1 = [];
end


% Weighting function W2 (K*S_e).
weights.exist_W2 = isfield(sys_struct.weights, 'W2'); 
if weights.exist_W2                 % Multiply by scalar weighting.
    w = opt_struct.obj_funct_w.mu2;
    weights.W2 = w * sys_struct.weights.W2;
else                                % Else, leave field empty.
    weights.W2 = [];
end

% Weighting function W3 (T_e).
weights.exist_W3 = isfield(sys_struct.weights, 'W3');
if weights.exist_W3                 % Multiply by scalar weighting.
    w = opt_struct.obj_funct_w.mu3;
    weights.W3 = w * sys_struct.weights.W3;
else                                % Else, leave field empty.
    weights.W3 = [];
end

% Weighting function Wd1 (S_c).
weights.exist_Wd1 = isfield(sys_struct.weights, 'Wd1'); 
if weights.exist_Wd1                % Multiply by scalar weighting.
    w = opt_struct.obj_funct_w.rho1;
    weights.Wd1 = w * sys_struct.weights.Wd1;
else                                % Else, leave field empty.
    weights.Wd1 = [];
end

% Weighting function Wd2 (P*S_c).
weights.exist_Wd2 = isfield(sys_struct.weights, 'Wd2'); 
if weights.exist_Wd2                % Multiply by scalar weighting.
    w = opt_struct.obj_funct_w.rho2;
    weights.Wd2 = w * sys_struct.weights.Wd2;
else                                % Else, leave field empty.
    weights.Wd2 = [];
end

% Weighting function Wd3 (T_c).
weights.exist_Wd3 = isfield(sys_struct.weights, 'Wd3'); 
if weights.exist_Wd3                % Multiply by scalar weighting.
    w = opt_struct.obj_funct_w.rho3;
    weights.Wd3 = w * sys_struct.weights.Wd3;
else                                % Else, leave field empty.
    weights.Wd3 = [];
end

% Weighting function Wni1 (S_ni).
weights.exist_Wni1 = isfield(sys_struct.weights, 'Wni1');
if weights.exist_Wni1               % Multiply by scalar weighting.
    w = opt_struct.obj_funct_w.eta1;
    weights.Wni1 = w * sys_struct.weights.Wni1;
else                                % Else, leave field empty.
    weights.Wni1 = [];
end

% Weighting function Wni2 (K*S_ni).
weights.exist_Wni2 = isfield(sys_struct.weights, 'Wni2');
if weights.exist_Wni2               % Multiply by scalar weighting.
    w = opt_struct.obj_funct_w.eta2;
    weights.Wni2 = w * sys_struct.weights.Wni2;
else                                % Else, leave field empty.
    weights.Wni2 = [];
end

% Weighting function Wni3 (T_ni).
weights.exist_Wni3 = isfield(sys_struct.weights, 'Wni3');
if weights.exist_Wni3               % Multiply by scalar weighting.
    w = opt_struct.obj_funct_w.eta3;
    weights.Wni3 = w * sys_struct.weights.Wni3;
else                                % Else, leave field empty.
    weights.Wni3 = [];
end


% *************************************************************************
%
% UNPACK CONSTRAINT FUNCTIONALS
%


% Constraint functional W1c (S_e).
weights.exist_W1c = isfield(sys_struct.weights, 'W1c');
if weights.exist_W1c                % Multiply by scalar weighting.
    weights.W1c = w * sys_struct.weights.W1c;
else                                % Else, leave field empty.
    weights.W1c = [];
end

% Constraint functional W2c (K*S_e).
weights.exist_W2c = isfield(sys_struct.weights, 'W2c'); 
if weights.exist_W2c                % Multiply by scalar weighting.
    weights.W2c = w * sys_struct.weights.W2c;
else                                % Else, leave field empty.
    weights.W2c = [];
end

% Constraint functional W3c (T_e).
weights.exist_W3c = isfield(sys_struct.weights, 'W3c');
if weights.exist_W3c                % Multiply by scalar weighting.
    weights.W3c = w * sys_struct.weights.W3c;
else                                % Else, leave field empty.
    weights.W3c = [];
end

% Constraint functional Wd1c (S_c).
weights.exist_Wd1c = isfield(sys_struct.weights, 'Wd1c');
if weights.exist_Wd1c               % Multiply by scalar weighting.
    weights.Wd1c = w * sys_struct.weights.Wd1c;
else                                % Else, leave field empty.
    weights.Wd1c = [];
end

% Constraint functional Wd2c (P*S_c).
weights.exist_Wd2c = isfield(sys_struct.weights, 'Wd2c');
if weights.exist_Wd2c               % Multiply by scalar weighting.
    weights.Wd2c = w * sys_struct.weights.Wd2c;
else                                % Else, leave field empty.
    weights.Wd2c = [];
end

% Constraint functional Wd3c (T_c).
weights.exist_Wd3c = isfield(sys_struct.weights, 'Wd3c');
if weights.exist_Wd3c              % Multiply by scalar weighting.
    weights.Wd3c = w * sys_struct.weights.Wd3c;
else                                % Else, leave field empty.
    weights.Wd3c = [];
end

% Constraint functional Wni1c (S_ni).
weights.exist_Wni1c = isfield(sys_struct.weights, 'Wni1c');
if weights.exist_Wni1c              % Multiply by scalar weighting.
    weights.Wni1c = w * sys_struct.weights.Wni1c;
else                                % Else, leave field empty.
    weights.Wni1c = [];
end

% Constraint functional Wni2c (K*S_ni).
weights.exist_Wni2c = isfield(sys_struct.weights, 'Wni2c');
if weights.exist_Wni2c              % Multiply by scalar weighting.
    weights.Wni2c = w * sys_struct.weights.Wni2c;
else                                % Else, leave field empty.
    weights.Wni2c = [];
end

%  Constraint functional Wni3c (T_ni). 
weights.exist_Wni3c = isfield(sys_struct.weights, 'Wni3c');
if weights.exist_Wni3c              % Multiply by scalar weighting.
    weights.Wni3c = w * sys_struct.weights.Wni3c;
else                                % Else, leave field empty.
    weights.Wni3c = [];
end


%%
% *************************************************************************
%
% UNPACK opt_struct
%
% *************************************************************************


% *************************************************************************
%
% UNPACK SETTINGS
%


% Bounding box for solution, initial condition.
xmin1 = opt_struct.settings.xmin1;          % Lower bound.
xmax1 = opt_struct.settings.xmax1;          % Upper bound.
x01 = opt_struct.settings.x01;              % Initial condition.

% Optimizer settings.
maxiter = opt_struct.settings.maxiter;      % Maximum number of iterations.
sum1_max0 = opt_struct.settings.sum1_max0;  % Sum or max formulation.
algo = opt_struct.settings.algo;            % Algorithm type.



% *************************************************************************
%
% UNPACK BASIS PARAMETERS
%

basis_type = opt_struct.basis.type;         % Basis type.
N = opt_struct.basis.N;                     % Basis dimension.
p = opt_struct.basis.p;                     % Basis pole.
z = opt_struct.basis.z;                     % Basis zero.



%%
% *************************************************************************
% *************************************************************************
%
% INTEGRAL AUGMENTATION AT PLANT OUTPUT
%
% *************************************************************************
% *************************************************************************


if aug_integ
    
    % Form augmented plant state space.
    Ad =    [   Ap      zeros(nx_p, no_p)
                Cp      zeros(no_p)         ];
    
    Bd =    [   Bp
                Dp  ];
    
    Cd =    [   zeros(no_p, nx_p)   eye(no_p)	];
    
    Dd =    zeros(no_p, ni_p);
    
    % Form augmented plant.
    P_d = ss(Ad,Bd,Cd,Dd);
    
else
    
    % Else, copy original plant without augmentation.
    Ad = Ap;
    Bd = Bp;
    Cd = Cp;
    Dd = Dp;
    
    P_d = ss(Ad,Bd,Cd,Dd);
    
end

        

        
%%
% *************************************************************************
% *************************************************************************
%
% EXTEND PLANT OUTPUT TO INCLUDE STATES FED BACK TO INNER LOOP (INNER-OUTER
% ONLY)
%
% *************************************************************************
% *************************************************************************
        

if strcmp(loop_type, 'in_out')      % Check if inner-outer structure used.
    
    % Extend plant output to include states fed back to inner loop.
    % NOTE: Plant C matrix will have increased in number of columns if
    % integral augmentation at output was selected. This is handled
    % directly below.
    Cd =    [               Cd      
                [ Mi  zeros(n_xi, size(Cd,2)-nx_p) ]   ];
            
    Dd =    [   Dd
                zeros(n_xi,ni_p) ];
           
    % Form extended plant.
    P_d = ss(Ad,Bd,Cd,Dd);        
  
end

% Size of design plant.
[n_e, n_u] = size(P_d);


%%
% *************************************************************************
% *************************************************************************
%
% BILINEAR TRANSFORMATION
%
% *************************************************************************
% *************************************************************************

if bilinear
    
    % Backup plant before bilinear transform.
    P_d_BeforeBilin = P_d;     
    
    % Perform (forward) bilinear transform.
    [Adb,Bdb,Cdb,Ddb] = bilin(Ad,Bd,Cd,Dd,1,'Sft_jw',[p2 p1]);
    
    % Form transformed plant
    P_d = ss(Adb,Bdb,Cdb,Ddb);
    
    
end


%%
% *************************************************************************
% *************************************************************************
%
% FORM NOMINAL CONTROLLER
%
% *************************************************************************
% *************************************************************************

if youla1_zames0

    % ***********************
    %
    % YOULA PARAMETERIZATION
    %
    
    [Ko,F,L]=f_KNominal(P_d);

else
    
    % ***********************
    %
    % ZAMES PARAMETERIZATION
    %
    % Ko = 0. F = 0, L = 0.
    
    Ko = zeros(n_e,n_u) * ss(0,0,0,0);
    
    F = zeros(n_u,n_e);
    L = zeros(n_e,n_u);
    
end


%%
% *************************************************************************
% *************************************************************************
%
% COPRIME FACTORIZATION
%
% *************************************************************************
% *************************************************************************

% ***********************
%
% SETTTINGS FOR COPRIME FACTORIZATION FUNCTION
%

% Loop type.
settings_coprfac.loop_type = loop_type;  

% Youla or Zames.
settings_coprfac.youla1_zames0 = youla1_zames0;   

% Number of states fed back to inner loop.
settings_coprfac.n_xi = n_xi;  

% If integral augmentation is to be performed.
settings_coprfac.aug_integ = aug_integ;


% ***********************
%
% CALL FUNCTION
%

[T11rz, T12rz, T21rz, T11dz, T12dz, T21dz, T11niz, T12niz, T21niz] = ...
    f_CoprFac(P_d, F, L, weights, settings_coprfac);

%% 
% *************************************************************************
% *************************************************************************
%
% INITIAL Q PARAMETER
%
% *************************************************************************
% *************************************************************************


% ***********************
%
% FORM VECTOR OF LOWER BOUND, UPPER BOUND, AND INITIAL CONDITION
%
xmin = xmin1*ones(N*n_u*n_e,1);         % LB vector.
xmax = xmax1*ones(N*n_u*n_e,1);      	% UB vector.
x0 = x01*ones(N*n_u*n_e,1);             % IC vector.

% ***********************
%
% FORM BASIS
%
q = f_Basis(N, p, z, basis_type);

% ***********************
%
% FORM INITIAL Q PARAMETER
%
Q = f_FormQN(x0, q, n_u, n_e, N);


%% 
% *************************************************************************
% *************************************************************************
%
% GENERATE PROBLEM DATA FOR OPTIMIZATION
%
% *************************************************************************
% *************************************************************************

[Datarz, Datadz, Dataniz] = f_GenData(P_d, weights, n_xi);


%% 
% *************************************************************************
% *************************************************************************
%
% VECTORIZATION
%
% *************************************************************************
% *************************************************************************

% r -> z
[Mrz, Mobjrz, Mconrz] = f_Vectorize ...
    (T11rz, T12rz, T21rz, q, n_u, n_e, Datarz);

% d -> z
[Mdz, Mobjdz, Mcondz] = f_Vectorize ...
    (T11dz, T12dz, T21dz, q, n_u, n_e, Datadz);

% ni -> z
[Mniz, Mobjniz, Mconniz] = f_Vectorize ...
    (T11niz, T12niz, T21niz, q, n_u, n_e, Dataniz);



%% 
% *************************************************************************
% *************************************************************************
%
% RUN OPTIMIZATION
%
% *************************************************************************
% *************************************************************************

% ***********************
%
% SETUP
%
dim = length(x0);         % Dimension of problem


if algo == 1
    % -------- ACCPM -------- %
    switch sum1_max0
        case 0
            % % Weighted Minmax
            
            [xk,fx,iter_cnt,perf_meas]=...
                f_ACCPM_GenMixSens_Optimizer(dim,N,x0,Mobjrz,Mobjdz,...
                Mconrz,Mcondz,T11rz, T12rz, T21rz,T11dz, T12dz, ...
                T21dz,Datarz,Datadz,Q,q,n_u,n_e,xmax,xmin,maxiter);
            
%             %--------- Added for HSV IO WITH Tniu ---------%
%             [xk,fx,iter_cnt,perf_meas]=...
%                 f_ACCPM_GenMixSens_Optimizer_With_Tniu...
%                 (N,NQ,xk,Mobjrz,Mobjdz,Mobjniz,Mconrz,Mcondz,...
%                 Mconniz,T11rz, T12rz, T21rz,T11dz, T12dz, T21dz,...
%                 T11niz,T12niz,T21niz, Datarz, Datadz,Dataniz, Q,q,...
%                 n_u,n_e,xmax,xmin,MaxIter);
%             
        case 1
            % % Weighted Sum
            [xk,fx,iter_cnt,perf_meas]=...
                f_ACCPM_GenMixSens_Optimizer_Sum...
                (dim,N,x0,Mobjrz,Mobjdz,Mconrz,Mcondz,...
                T11rz, T12rz, T21rz,T11dz, T12dz, T21dz,...
                Datarz, Datadz, Q,q,n_u,n_e,xmax,xmin,maxiter);
    end
    
elseif algo == 2
    
    % -------- KELLEY'S CPM -------- %
    
    switch sum1_max0
        case 0
            %         Weighted Minmax
            
            [xk, frz, fdz, fnz] = f_KelleyCPM_GenMix_Optimizer ...
                (dim, N, x0, Mobjrz, Mobjdz, Mobjniz, ...
                Mconrz, Mcondz, Mconniz, ...
                T11rz, T12rz, T21rz, T11dz, T12dz, T21dz, ...
                T11niz, T12niz, T21niz, ...
                Datarz, Datadz, Dataniz, ...
                Q, q, n_u, n_e, maxiter, xmax, xmin);
            
%             %--------- Added for HSV IO WITH Tniu ---------%
%             [xk,frz,fdz]=f_KelleyCPM_GenMix_Optimizer_With_Tniu...
%                 (N,NQ,xk,Mobjrz,Mobjdz,Mobjniz,Mconrz,Mcondz,T11rz,...
%                 T12rz, T21rz,T11dz, T12dz, T21dz,T11niz,T12niz,...
%                 T21niz, Datarz, Datadz,Dataniz, Q,q,n_u,n_e,MaxIter,...
%                 xmax,xmin);
            
        case 1
            %         Weighted sum
            [xk,fo]=f_KelleyCPM_GenMix_Optimizer_Sum...
                (dim,N,x0,Mobjrz,Mobjdz,Mconrz,Mcondz,T11rz,...
                T12rz, T21rz,T11dz, T12dz, T21dz, Datarz,Datadz,...
                Q,q,n_u,n_e,maxiter,xmax,xmin);
    end
    
elseif algo == 3
    % -------- SOLVOPT -------- %
    
    opts(1) = -1;       % negative => minimization
    opts(2) = 1e-4;
    opts(3) = 1e-4;
    opts(4) = maxiter;  % default num iter 15000
    opts(5) = 0;        % 1->verbose, 0->silent
    N = N/(n_u*n_e);
    [xk_solvopt,fx_solvopt,opts_solvopt] = solvopt(x0,...
        @(x)solvopt_fval(x,N,Mobjrz,Mobjdz,Mconrz,Mcondz,T11rz,...
        T12rz, T21rz,T11dz, T12dz, T21dz, Datarz, Datadz, q,n_u,...
    n_e),@(x)solvopt_sg(x,N,Mobjrz,Mobjdz,Mconrz,Mcondz,T11rz,...
    T12rz, T21rz,T11dz, T12dz, T21dz, Datarz, Datadz, q,n_u,n_e),opts);
    
else
    disp('Choose a valid algorithm')
end



%% 
% *************************************************************************
% *************************************************************************
%
% FORM OPTIMAL Q PARAMETER
%
% *************************************************************************
% *************************************************************************


Q = f_FormQN(xk, q, n_u, n_e, N);


% %%
% [forz, Gfo] =feval('f_Hinf', Mobjrz, xk, T11rz, T12rz, T21rz, Q, Datarz.ObjVec);
% [fodz, Gfo] =feval('f_Hinf', Mobjdz, xk, T11dz, T12dz, T21dz, Q, Datadz.ObjVec);
% fx=max([forz,fodz]);


%% 
% *************************************************************************
% *************************************************************************
%
% FORM OPTIMAL CONTROLLER FROM OPTIMAL Q PARAMETER
%
% *************************************************************************
% *************************************************************************


if youla1_zames0
    
    K = f_FormK(P_d,Q,F,L);                % Youla parameterization
    
else
    
    K = f_FormK_ZamesParam(P_d,Q,F,L);     % Zames parameterization 
    
end

% T_copr.T11rz = T11rz;
% T_copr.T12rz = T12rz;
% T_copr.T21rz = T21rz;
% 
% T_copr.T11dz = T11dz;
% T_copr.T12dz = T12dz;
% T_copr.T21dz = T21dz;

% --------- Added for Hypersonic inner-outer WITH Tniu ---------%
% T_copr.T11niz = T11niz;
% T_copr.T12niz = T12niz;
% T_copr.T21niz = T21niz;


%% 
% *************************************************************************
% *************************************************************************
%
% INVERSE BILINEAR TRANSFORMATION
%
% *************************************************************************
% *************************************************************************


if bilinear
    
    % Get state space data from pre inverse bilin controller K.
    [Acp1,Bcp1,Ccp1,Dcp1] = ssdata(K);
    
    % Backup K before inverse bilin transformation.
    K_BeforeInvBilin = K; 
    
    % Perform inverse bilin transformation.
    [Atk1,Btk1,Ctk1,Dtk1] = bilin(Acp1,Bcp1,Ccp1,Dcp1,-1,'Sft_jw',[p2 p1]);
    
    % Post inverse bilin transform controller K.
    K = ss(Atk1,Btk1,Ctk1,Dtk1);
    
end


%% 
% *************************************************************************
% *************************************************************************
%
% INTEGRAL AUGMENTATION
%
% *************************************************************************
% *************************************************************************

% Take minimum realization of K.
% K = minreal(K);


if aug_integ
    
    % Augment K with integrators at the input in the first n_e - n_xi
    % channels. 
    
    % NOTE: If a standard P-K feedback structure is used, there is no inner
    % loop (i.e., no states fed back in the inner loop, so n_xi = 0). In
    % this case, n_xi has already been declared as 0, so K will be
    % augmented in all output channels, as desired.
    
    % Augmentation for first n_e - n_xi channels.
    aug = 1/s * eye(n_e - n_xi);    
    
    % Include last n_xi channels (no augmentation in these channels).
    aug = blkdiag(aug, eye(n_xi));
    
    % Augment K at its input.
    K = series(aug, K);
        
end



%% 
% *************************************************************************
% *************************************************************************
%
% OPTIMAL PERFORMANCE MEASURE gamma_opt
%
% Requires forming closed loop maps, taking infinity norm of stacked closed
% loop maps (with wieghting functions) at each loop breaking point, then
% taking the maximum.
%
% *************************************************************************
% *************************************************************************


% *************************************************************************
%
% FORM CLOSED LOOP MAPS
%

switch loop_type
    
    case 'std_pk'
        
        [Lo,Li,So,Si,To,Ti,KS,PS] = f_CLTFM(P_ss,K);
        
    case 'in_out'
        
        [Lo,Li,So,Si,To,Ti,KS,PS,Tniy,Tniu] = f_CLMapInnerOuter_BigK...
            (P_ss,K,Mi);
    
end


% *************************************************************************
%
% FORM STACKED CL MAPS (WITH WEIGHTING FUNCTIONS)
%

WTrz = [];
WTdz = [];
WTniz = [];


% Weighting function W1 (S_e).
if weights.exist_W1                 
    WTrz = [ WTrz ; series(So, weights.W1) ];
end

% Weighting function W2 (K*S_e).
if weights.exist_W2                 
    WTrz = [ WTrz ; series(KS, weights.W2) ];
end

% Weighting function W3 (T_e).
if weights.exist_W3                
    WTrz = [ WTrz ; series(To, weights.W3) ];
end

% Weighting function Wd1 (S_c).
if weights.exist_Wd1                
    WTdz = [ WTdz ; series(Si, weights.Wd1) ];
end

% Weighting function Wd2 (P*S_c).
if weights.exist_Wd2                
    WTdz = [ WTdz ; series(PS, weights.Wd2) ];
end

% Weighting function Wd3 (T_c).
if weights.exist_Wd3                
    WTdz = [ WTdz ; series(Ti, weights.Wd3) ];
end

% Weighting function Wni1 (S_ni).
if weights.exist_Wni1   
    error('NOT IMPLEMENTED!!!!')
    WTniz = [ WTniz ; series(Sni, weights.Wni1) ];
end

% Weighting function Wni2 (K*S_ni).
if weights.exist_Wni2   
    WTniz = [ WTniz ; series(Tniu, weights.Wni2) ];
end

% Weighting function Wni3 (T_ni).
if weights.exist_Wni2   
    WTniz = [ WTniz ; series(Tniy, weights.Wni3) ];
end

% *************************************************************************
%
% FIND OPTIMAL PERFORMANCE AT EACH LOOP BREAKING POINT
%

gamma_rz = norm(WTrz,inf);          % H-inf norm r -> z.
gamma_dz = norm(WTdz,inf);          % H-inf norm d -> z.
gamma_niz = norm(WTniz,inf);        % H-inf norm ni -> z.

% Optimal performance measure (max of r -> z, d -> z, ni -> z).
gamma_opt = max([gamma_rz, gamma_dz, gamma_niz]);



end
