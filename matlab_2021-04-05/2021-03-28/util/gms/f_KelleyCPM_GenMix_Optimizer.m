function [xk, frz, fdz, fnz] = f_KelleyCPM_GenMix_Optimizer ...
    (dim, N, x0, Mobjrz, Mobjdz, Mobjniz, Mconrz, Mcondz, Mconniz, ...
    T11rz, T12rz, T21rz, T11dz, T12dz, T21dz, T11niz, T12niz, T21niz, ...
    Datarz, Datadz, Dataniz, Q, q, n_u, n_e, MaxIter, xmax, xmin)
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% KELLY'S CUTTING PLANE METHOD (CPM) ALGORITHM    
%
% Karan Puttannaiah
% Brent Wallace  
%
% 12/14/2020
%
% This program implements the Kelly's CPM algorithm for the convex program:
%
%             min           f_0(x)
%           x \in R^n
%
%             s.t.
%                           f_i(x) <= 0,   i = 1,...,m
%
% Where f_i: R^n -> R, i = 0,...,m, are convex and subdifferentiable.
%                           
%
% ***** CALL SYNTAX:
%
% [xk, frz, fdz, fnz] = f_KelleyCPM_GenMix_Optimizer ...
%     (dim, N, x0, Mobjrz, Mobjdz, Mobjniz, Mconrz, Mcondz, Mconniz, ...
%     T11rz, T12rz, T21rz, T11dz, T12dz, T21dz, T11niz, T12niz, T21niz, ...
%     Datarz, Datadz, Dataniz, Q, q, n_u, n_e, MaxIter, xmax, xmin)
%
% ***** INPUTS:
%
%   dim                 : Dimension of the program decision variable space.
%   N                   : Dimension of Q parameter basis.
%   x0                  : Feasible initial condition (N*n_u*n_e vector).
%   Mobjrz              : Vectorized r -> z Q parameterization transfer
%                           function matrices. See f_Vectorize.m for
%                           further documentation.
%   Mobjdz              : Vectorized d -> z Q parameterization transfer
%                           function matrices. See f_Vectorize.m for
%                           further documentation.
%   Mobjniz             : Vectorized ni -> z Q parameterization transfer
%                           function matrices. See f_Vectorize.m for
%                           further documentation.
%   Mconrz              : Vectorized r -> z Q parameterization transfer
%                           function matrices for constraint functionals.
%                           See f_Vectorize.m for further documentation.
%   Mcondz              : Vectorized d -> z Q parameterization transfer
%                           function matrices for constraint functionals.
%                           See f_Vectorize.m for further documentation.
%   Mconniz              : Vectorized ni -> z Q parameterization transfer
%                           function matrices for constraint functionals.
%                           See f_Vectorize.m for further documentation.
%   T11rz               : (1,1) transfer function matrix of r -> z Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   T12rz               : (1,2) transfer function matrix of r -> z Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   T21rz               : (2,1) transfer function matrix of r -> z Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   T11dz               : (1,1) transfer function matrix of d -> z Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   T12dz               : (1,2) transfer function matrix of d -> z Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   T21dz               : (2,1) transfer function matrix of d -> z Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   T11niz              : (1,1) transfer function matrix of ni -> z Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   T12niz              : (1,2) transfer function matrix of ni -> z Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   T21niz              : (2,1) transfer function matrix of ni -> z Q
%                           parameterization. See f_CoprFac.m for further
%                           documentation.
%   Datarz              : r -> z relevant optimization data (data
%                           structure). See f_GenData.m for further
%                           documentation.
%   Datadz              : d -> z relevant optimization data (data
%                           structure). See f_GenData.m for further
%                           documentation.
%   Dataniz             : ni -> z relevant optimization data (data
%                           structure). See f_GenData.m for further
%                           documentation.
%   Q                   : Initial Q parameter formed from basis qk and
%                           initial condition vector x0.
%   q                   : 1 by N cell array of zpk objects holding basis
%                           transfer function terms. 
%   n_u                 : Number of control signals (integer).
%   n_e                 : Number of regulated signals (integer).
%   MaxItern            : Maximum number of iterations before algorithm
%                           times out (integer).
%   xmax                : N*n_u*n_e upper bound vector for optimum. Optimum
%                           must lie in the box [xmin, xmax].
%   xmin                : N*n_u*n_e lower bound vector for optimum. Optimum
%                           must lie in the box [xmin, xmax].
%
% ***** OUTPUTS:
%
%   xk                  : Suboptimal solution (within tolerance). N*n_u*n_e
%                           dimensional vector.
%   frz                 : mu * ||Trz(K)||_{Hoo} optimal value (double).
%   fdz                 : rho * ||Tdz(K)||_{Hoo} optimal value (double).
%   fnz                 : eta * ||Tniz(K)||_{Hoo} optimal value (double).
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
% DETERMINE WHICH LOOP BREAKING POINTS ARE INVOLVED IN THE OPTIMIZATION
%

do_rz = ~ isempty(T11rz);       % r -> z
do_dz = ~ isempty(T11dz);       % d -> z
do_nz = ~ isempty(T11niz);      % ni -> z

% ***********************
%
% OBJECTIVE FUNCTION AND FEASIBLE SPACE TOLERANCES
%
tol_obj = 1e-4;
tol_feas = 1e-4;

% Turn of warnings.
warning off; 

% Display off for linprog.m.
options = optimset('Display','off'); 

% ***********************
%
% INITIALIZE ALGORITHM PARAMETERS
%

iter = 0;     	% Iteration count.
xk  = x0;    	% Initial query point.
cont = 1;       % Boolean variable which controls algorithm loop.

% % Initialize objective function oputputs.
frz{1} = -1;
fdz{1} = -1;
fnz{1} = -1;

% Number of constraints.
nConrz = Datarz.ConNum; 
nCondz = Datadz.ConNum;
nConnz = Dataniz.ConNum;

% Below matrices are used in solving the LP: min c'x   s.t. Aw<b.
Ao = [];                    % A matrix associated with objective function.
bo = [];                    % b vector associated with objective function.
c = [zeros(dim,1); 1];      % cvector associated with the variable x.

% Initialize optimal objective function lower bound/upper bound estimate
% differences for each loop breaking point. If the current loop breaking
% point is not involved in the optimization, make the difference negative
% to make the condition check inactive.

% r -> z
if do_rz                
    UkminLkrz = 1000;
else
    UkminLkrz = -1;
end

% d -> z
if do_dz                
    UkminLkdz = 1000;
else
    UkminLkdz = -1;
end

% ni -> z
if do_nz                
    UkminLknz = 1000;
else
    UkminLknz = -1;
end

% Set flag for constrained optimization at each loop breaking point.
constraint_flagrz = 1; 
constraint_flagdz = 1;
constraint_flagnz = 1;

% Initialize empty N*n_u*n_e + 1 dimensional vector which will hold optimum
% for the associated LP.
w = zeros(dim+1,1);

% ***********************
%
% DATA STORAGE
%

% Stores xk at each iteration.
xkStore = NaN*ones(length(xk),MaxIter); 
  
% Stores exit flag at each iteration.
ExitFlagStore = NaN*ones(1,MaxIter);   
   
% Stores objective function value at each iteration.
foStore = NaN*ones(3,MaxIter);


% *************************************************************************
% *************************************************************************
%
% BEGIN ALGORITHM
%
% *************************************************************************
% *************************************************************************


while cont
    
    % For constraints. Initialized empty each iteration, then appended to
    % main LP matrices Ao and bo at the end of the iteration.
    Ac=[]; bc=[];

    % *********************************************************************
    %
    % r -> z
    %

    if do_rz

    % ***********************
    %
    % WEIGHTED MIXED SENSITIVITY r -> z
    %

    % Evaluate r -> z objective function (H-Infty norm) f_{0,rz}(xk) and
    % subgradient G_{0,rz}(xk) at current query point xk.
    [forz, Gfo] = ...
        feval('f_Hinf', Mobjrz, xk, T11rz, T12rz, T21rz, Q, Datarz.ObjVec);
    
    % If U_{k,rz} - L_{k,rz} > tol_obj, then add to relaxed LP.
    if UkminLkrz > tol_obj
        Ao = [Ao; Gfo' -1];
        bo = [bo; Gfo'*xk-forz];
    end
    
    % ***********************
    %
    % CONSTRAINTS r -> z
    %

    % Evaluate r -> z constraint functionals (H-Infty or L-Infty norm)
    % f_{i,rz}(xk) and subgradient G_{i,rz}(xk) at current query point xk.
    frz{1}=[];
    for ii = 1:nConrz
        Mrz = Mconrz(ii,:);
        [frz{ii}, Gf_ii, ConValVec] = ...
            feval(Datarz.ConNam{ii}, Mrz, xk, T11rz, T12rz, T21rz, Q, ...
                Datarz.ConVec{ii}, Datarz.ConVal{ii});
        % In AllStep this is changed. 
        % Originally it was frz{ii} - Datarz.ConVal{ii}
        frz{ii} = frz{ii} - ConValVec'; 
        
        % If f_{i,rz}(xk) < tol_feas for some i, then add to relaxed LP.
        if constraint_flagrz > 0
            Ac = [Ac; Gf_ii' zeros(size(Gf_ii',1),1)];
            bc = [bc; Gf_ii'*xk-frz{ii}];
        end
    end

    end

    % *********************************************************************
    %
    % d -> z
    %

    if do_dz

    % ***********************
    %
    % WEIGHTED MIXED SENSITIVITY d -> z
    %

    % Evaluate d -> z objective function (H-Infty norm) f_{0,dz}(xk) and
    % subgradient G_{0,dz}(xk) at current query point xk.
    [fodz, Gfo] = ...
        feval('f_Hinf', Mobjdz, xk, T11dz, T12dz, T21dz, Q, Datadz.ObjVec);
    
    % If U_{k,dz} - L_{k,dz} > tol_obj, then add to relaxed LP.
    if UkminLkdz > tol_obj
        Ao = [Ao; Gfo' -1];
        bo = [bo; Gfo'*xk-fodz];
    end
    
    % ***********************
    %
    % CONSTRAINTS d -> z
    %

    % Evaluate d -> z constraint functionals (H-Infty or L-Infty norm)
    % f_{i,dz}(xk) and subgradient G_{i,dz}(xk) at current query point xk.
    fdz{1}=[];
    for ii = 1:nCondz
        Mdz = Mcondz(ii,:);
        [fdz{ii}, Gf_ii, ConValVec] = ...
            feval(Datadz.ConNam{ii}, Mdz, xk, T11dz, T12dz, T21dz, Q, ...
                Datadz.ConVec{ii}, Datadz.ConVal{ii});
        % In AllStep this is changed. 
        % Originally it was fdz{ii} - Datadz.ConVal{ii}
        fdz{ii} = fdz{ii} - ConValVec'; 
        
        % If f_{i,dz}(xk) < tol_feas for some i, then add to relaxed LP.
        if constraint_flagrz > 0
            Ac = [Ac; Gf_ii' zeros(size(Gf_ii',1),1)];
            bc = [bc; Gf_ii'*xk-fdz{ii}];
        end
    end
    
    end
   
    % *********************************************************************
    %
    % ni -> z
    %

    if do_nz

    % ***********************
    %
    % WEIGHTED MIXED SENSITIVITY ni -> z
    %

    % Evaluate ni -> z objective function (H-Infty norm) f_{0,niz}(xk) and
    % subgradient G_{0,niz}(xk) at current query point xk.
    [fonz, Gfo] = ...
       feval('f_Hinf', Mobjniz, xk, T11niz, T12niz, T21niz, ...
                Q, Dataniz.ObjVec);
    
    % If U_{k,niz} - L_{k,niz} > tol_obj, then add to relaxed LP.
    if UkminLknz > tol_obj
        Ao = [Ao; Gfo' -1];
        bo = [bo; Gfo'*xk-fonz];
    end
    
    % ***********************
    %
    % CONSTRAINTS ni -> z
    %

    % Evaluate ni -> z constraint functionals (H-Infty or L-Infty norm)
    % f_{i,niz}(xk) and subgradient G_{i,niz}(xk) at current query point
    % xk.
    fnz{1}=[];
    for ii = 1:nConnz
        Mnz = Mconniz(ii,:);
        [fnz{ii}, Gf_ii, ConValVec] = ...
            feval(Dataniz.ConNam{ii}, Mnz, xk, ...
                T11niz, T12niz, T21niz, Q, ...
                Dataniz.ConVec{ii}, Dataniz.ConVal{ii});
        % In AllStep this is changed. 
        % Originally it was fniz{ii} - Dataniz.ConVal{ii}
        fnz{ii} = fnz{ii} - ConValVec'; 
        
        % If f_{i,niz}(xk) < tol_feas for some i, then add to relaxed LP.
        if constraint_flagnz > 0
            Ac = [Ac; Gf_ii' zeros(size(Gf_ii',1),1)];
            bc = [bc; Gf_ii'*xk-fnz{ii}];
        end
    end

    end

    % *********************************************************************
    %
    % SOLVE ASSOCIATED LP
    %

    % Append mixed sensitivity LP matrices to constraint functional LP
    % matrices.
    Ao = [  Ao
            Ac  ]; 
    bo = [  bo
            bc  ];
    
    % Solve LP (used optimization toolbox function: linprog)
    [w,fval,exitflag] = linprog(c,Ao,bo,[],[],xmin,xmax,xk,options);

    % Check if problem is giving empty w. Try using other algorithms
    optionstemp = options; % temporary option
    if exitflag == -4
        optionstemp.Algorithm='dual-simplex';
        [w,fval,exitflag] = linprog(c,Ao,bo,[],[],xmin,xmax,xk,optionstemp);
    end
    if exitflag == -4
        optionstemp.Algorithm='active-set';
        [w,fval,exitflag] = linprog(c,Ao,bo,[],[],xmin,xmax,xk,optionstemp);
    end
    

    % *********************************************************************
    %
    % PREPARE FOR NEXT ITERATION
    %

    % ***********************
    %
    % EVALUATE UPPER BOUND / LOWER BOUND DIFFERENCE U_k - L_k FOR EACH LOOP
    % BREAKING POINT
    %

    if do_rz
        UkminLkrz = forz - c'*w;            % r -> z.
    end
    if do_dz
        UkminLkdz = fodz - c'*w;            % d -> z.
    end
    if do_nz
        UkminLknz = fonz - c'*w;            % ni -> z.
    end

    % ***********************
    %
    % VARIABLE UPDATES
    %

    % Increment iteration counter.
    iter = iter + 1;

    % Update xk.
    xk = w(1:dim); 

    % Update Q parameter QN(xk).
    Q = f_FormQN(xk, q, n_u, n_e, N); 

    % ***********************
    %
    % DATA STORAGE
    %
   
    % Store objective function values at the current iteration.
    if do_rz
        foStore(1,iter) = forz;            % r -> z.
    end
    if do_dz
        foStore(2,iter) = fodz;            % d -> z.
    end
    if do_nz
        foStore(3,iter) = fonz;            % ni -> z.
    end

    % Store xk.
    xkStore(:,iter) = xk; 

    % Store exit flag.
    ExitFlagStore(1,iter) = exitflag;
    
    % ***********************
    %
    % CHECK IF CONSTRAINT FUNCTION TOLERANCES ARE MET
    %
    % Check if fi(xk) < eps_feas for all i.
    %

    % r -> z
    constraint_flagrz = 0;
    for ii = 1:nConrz
        if frz{ii} > tol_feas
            constraint_flagrz = 1;
        end
    end

    % d -> z
    constraint_flagdz = 0;
    for ii = 1:nCondz
        if fdz{ii} > tol_feas
            constraint_flagdz = 1;
        end
    end

    % ni -> z
    constraint_flagnz = 0;
    for ii = 1:nConnz
        if fnz{ii} > tol_feas
            constraint_flagnz = 1;
        end
    end
    

    % ***********************
    %
    % DETERMINE IF LOOP SHOULD CONTINUE TO NEXT ITERATION
    % 
    % Loop must continue if U_k - L_k > eps_obj or fi(xk) > eps_feas for
    % any i, at any loop breaking point.
    %
    cont = (UkminLkrz > tol_obj) || (UkminLkdz > tol_obj) || ...
            (UkminLknz > tol_obj) || constraint_flagrz || ...
            constraint_flagdz || constraint_flagnz;

    % ***********************
    %
    % DISPLAY OUTPUT DATA FOR THIS ITERATION
    %
    prntstrng = [ 'ITERATION: ' num2str(iter) ];

    if do_rz
        prntstrng = [ prntstrng  ...
            '  ///  mu * ||Trz(K)||_{Hoo} = ' num2str(forz,6) ...
            '  ///  U_{rz} - L{rz} = ' num2str(UkminLkrz,6) ]; 
    end
    if do_dz
        prntstrng = [ prntstrng  ...
            '  ///  rho * ||Tdz(K)||_{Hoo} = ' num2str(fodz,6) ...
            '  ///  U_{dz} - L{dz} = ' num2str(UkminLkdz,6) ]; 
    end
    if do_nz
        prntstrng = [ prntstrng  ...
            '  ///  eta * ||Tniz(K)||_{Hoo} = ' num2str(fonz,6) ...
            '  ///  U_{niz} - L{niz} = ' num2str(UkminLknz,6) ];
    end

    disp(prntstrng)     % Display string in command window.

    % ***********************
    %
    % IF MAX NUMBER OF ITERATIONS REACHED, THEN BREAK AND END PROGRAM
    %
    if iter == MaxIter
        
        % Display error.
        disp(' ***** ERROR: Max number of iterations exceeded.')
        
        % Break and end program.
        break;

    end

end