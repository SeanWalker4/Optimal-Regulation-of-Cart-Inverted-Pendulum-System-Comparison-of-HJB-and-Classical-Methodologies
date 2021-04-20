% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% CART INVERTED PENDULUM MAIN PROGRAM
%
% Brent Wallace
% Sean Walker
%
% 3/18/2021
%
% Main program for cart inverted pendulum model (with translational
% damping, no rotational damping, pendulum mass concentrated at end). See
% Ogata, K. "Modern Conrol Engineering" 3rd ed. pp. 106.
%
%
% *************************************************************************
%                               PROGRAM FUNCTIONALITY:
% *************************************************************************
%
% This program runs designs for different control laws in the nonlinear
% simulation hosted on the simulink model "cip_sim.slx". The respective
% sections of the program can be run to yield the results for each
% controller.
%
% NOTE:
%
% In order to know which plant should be used in the simulation, the
% variable p_select (an integer) is used to select the proper plant. The
% proper value of p_select and its corresponding plant are listed below:
%
%       Linear Plant             1
%       Nonlinear Plant          2
%
% NOTE:
%
% In order to know which control law should be used in the simulation, the
% variable k_select (an integer) is used to select the proper controller.
% The proper value of k_select and its corresponding control law are listed
% below:
%
%       Linear HJB                  1
%       Nonlinear HJB               2
%       Classical GMS H-Infinity    3
%
% *************************************************************************
% *************************************************************************
% *************************************************************************


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% PROGRAM INITIALIZATION
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% ***********************
%
% FIGURES
%

% Close figures, clear command window.
clear
close all
clc

% Figure saving controls.
savefigs = 0;           % Save figures control.
relpath = 'figures/';  	% Relative file path for saved figures.
figcount = 1;           % Initialize figure counter.

% Create save directory if figures are to be saved.
if savefigs
    time = fix(clock);
    timestamp = '';
    for i = 1:length(time)
       timestamp = [timestamp , num2str(time(i))];
       if i < length(time)
           timestamp = [timestamp , '_'];
       end
    end
    timestamp = [timestamp, '\'];
    relpath = [relpath, timestamp];   % Update relative path.
    mkdir(relpath);                   % Create directory for relative path.
end


% ***********************
%
% ADD SIMULATOR
%
addpath('sim');

% ***********************
%
% ADD UTIL
%
addpath('util');

% ***********************
%
% ADD GMS OPTIMIZER
%
addpath('util/gms');


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% VARIABLE INITIALIZATION
% 
% *************************************************************************
% *************************************************************************
% *************************************************************************

% Laplace transform variable.
s = tf('s');

% Radians to degrees conversion.
r2d = 180/pi;  

% Frequency vector for frequency response plots.
wmin = 1e-4;
wmax = 1e4;
numwpts = 1000;
wvec = logspace(log10(wmin),log10(wmax),numwpts);

% x-axis step reference command.
r = 0;

% ***********************
%
% TYPE OF SIMULATION TO RUN
%
%   1 : Initial condition response.
%   2 : Disturbance response.
%

sim_type = 1;


% ***********************
%
% SIMULATION TYPE SPECIFIC SETTINGS
%
switch sim_type
    
    case 1      % Initial condition response.
        
        % Initial condition.
        %   STATE : [ x(m)  x'(m)  theta(rad)  theta'(rad/s)  ]

        x0 = [  0   
                0 
                30 / r2d 
                0 / r2d     ];

        
        % Force disturbance parameters.
        %
        %   dF(t) = { dF_amp,   0 <= t <= dF_tf
        %           { 0,        else

        dF_amp = 0;
        dF_tf = 0;
        
        % Plot x-axis limits
        tf_x = 30;
        tf_theta = 5;
        tf_u = 5;
        
        
    case 2      % Disturbance response.
        
        % Initial condition.
        %   STATE : [ x(m)  x'(m)  theta(rad)  theta'(rad/s)  ]

        x0 = zeros(4,1);
        
        % Force disturbance parameters.
        %
        %   dF(t) = { dF_amp,   0 <= t <= dF_tf
        %           { 0,        else

        dF_amp = 5;
        dF_tf = 0.25;
        
        % Plot x-axis limits
        tf_x = 30;
        tf_theta = 5;
        tf_u = 5;
    
    
end


% % % Plot x-axis limits
% %         tf_x = 50;
% %         tf_theta = 20;
% %         tf_u = 10;

% Length of simulation (sec).
t_f = max([tf_x tf_theta tf_u]);   % Final time (s).
% t_f = 8;                         % Final time (s).



% Plant selection.
%       Linear Plant             1
%       Nonlinear Plant          2

p_select = 2;

% Plot individual design results.
plot_individual = 0;

% Plot frequency responses (GMS only).
plot_freq_resp = 0;

% ***********************
%
% SIMULINK VARIANTS
%
% See:      https://www.youtube.com/watch?v=PE0IsVk3gy4
%
% Used for selecting between different subsystems (e.g., linear vs.
% nonlinear plant).
%

% Plant selection.
v_linear_plant = Simulink.Variant('p_select == 1');
v_nonlinear_plant = Simulink.Variant('p_select == 2');

% Controller selection.
v_lin_hjb = Simulink.Variant('k_select == 1');
v_nonlin_hjb = Simulink.Variant('k_select == 2');
v_gms = Simulink.Variant('k_select == 3');

%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% SYSTEM INITIALIZATION
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
% LINEARIZED SYSTEM STATE SPACE
%
%   STATE : [ x(m)  x'(m)  theta(rad)  theta'(rad/s)  ]^T
%   
%   INPUT : u(N)
%
%   OUTPUT: x(m)
% 

A = [   0       1       0               0
        0       -b/M    -(m*g/M)        0
        0       0       0               1
        0       0       (g/l)*(1+(m/M)) 0  ];

B = [   0
        1/M
        0
        -1/(M*l) ];
    
C = [ 1   0   0   0];

D = 0;

P_ss = ss(A,B,C,D);




%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
% 
% LQR
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

% *************************************************************************
% 
% DESIGN
%
% *************************************************************************

% ***********************
% 
% PARAMETER INITIALIZATION
%

% Stage-cost state penalty.
Q = diag([0.25;1;2;2]);
    
% Control penalty.
R = 0.2;


% ***********************
% 
% SOLVE THE ASSOCIATED ALGEBRAIC RICATTI EQUATION
%
%   0 = Q - K * B * R^{-1} * B^T * K + K * A + A^T * K
%

% ***********************
%
% ALGEBRAIC RICATTI EQUATION
%

K_icare = icare(A,B,Q,R,[],[],[]);

% ***********************
%
% OPTIMAL CONTROL LAW FULL STATE FEEDBACK MATRIX
%
%   u*(t) = - R^{-1} * B^T * K * x(t)
%

G_icare = - inv(R) * B' * K_icare;

% *************************************************************************
% 
% SIMULATION
%
% *************************************************************************

% ***********************
%
% SET SIMULATION STOP TIME
%

set_param('cip_sim', 'StopTime', num2str(t_f))


% ***********************
% 
% SELCT LINEAR HJB FOR SIMULATOR
%

k_select = 1;


% ***********************
% 
% RUN SIMULATION
%

simout = sim('cip_sim');


% ***********************
% 
% EXTRACT SIMULATION DATA
%

tvec_lin_hjb = simout.tout;                             % Time vector.
x_t_lin_hjb = simout.xvec.x.Data;                       % x(t).
xdot_t_lin_hjb = simout.xvec.xdot.Data;                 % x^{dot}(t).
theta_t_lin_hjb = r2d * simout.xvec.theta.Data;         % theta(t).
thetadot_t_lin_hjb = r2d * simout.xvec.thetadot.Data;   % theta^{dot}(t).
u_t_lin_hjb = simout.uvec.Data;                         % u(t).


% *************************************************************************
% 
% PLOT RESULTS
%
% *************************************************************************

if plot_individual

    % ***********************
    %
    % CART POSITION
    %

    figure(figcount)

    h_fig = plot(tvec_lin_hjb,x_t_lin_hjb);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Cart Position x(t)')
    xlabel('Time (s)')
    ylabel('x(t) (m)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['x_t_lin_hjb'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % CART VELOCITY
    %

    figure(figcount)

    h_fig = plot(tvec_lin_hjb,xdot_t_lin_hjb);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Cart Velocity x^{dot}(t)')
    xlabel('Time (s)')
    ylabel('x^{dot}(t) (m/s)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['xdot_t_lin_hjb'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % PENDULUM ANGLE
    %

    figure(figcount)

    h_fig = plot(tvec_lin_hjb, theta_t_lin_hjb);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Pendulum Angular Position \theta(t)')
    xlabel('Time (s)')
    ylabel('\theta(t) (deg)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['theta_t_lin_hjb'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % PENDULUM ANGLULAR RATE
    %

    figure(figcount)

    h_fig = plot(tvec_lin_hjb, thetadot_t_lin_hjb);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Pendulum Angular Rate \theta^{dot}(t)')
    xlabel('Time (s)')
    ylabel('\theta^{dot}(t) (deg/s)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['thetadot_t_lin_hjb'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % CONTROL SIGNAL u
    %

    figure(figcount)

    h_fig = plot(tvec_lin_hjb, u_t_lin_hjb);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Control Signal u(t)')
    xlabel('Time (s)')
    ylabel('u(t) (Newtons)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['u_t_lin_hjb'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.

end

%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
% 
% NONLINEAR HJB
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

% *************************************************************************
% 
% SIMULATION
%
% *************************************************************************

% ***********************
% 
% SELCT NONLINEAR HJB FOR SIMULATOR
%

k_select = 2;


% ***********************
% 
% RUN SIMULATION
%

simout = sim('cip_sim');


% ***********************
% 
% EXTRACT SIMULATION DATA
%

tvec_nln_hjb = simout.tout;                             % Time vector.
x_t_nln_hjb = simout.xvec.x.Data;                       % x(t).
xdot_t_nln_hjb = simout.xvec.xdot.Data;                 % x^{dot}(t).
theta_t_nln_hjb = r2d * simout.xvec.theta.Data;         % theta(t).
thetadot_t_nln_hjb = r2d * simout.xvec.thetadot.Data;   % theta^{dot}(t).
u_t_nln_hjb = simout.uvec.Data;                         % u(t).


% *************************************************************************
% 
% PLOT RESULTS
%
% *************************************************************************

if plot_individual

    % ***********************
    %
    % CART POSITION
    %

    figure(figcount)

    h_fig = plot(tvec_nln_hjb,x_t_nln_hjb);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Cart Position x(t)')
    xlabel('Time (s)')
    ylabel('x(t) (m)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['x_t_nln_hjb'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % CART VELOCITY
    %

    figure(figcount)

    h_fig = plot(tvec_nln_hjb,xdot_t_nln_hjb);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Cart Velocity x^{dot}(t)')
    xlabel('Time (s)')
    ylabel('x^{dot}(t) (m/s)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['xdot_t_nln_hjb'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % PENDULUM ANGLE
    %

    figure(figcount)

    h_fig = plot(tvec_nln_hjb, theta_t_nln_hjb);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Pendulum Angular Position \theta(t)')
    xlabel('Time (s)')
    ylabel('\theta(t) (deg)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['theta_t_nln_hjb'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % PENDULUM ANGLULAR RATE
    %

    figure(figcount)

    h_fig = plot(tvec_nln_hjb, thetadot_t_nln_hjb);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Pendulum Angular Rate \theta^{dot}(t)')
    xlabel('Time (s)')
    ylabel('\theta^{dot}(t) (deg/s)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['thetadot_t_nln_hjb'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % CONTROL SIGNAL u
    %

    figure(figcount)

    h_fig = plot(tvec_nln_hjb, u_t_nln_hjb);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Control Signal u(t)')
    xlabel('Time (s)')
    ylabel('u(t) (Newtons)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['u_t_nln_hjb'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.

end

%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
% 
% GMS
%
% *************************************************************************
% *************************************************************************
% *************************************************************************


% *************************************************************************
% 
% CONFIGURE OPTIMIZER
%
% *************************************************************************


% *************************************************************************
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
% 
% *************************************************************************

settings.loop_type = 'in_out';

settings.aug_integ = 0;

% % p1 = -0.001;
% p1 = -0.1;
% % p1 = -1;
% % p1 = -10;
% p2 = -1e20;

% settings.bilin_param.p1 = p1;
% settings.bilin_param.p2 = p2;

youla1_zames0 = 1;
settings.youla1_zames0 = youla1_zames0;

sys_struct.settings = settings;


% *************************************************************************
%
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
% 
% *************************************************************************


P.P_ss = P_ss;

Mi = [  0   0   1   0   ];      % Feed back x3 (theta) in inner loop.

P.Mi = Mi;

sys_struct.P = P;


% *************************************************************************
%
%       weights         : Closed-loop map weighting functions (state-space
%                           objects). NOTE: declare only if the
%                           respective closed-loop map is involved in the
%                           optimization.
%           W1          : Weighting function W1 (S_e).
%           W2          : Weighting function W2 (K*S_e).
%           W3          : Weighting function W3 (T_e).
%           Wd1         : Weighting function W4 (S_c).
%           Wd2         : Weighting function W5 (P*S_c).
%           Wd3         : Weighting function W6 (T_c).
%           Wni1        : Weighting function W7 (S_ni).
%           Wni2        : Weighting function W8 (K*S_ni).
%           Wni3        : Weighting function W9 (T_ni).
% 
% *************************************************************************


% *************************************************************************
%
% WEIGHTING W1
%
%   W1(s) = 1/Ms * (s+Ms*wb)/(s+eps*wb)
%
% Where wb > 0, 0 < eps << 1 < Ms. 
% 

eps1 = 0.1; 
Ms1 = 5; 
wb1 = 1;

W1 = 1/Ms1 * tf([1 Ms1*wb1],[1 eps1*wb1]);

W1 = ss(W1);


% *************************************************************************
%
% WEIGHTING W2
%
%   W2(s) = Mu
%
% Where Mu > 0. 
% 

Mu = 0.05;

W2 = Mu;

W2 = ss(W2);


% *************************************************************************
%
% WEIGHTING W3
%
%   W3(s) = (s + wbc/My) / wbc
%
% Where My > 0, wbc > 0. 
% 

My = 0.01;

% wbc = 1;
% 
% W3 = tf([1 wbc/My],[wbc]);
% 
% W3 = ss(W3);

W3 = ss(My);


% *************************************************************************
%
% WEIGHTING Wd1
%
%   Wd1(s) = 1/Ms * (s+Ms*wb)/(s+eps*wb)
%
% Where wb > 0, 0 < eps << 1 < Ms. 
% 

% epsd1 = 0.1; 
% Msd1 = 5; 
% wbd1 = 1;
% Wd1 = 1/Msd1 * tf([1 Msd1*wbd1],[1 epsd1*wbd1]);

Wd1 = 0.3;

Wd1 = ss(Wd1);

% Wd1 = W1;


% *************************************************************************
%
% WEIGHTING Wd2
%
%   Wd2(s) = Mu
%
% Where Mu > 0. 
% 

Mdu = 4;
Wd2 = Mdu;
Wd2 = ss(Wd2);

% Wd2 = W2;


% *************************************************************************
%
% WEIGHTING Wd3
%
%   Wd3(s) = (s + wbc/My) / wbc
%
% Where My > 0, wbc > 0. 
% 

Mdy = 0.1;

% wbcd = 1;
% 
% Wd3 = tf([1 wbcd/My],[wbcd]);
% 
% Wd3 = ss(Wd3);

Wd3 = ss(Mdy);


% *************************************************************************
%
% STORE WEIGHTINGS
%


if exist('W1','var')
    weights.W1 = W1;
end
if exist('W2','var')
    weights.W2 = W2;
end
if exist('W3','var')
    weights.W3 = W3;
end

if exist('Wd1','var')
    weights.Wd1 = Wd1;
end
if exist('Wd2','var')
    weights.Wd2 = Wd2;
end
if exist('Wd3','var')
    weights.Wd3 = Wd3;
end

sys_struct.weights = weights;

% *************************************************************************
%
% LEFT EMPTY: CONSTRAINT FUNCTIONALS
%
%       weights_c       : Closed-loop map constraint functionals. NOTE:
%                           declare only if the respective closed-loop map
%                           functional is involved in the optimization.
%                           Each is a struct with the following fields:
%               tfm     : Transfer function matrix of constraint functional
%                           (state space object)
%               Fun     : Functional type (string). Has following options:
%                           'f_Hinf' = H-Infty constraint functional.
%                           'f_Linf' = L-Infty constraint functional.
%               Val     : Functional constraint in optimization (double).                    
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
% *************************************************************************


% *************************************************************************
%
%   opt_struct          : Data structure with the following fields:
%       settings        : Optimization settings. Data structure with the
%                           following fields:
%           xmax1       : Upper bound on optimization solution (double).
%                           i.e., if nu is the number of controls, ny the
%                           number of measured signals, N the dimension of
%                           the basis, then nu*ne*N optimum xopt lies in
%                           the box [xmin1*ones(nu*ne*N),
%                           xmax1*ones(nu*ne*N)]
%           xmin1       : Lower bound on optimization solution (double).
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
%
% *************************************************************************


xmin1 = -1000;
xmax1 = 1000;
x01 = 1;

settings_opt.xmin1 = xmin1;
settings_opt.xmax1 = xmax1;
settings_opt.x01 = x01;

maxiter = 400;
sum1_max0 = 0;
algo = 2;

settings_opt.maxiter = maxiter;
settings_opt.sum1_max0 = sum1_max0;
settings_opt.algo = algo;

opt_struct.settings = settings_opt;


% *************************************************************************
%
%       basis           : Basis parameters. Data structure with the
%                           following fields:
%           type        : Type of basis used (integer). Has the following
%                           options:
%                           1 = fixed pole low pass
%                           2 = fixed pole all pass
%                           3 = variable pole low pass
%                           4 = variable pole all pass
%                           5 = pole-zero
%                           6 = Lauguerre
%           N           : Number of basis terms (integer).
%           p           : Pole location in basis (double). 
%           z           : Zero location in basis (double). 
%
% *************************************************************************

basis_type = 2;
N = 8;
p = 10;
z = 1;

basis.type = basis_type;
basis.N = N;
basis.p = p;
basis.z = z;

opt_struct.basis = basis;


% *************************************************************************
%
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
% *************************************************************************

mu1=1; 
mu2=1; 
mu3=1; 

obj_funct_w.mu1 = mu1;
obj_funct_w.mu2 = mu2;
obj_funct_w.mu3 = mu3;

rho1=1; 
rho2=1; 
rho3=1;

obj_funct_w.rho1 = rho1;
obj_funct_w.rho2 = rho2;
obj_funct_w.rho3 = rho3;


opt_struct.obj_funct_w = obj_funct_w;


% *************************************************************************
%
% CALL GMS OPTIMIZATION
% 

[K_gms, gamma_opt] = gms_main(sys_struct, opt_struct);

%%

% *************************************************************************
% 
% SIMULATION
%
% *************************************************************************

% ***********************
% 
% SELCT GMS FOR SIMULATOR
%

k_select = 3;

% ***********************
% 
% RUN SIMULATION
%

simout = sim('cip_sim');


% ***********************
% 
% EXTRACT SIMULATION DATA
%

tvec_gms = simout.tout;                             % Time vector.
x_t_gms = simout.xvec.x.Data;                       % x(t).
xdot_t_gms = simout.xvec.xdot.Data;                 % x^{dot}(t).
theta_t_gms = r2d * simout.xvec.theta.Data;         % theta(t).
thetadot_t_gms = r2d * simout.xvec.thetadot.Data;   % theta^{dot}(t).
u_t_gms = simout.uvec.Data;                         % u(t).


% *************************************************************************
% 
% PLOT RESULTS
%
% *************************************************************************

if plot_individual

    % ***********************
    %
    % CART POSITION
    %

    figure(figcount)

    h_fig = plot(tvec_gms,x_t_gms);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Cart Position x(t)')
    xlabel('Time (s)')
    ylabel('x(t) (m)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['x_t_gms'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % CART VELOCITY
    %

    figure(figcount)

    h_fig = plot(tvec_gms,xdot_t_gms);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Cart Velocity x^{dot}(t)')
    xlabel('Time (s)')
    ylabel('x^{dot}(t) (m/s)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['xdot_t_gms'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % PENDULUM ANGLE
    %

    figure(figcount)

    h_fig = plot(tvec_gms, theta_t_gms);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Pendulum Angular Position \theta(t)')
    xlabel('Time (s)')
    ylabel('\theta(t) (deg)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['theta_t_gms'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % PENDULUM ANGLULAR RATE
    %

    figure(figcount)

    h_fig = plot(tvec_gms, thetadot_t_gms);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Pendulum Angular Rate \theta^{dot}(t)')
    xlabel('Time (s)')
    ylabel('\theta^{dot}(t) (deg/s)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['thetadot_t_gms'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.


    % ***********************
    %
    % CONTROL SIGNAL u
    %

    figure(figcount)

    h_fig = plot(tvec_gms, u_t_gms);

    lgd = {};

    set(h_fig, 'linewidth', 2);
    title('Control Signal u(t)')
    xlabel('Time (s)')
    ylabel('u(t) (Newtons)')
    grid on
    % legend(lgd)

    % SAVE PLOT
    if savefigs
        filename = ['u_t_gms'];
        savepdf(figcount, relpath, filename); 
    end

    figcount = figcount + 1;        % Increment figure counter.

end

%%
% *************************************************************************
%
% FREQUENCY DOMAIN RESULTS (GMS ONLY)
% 
% *************************************************************************


if plot_freq_resp

% ***********************
%
% FORM CLOSED LOOP MAPS
% 


[Lo,Li,So,Si,To,Ti,KS,PS,Tniy,Tniu] = f_CLMapInnerOuter_BigK...
            (P_ss,K_gms,Mi);


        
% ***********************
%
% PLANT FREQUENCY RESPONSE
%

axes_vec_mag = [wmin,wmax,-40,50];      % Window to be used for magnitude
axes_vec_ph = [wmin,wmax,-200,0];      % Window to be used for phase

% See plotbode.m under the folder "SDToolbox" for documentation

sys_cell = {P_ss};
lgd_text = [];
wvec_cell = wvec;

plottypes = [1;1];  
ttl_cell = {['Plant P Magnitude Response'],['Plant P Phase Response']};
axes_cell = {axes_vec_mag, axes_vec_ph};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);

% SAVE FIGURE
filename = ['P_mag'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter

% SAVE FIGURE
filename = ['P_ph'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter


% ***********************
%
% CONTROLLER FREQUENCY RESPONSE
%

axes_vec_mag = [wmin,wmax,-40,50];      % Window to be used for magnitude
axes_vec_ph = [wmin,wmax,100,400];      % Window to be used for phase

% See plotbode.m under the folder "SDToolbox" for documentation

sys_cell = {K_gms};
lgd_text = [];
wvec_cell = wvec;

plottypes = [1;1];  
ttl_cell = {['Controller K Magnitude Response'],['Controller K Phase Response']};
axes_cell = {axes_vec_mag, axes_vec_ph};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);

% SAVE FIGURE
filename = ['K_mag'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter

% SAVE FIGURE
filename = ['K_ph'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter

% ***********************
%
% LOOP FREQUENCY RESPONSE
%

axes_vec_mag = [wmin,wmax,-40,40];      % Window to be used for magnitude
axes_vec_ph = [wmin,wmax,-50,190];      % Window to be used for phase

% See plotbode.m under the folder "SDToolbox" for documentation

sys_cell = {Lo};
lgd_text = [];
wvec_cell = wvec;

plottypes = [1;1];  
ttl_cell = {['Loop L Magnitude Response'],['Loop L Phase Response']};
axes_cell = {axes_vec_mag, axes_vec_ph};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);

% SAVE FIGURE
filename = ['L_mag'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter

% SAVE FIGURE
filename = ['L_ph'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter


% ***********************
%
% SENSITIVITY AT OUTPUT, gamma * |W_1^(-1)|
%

axes_vec_mag = [wmin,wmax,-40,100];      % Window to be used for magnitude

% See plotbode.m under the folder "SDToolbox" for documentation

sys_cell = {So, gamma_opt * inv(W1)};
lgd_text = {'S_e', '\gamma_{opt} |W_1^{-1}|'};
wvec_cell = wvec;

plottypes = [1;0];  
ttl_cell = {['Sensitivity at Error S_e']};
axes_cell = {axes_vec_mag};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);

% SAVE FIGURE
filename = ['S_e'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter

% ***********************
%
% COMPLEMENTARY SENSITIVITY AT OUTPUT, gamma * |W_3^(-1)|
%

axes_vec_mag = [wmin,wmax,-40,100];      % Window to be used for magnitude

% See plotbode.m under the folder "SDToolbox" for documentation
sys_cell = {To, gamma_opt * inv(W3)};
lgd_text = {'T_e', '\gamma_{opt} |W_3^{-1}|'};
wvec_cell = wvec;

plottypes = [1;0];  
ttl_cell = {['Complementary Sensitivity at Error T_e']};
axes_cell = {axes_vec_mag};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);


% SAVE FIGURE
filename = ['T_e'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter


% ***********************
%
% SENSITIVITY AND COMPLEMENTARY SENSITIVITY AT OUTPUT
%

axes_vec_mag = [wmin,wmax,-40,20];      % Window to be used for magnitude

% See plotbode.m under the folder "SDToolbox" for documentation
sys_cell = {So, To};
lgd_text = {'S_e', 'T_e'};
wvec_cell = wvec;

plottypes = [1;0];  
ttl_cell = {['Sensitvity S_e and Complementary Sensitivity T_e at Error']};
axes_cell = {axes_vec_mag};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);


% SAVE FIGURE
filename = ['S_e_and_T_e'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter


% ***********************
%
% REFERENCE TO CONTROL Tru, gamma * |W_2^(-1)|
%

axes_vec_mag = [wmin,wmax,-40,100];      % Window to be used for magnitude

% See plotbode.m under the folder "SDToolbox" for documentation
sys_cell = {KS, gamma_opt * inv(W2)};
lgd_text = {'KS_e', '\gamma_{opt} |W_2^{-1}|'};
wvec_cell = wvec;

plottypes = [1;0];  
ttl_cell = {['Reference to Control: T_{ru} = KS_e']};
axes_cell = {axes_vec_mag};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);

% SAVE FIGURE
filename = ['Tru'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter


% ***********************
%
% SENSITIVITY AT INPUT, gamma * |W_d1^(-1)|
%

axes_vec_mag = [wmin,wmax,-40,100];      % Window to be used for magnitude

% See plotbode.m under the folder "SDToolbox" for documentation

sys_cell = {Si, gamma_opt * inv(Wd1)};
lgd_text = {'S_u', '\gamma_{opt} |W_4^{-1}|'};
wvec_cell = wvec;

plottypes = [1;0];  
ttl_cell = {['Sensitivity at Control S_u']};
axes_cell = {axes_vec_mag};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);

% SAVE FIGURE
filename = ['S_c'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter

% ***********************
%
% COMPLEMENTARY SENSITIVITY AT INPUT, gamma * |W_d3^(-1)|
%

axes_vec_mag = [wmin,wmax,-40,100];      % Window to be used for magnitude

% See plotbode.m under the folder "SDToolbox" for documentation
sys_cell = {Ti, gamma_opt * inv(Wd3)};
lgd_text = {'T_u', '\gamma_{opt} |W_6^{-1}|'};
wvec_cell = wvec;

plottypes = [1;0];  
ttl_cell = {['Complementary Sensitivity at Control T_u']};
axes_cell = {axes_vec_mag};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);


% SAVE FIGURE
filename = ['T_c'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter


% ***********************
%
% SENSITIVITY AND COMPLEMENTARY SENSITIVITY AT INPUT
%

axes_vec_mag = [wmin,wmax,-40,20];      % Window to be used for magnitude

% See plotbode.m under the folder "SDToolbox" for documentation
sys_cell = {Si, Ti};
lgd_text = {'S_u', 'T_u'};
wvec_cell = wvec;

plottypes = [1;0];  
ttl_cell = {['Sensitvity S_u and Complementary Sensitivity T_u at Control']};
axes_cell = {axes_vec_mag};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);


% SAVE FIGURE
filename = ['S_c_and_T_c'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter



% ***********************
%
% INPUT DISTURBANCE TO OUTPUT Tdiy, gamma * |W_d2^(-1)|
%

axes_vec_mag = [wmin,wmax,-40,100];      % Window to be used for magnitude

% See plotbode.m under the folder "SDToolbox" for documentation
sys_cell = {PS, gamma_opt * inv(Wd2)};
lgd_text = {'PS_c', '\gamma_{opt} |W_5^{-1}|'};
wvec_cell = wvec;

plottypes = [1;0];  
ttl_cell = {['Input Disturbance to Plant Output: T_{d_iy} = PS_c']};
axes_cell = {axes_vec_mag};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);

% SAVE FIGURE
filename = ['Tdiy'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter


% ***********************
%
% INNER-LOOP DISTURBANCE TO CONTROL Tniu
%

axes_vec_mag = [wmin,wmax,-40,100];      % Window to be used for magnitude

% See plotbode.m under the folder "SDToolbox" for documentation
sys_cell = {Tniu};
lgd_text = {};
wvec_cell = wvec;

plottypes = [1;0];  
ttl_cell = {['Inner-Loop Disturbance to Control: T_{n_iu}']};
axes_cell = {axes_vec_mag};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);

% SAVE FIGURE
filename = ['Tniu'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter


% ***********************
%
% INNER-LOOP DISTURBANCE TO PLANT OUTPUT Tniy
%

axes_vec_mag = [wmin,wmax,-40,100];      % Window to be used for magnitude

% See plotbode.m under the folder "SDToolbox" for documentation
sys_cell = {Tniy};
lgd_text = {};
wvec_cell = wvec;

plottypes = [1;0];  
ttl_cell = {['Inner-Loop Disturbance to Output: T_{n_iy}']};
axes_cell = {axes_vec_mag};

plotbode(sys_cell, wvec_cell, plottypes, ttl_cell, lgd_text, axes_cell, figcount);

% SAVE FIGURE
filename = ['Tniy'];
if savefigs
    savepdf(figcount,relpath,filename);
end

figcount = figcount + 1;        % Increment figure counter

end


%%
% *************************************************************************
% *************************************************************************
% *************************************************************************
% 
% PLOT RESULTS FOR ALL 3 DESIGNS
%
% *************************************************************************
% *************************************************************************
% *************************************************************************

% ***********************
%
% LEGEND
%
lgd = {'LQR', 'HJB', 'GMS'};


% ***********************
%
% CART POSITION
%

figure(figcount)

h_fig = plot(tvec_lin_hjb,x_t_lin_hjb, ...
        tvec_nln_hjb,x_t_nln_hjb, tvec_gms,x_t_gms);

set(h_fig, 'linewidth', 2);
title('Cart Position x(t)')
xlabel('Time (s)')
ylabel('x(t) (m)')
grid on
xlim([0,tf_x])
legend(lgd)

% SAVE PLOT
if savefigs
    filename = ['x_t_comp'];
    savepdf(figcount, relpath, filename); 
end

figcount = figcount + 1;        % Increment figure counter.


% ***********************
%
% CART VELOCITY
%

figure(figcount)

h_fig = plot(tvec_lin_hjb,xdot_t_lin_hjb, ...
        tvec_nln_hjb,xdot_t_nln_hjb, tvec_gms,xdot_t_gms);

set(h_fig, 'linewidth', 2);
title('Cart Velocity x^{dot}(t)')
xlabel('Time (s)')
ylabel('x^{dot}(t) (m/s)')
grid on
legend(lgd)

% SAVE PLOT
if savefigs
    filename = ['xdot_t_comp'];
    savepdf(figcount, relpath, filename); 
end

figcount = figcount + 1;        % Increment figure counter.


% ***********************
%
% PENDULUM ANGLE
%

figure(figcount)

h_fig = plot(tvec_lin_hjb,theta_t_lin_hjb, ...
        tvec_nln_hjb,theta_t_nln_hjb, tvec_gms,theta_t_gms);

set(h_fig, 'linewidth', 2);
title('Pendulum Angular Position \theta(t)')
xlabel('Time (s)')
ylabel('\theta(t) (deg)')
grid on
xlim([0,tf_theta])
legend(lgd)

% SAVE PLOT
if savefigs
    filename = ['theta_t_comp'];
    savepdf(figcount, relpath, filename); 
end

figcount = figcount + 1;        % Increment figure counter.


% ***********************
%
% PENDULUM ANGLULAR RATE
%

figure(figcount)

h_fig = plot(tvec_lin_hjb,thetadot_t_lin_hjb, ...
        tvec_nln_hjb,thetadot_t_nln_hjb, tvec_gms,thetadot_t_gms);


set(h_fig, 'linewidth', 2);
title('Pendulum Angular Rate \theta^{dot}(t)')
xlabel('Time (s)')
ylabel('\theta^{dot}(t) (deg/s)')
grid on
legend(lgd)

% SAVE PLOT
if savefigs
    filename = ['thetadot_t_comp'];
    savepdf(figcount, relpath, filename); 
end

figcount = figcount + 1;        % Increment figure counter.


% ***********************
%
% CONTROL SIGNAL u
%

figure(figcount)

h_fig = plot(tvec_lin_hjb,u_t_lin_hjb, ...
        tvec_nln_hjb,u_t_nln_hjb, tvec_gms,u_t_gms);

set(h_fig, 'linewidth', 2);
title('Control Signal u(t)')
xlabel('Time (s)')
ylabel('u(t) (Newtons)')
grid on
xlim([0,tf_u])
legend(lgd)

% SAVE PLOT
if savefigs
    filename = ['u_t_comp'];
    savepdf(figcount, relpath, filename); 
end

figcount = figcount + 1;        % Increment figure counter.