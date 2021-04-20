% *************************************************************************
% *************************************************************************
% *************************************************************************
%
% EEE 587: OPTIMAL CONTROL
% 
% HOMEWORK ASSIGNMENT 4
%
% Brent Wallace
%
% 2021/03/08
%
%   
% *************************************************************************
% *************************************************************************
% *************************************************************************

%%
% *************************************************************************
% *************************************************************************
%
% INITIALIZATION
% 
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
savefigs = 1;           % Save figures control.
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




%%
% *************************************************************************
% *************************************************************************
%
% INITIALIZATION
% 
% *************************************************************************
% *************************************************************************


% *************************************************************************
%
% VARIABLE INITIALIZATION
%

% Time.
t0 = 0;             % Initial time.
tf = 10;            % Final time.
dT = 0.02;          % Time step. NOTE: Make sure (tf-t0)/dT is an integer!
N = (tf-t0) / dT + 1;     % Number of stages.

tvec_pretrans = 0:dT:tf-t0;      % Time vector (before translation).
tvec = t0:dT:tf;                 % Time vector (after translation).

% Initial state.
x_0 = [1 ; -1];    


% *************************************************************************
%
% INITIALIZE OPTIMIZATION PARAMETERS
% 

% Terminal cost state penalty.
H = [   20   0
        0   0   ];

% Stage-cost state penalty.
Q =  [   1   0 
         0   2   ];
    
% Control penalty.
R = 1;

% *************************************************************************
%
% INITIALIZE SYSTEM
% 


A = [   0   1
        -1   -2   ];

B = [   0
        1   ];  
    
C = [ 0   0 ];

D = 0;

sys_ct = ss(A, B, C, D);

% ***********************
%
% ZERO ORDER HOLD DISCRETIZATION
%

sys_zoh = ss(c2d(sys_ct, dT));

A_dt = sys_zoh.a;
B_dt = sys_zoh.b;

% Dimension of state vector.
n = size(A,1);


%%
% *************************************************************************
% *************************************************************************
%
% SOLVE DIFFERENTIAL RICATTI EQUATION
%
% To solve the differential Riccati equation with final value:
%
% 0 = K^{dot}(t) + Q - K(t) * B * R^{-1} * B^T * K(t) + K(t) * A + A^T *
% K(t),
%   t_0 <= t <= t_f,   K(t_f) = H.
%
% We solve the following differential Riccati equation with initial value:
%
% 0 = X^{dot}(t) - Q + X(t) * B * R^{-1} * B^T * X(t) - X(t) * A - A^T *
% X(t),
%   0 <= t <= t_f - t_0,   X(0) = H.
%
% Then 
%           K(t) = X( - (t - (t_f-t0)) ),   0 <= t <= t_f - t_0 
%
% satisfies the differential Riccati equation with final value specified
% above.
%
% *************************************************************************
% *************************************************************************

% *************************************************************************
%
% SOLVE:
%
% 0 = X^{dot}(t) - Q + X(t) * B * R^{-1} * B^T * X(t) - X(t) * A - A^T *
% X(t),
%   0 <= t <= t_f - t_0,   X(0) = H.
%
% NOTE:
%
% ode45 returns X(t) as an N x n^2 matrix, where N is the length of the
% time vector and n is the dimension of the state vector.
%

% % VARIABLE LENGTH TIME STEP
% [T, X] = ode45(@(t,X)mric(t, X, A, B, Q, R), [0 tf-t0], H);

% FIXED LENGTH TIME STEP
[T, X] = ode45(@(t,X)mric(t, X, A, B, Q, R), tvec_pretrans, H);


% *************************************************************************
%
% RESHAPE X(t):
%
% ode45 returns X(t) as an N x n^2 matrix, where N is the length of the
% time vector and n is the dimension of the state vector. We want X(t) as
% an N x 1 cell array of n x n matrices. Thus, for each k = 1,...,N, we
% make the following assignment:
%
%   X_k = reshape(X(k,:),size(A))
%   X_cell
%

% N = length(T);          % Length of time vector.

% N x 1 cell array, the k-th entry of which will contain the n x n matrix
% X(k).
X_cell = cell(N,1);     

for k = 1:N
    
    % Reshape X(k) into the proper n x n dimensions
    X_k = reshape(X(k,:), size(A));
    
    % Store X(k)
    X_cell{k} = X_k;
    
end


% *************************************************************************
%
% DEFINE K(t):
%
%           K(t) = X( - (t - (t_f-t0)) ) 
%

% % % % ***********************
% % % %
% % % % TIME VECTOR CORRESPONDING TO K(t)
% % % %
% % % % For 0 <= T <= t_f - t_0,
% % % %
% % % % t =  - (T - (t_f-t_0)) 
% % % %
% % % 
% % % tvec = - (T - (tf-t0));

% ***********************
%
% GET K(t) FROM X(t)
%
% Time indexing is reversed, so to get K_cell, invert the ordering of
% X_cell.
%

% N x 1 cell array, the k-th entry of which will contain the n x n matrix
% K(k).
K_cell = cell(N,1);   

% N x n^2 matrix, the k-th row of which will contain the vectorized version
% of K(k)
Kmat = zeros(N,n^2);

for k = 1:N
    
    % Store K(k).
    K_cell{k} = X_cell{N - k + 1};
    
    % Store vectorized K(k).
    Kmat(N,:) = X(N - k + 1, :);
    
end

%%
% *************************************************************************
% *************************************************************************
%
% SOLVE THE ASSOCIATED ALGEBRAIC RICATTI EQUATION
%
% *************************************************************************
% *************************************************************************

K_alg = icare(A,B,Q,R,[],[],[]);


%%
% *************************************************************************
% *************************************************************************
%
% OPTIMAL CONTROL, STATE, COST-TO-GO SEQUENCE
%
% *************************************************************************
% *************************************************************************

% *************************************************************************
%
% DATA STORAGE
% 

u_mat = zeros(1,N);         % Optimal control sequence.
x_mat = zeros(n,N);         % Optimal state trajectory.
% J_vec = zeros(1,N);         % Optimal cost-to-go at each stage.


% ***********************
%
% STATE, CONTROL
%

% Set initial condition.
x_mat(:,1) = x_0;

for k = 1:N-1
    
    x = x_mat(:,k);
    
    K = K_cell{k};
    
    u = - inv(R) * (B') * K * x;
    
    u_mat(:,k) = u;
    
    % Daynamics (approximated by ZOH transformed system).
    x_mat(:,k + 1) = A_dt * x + B_dt * u;
    
end



%%
% *************************************************************************
% *************************************************************************
%
% PLOT RESULTS
%
% *************************************************************************
% *************************************************************************


% ***********************
%
% STATE TRAJECTORY
%

% Title
ttl = 'Optimal State Trajectory x(t) = [x_0(t) x_1(t)]^T';

% Legend
lgd = [];
lgd{1} = 'x_0(t)';
lgd{2} = 'x_1(t)';


% Plot
figure(figcount)
h_fig = plot(tvec, x_mat(1,:), tvec, x_mat(2,:));
set(h_fig,'LineWidth',2)
hold on

% Formatting
grid on
title(ttl)
lgnd = legend(lgd);
set(lgnd, 'Location', 'Best')
xlabel('Time t (sec)')
ylabel('x(t)')


% SAVE PLOT
if savefigs
    filename = ['x_t_x0_eq_' num2str(x_0(1)) '_' num2str(x_0(2))];
    savepdf(figcount, relpath, filename); 
end

figcount = figcount + 1;        % Increment figure counter.



% ***********************
%
% CONTROL SEQUENCE
%

% Title
ttl = 'Optimal Control Sequence';

% Legend
lgd = [];

% Plot
figure(figcount)
h_fig = plot(tvec, u_mat);
set(h_fig,'LineWidth',2)
hold on

% Formatting
grid on
title(ttl)
% lgnd = legend(lgd);
set(lgnd, 'Location', 'Best')
xlabel('Time t (sec)')
ylabel('u(t)')


% SAVE PLOT
if savefigs
    filename = ['u_t_x0_eq_' num2str(x_0(1)) '_' num2str(x_0(2))];
    savepdf(figcount, relpath, filename); 
end

figcount = figcount + 1;        % Increment figure counter.


% ***********************
%
% K(t)
%

% Title
ttl = 'K(t)';

% Legend
lgd = [];
lgd{1} = 'K_{1,1}(t)';
lgd{2} = 'K_{1,1} (iCARE)';
lgd{3} = 'K_{1,2}(t)';
lgd{4} = 'K_{1,2} (iCARE)';
lgd{5} = 'K_{2,2}(t)';
lgd{6} = 'K_{2,2} (iCARE)';

% Plot
figure(figcount)

% K(1,1)
h_fig = plot(tvec, Kmat(:,1), tvec, ones(N,1)*K_alg(1,1));
set(h_fig,'LineWidth',2)

hold on

% K(1,2)
h_fig = plot(tvec, Kmat(:,2), tvec, ones(N,1)*K_alg(1,2));
set(h_fig,'LineWidth',2)

% K(2,2)
h_fig = plot(tvec, Kmat(:,4), tvec, ones(N,1)*K_alg(2,2));
set(h_fig,'LineWidth',2)

% Formatting
grid on
title(ttl)
lgnd = legend(lgd);
set(lgnd, 'Location', 'Best')
xlabel('Time t (sec)')
ylabel('K(t)')


% SAVE PLOT
if savefigs
    filename = ['K_t_x0_eq_' num2str(x_0(1)) '_' num2str(x_0(2))];
    savepdf(figcount, relpath, filename); 
end

figcount = figcount + 1;        % Increment figure counter.




%%
% *************************************************************************
% *************************************************************************
%
% MATRIX RICATTI DERIVATIVE FUNCTION
%
% To solve the differential Riccati equation with final value:
%
% 0 = K^{dot}(t) + Q - K(t) * B * R^{-1} * B^T * K(t) + K(t) * A + A^T *
% K(t),
%   t_0 <= t <= t_f,   K(t_f) = H.
%
% We solve the following differential Riccati equation with initial value:
%
% 0 = X^{dot}(t) - Q + X(t) * B * R^{-1} * B^T * X(t) - X(t) * A - A^T *
% X(t),
%   0 <= t <= t_f - t_0,   X(0) = H.
%
% Then 
%           K(t) = X( - (t - (t_f-t0)) ) 
%
% satisfies the differential Riccati equation with final value specified
% above.
% 
% *************************************************************************
% *************************************************************************



function dXdt = mric(t, X, A, B, Q, R)

% ***********************
%
% RESHAPE X FROM AN n^2 x 1 VECTOR TO AN n x n MATRIX
%
X = reshape(X, size(A)); 

% ***********************
%
% CALCULATE DERIVAVIVE
%
dXdt = Q - X * B * inv(R) * (B') * X + X * A + (A') * X;

% ***********************
%
% RESHAPE X FROM AN n x n MATRIX TO AN n^2 x 1 VECTOR
%
dXdt = dXdt(:); 

end


