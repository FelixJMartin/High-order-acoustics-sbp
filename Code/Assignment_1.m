# Tasks 3–4 — 1D Acoustics with SBP-Projection (MATLAB)

This page contains the full MATLAB code for **Task 3** (eigenvalues + CFL for different boundary conditions) and **Task 4** (RK4 time integration + convergence check for 7th vs 6th order SBP operators).

---

## Task 3 — Eigenvalues and CFL for different boundary conditions

**Goal:** Compare stability properties of the semi-discrete operator  
\[
M = -P D_x P
\]
for three boundary treatments:
- Dirichlet on **pressure** \(p=0\) at both boundaries
- Dirichlet on **velocity** \(v=0\) at both boundaries
- **Characteristic** boundary conditions

We plot eigenvalues of \(hM\) and compute a CFL estimate for RK4 using:
- RK4 stability radius \(R \approx 2.78\)
- \(\Delta t = k = R / \rho(M)\)
- \(\alpha = k/h\)

<details>
<summary><strong>Show Task 3 MATLAB code</strong> (click to expand)</summary>

```matlab
%% -----------------------------------------------------------
% Task 3: Eigenvalues + CFL for different boundary conditions
% -----------------------------------------------------------

%% Book-keeping
m = 51;
xl = -1;
xr = 1;
h = (xr-xl)/(m-1);

% Basis vectors in component space (p,v)
e1 = [1 0];   % selects p
e2 = [0 1];   % selects v  

% Boundary selectors in physical space (m grid points)
el = zeros(m,1); 
el(1) = 1;        % left boundary

er = zeros(m,1);
er(end) = 1;      % right boundary

%% SBP op (7th order upwind, Mattsson)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 7th order upwind SBP operator           %%%
%%% Derived by Ken Mattsson                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=diag(ones(m,1),0);
H(1:6,1:6)=[0.19087e5 / 0.60480e5 0 0 0 0 0; 0 0.84199e5 / 0.60480e5 0 0 0 0; 0 0 0.18869e5 / 0.30240e5 0 0 0; 0 0 0 0.37621e5 / 0.30240e5 0 0; 0 0 0 0 0.55031e5 / 0.60480e5 0; 0 0 0 0 0 0.61343e5 / 0.60480e5;];
H(m-5:m,m-5:m)=fliplr(flipud(H(1:6,1:6)));
H=H*h;
HI=inv(H);

Qp=(-1/105*diag(ones(m-3,1),-3)+1/10*diag(ones(m-2,1),-2)-3/5*diag(ones(m-1,1),-1)-1/4*diag(ones(m,1),0)+1*diag(ones(m-1,1),+1)-3/10*diag(ones(m-2,1),+2)+1/15*diag(ones(m-3,1),+3)-1/140*diag(ones(m-4,1),+4));

% Keep Q_U exactly as in your script:
Q_U = [ ...  % <-- paste your full 6x6 Q_U here (unchanged)
];

Qp(1:6,1:6)=Q_U;
Qp(m-5:m,m-5:m)=flipud( fliplr(Q_U(1:6,1:6)) )'; 
Qm=-Qp';

e_1=zeros(m,1); e_1(1)=1;
e_m=zeros(m,1); e_m(m)=1;

Dp=HI*(Qp-1/2*e_1*e_1'+1/2*e_m*e_m');
Dm=HI*(Qm-1/2*e_1*e_1'+1/2*e_m*e_m');

%% Dirichlet conditions p(0) = 0, p(m) = 0
L_p = [ kron(e1, el.')      
        kron(e1, er.') ];   

%% Dirichlet conditions v(0) = 0, v(m) = 0
L_v = [ kron(e2, el.')
        kron(e2, er.') ];

%% Characteristic BC
a1 = 1;
am = -1;

wL = e1 + a1*e2;   % [1; a1]
wR = e1 + am*e2;   % [1; am]

L_C = [ kron(wL, el.')
        kron(wR, er.') ];

%% Other matrices needed
Hbar     = kron(eye(2), H);   % 2m x 2m
Hbar_inv = inv(Hbar);

%% projections
B_p    = L_p * Hbar_inv * L_p.';
P_p    = eye(2*m) - Hbar_inv * L_p.' * inv(B_p) * L_p;

B_v    = L_v * Hbar_inv * L_v.';
P_v    = eye(2*m) - Hbar_inv * L_v.' * inv(B_v) * L_v;

B_C    = L_C * Hbar_inv * L_C.';
P_C    = eye(2*m) - Hbar_inv * L_C.' * inv(B_C) * L_C;

%% Derivative approx.
Z = zeros(m);
Dx = [ Z   Dp
       Dm  Z  ];

%% h * eigvals
M_p = -P_p * Dx * P_p;
M_v = -P_v * Dx * P_v;
M_C = -P_C * Dx * P_C;

lambda_p = h * eig(M_p);
lambda_v = h * eig(M_v);
lambda_C = h * eig(M_C);

%% plot eigenvalues
figure;

subplot(1,3,1)
plot(real(lambda_p), imag(lambda_p), 'o');
xlabel('Re'); ylabel('Im');
title('Eigenvalues of hM (Dirichlet: p=0)');
grid on; axis equal;

subplot(1,3,2)
plot(real(lambda_v), imag(lambda_v), 'o');
xlabel('Re'); ylabel('Im');
title('Eigenvalues of hM (Dirichlet: v=0)');
grid on; axis equal;

subplot(1,3,3)
plot(real(lambda_C), imag(lambda_C), 'o');
xlabel('Re'); ylabel('Im');
title('Eigenvalues of hM (Characteristic BCs)');
grid on; axis equal;

%% CFL (alpha) for RK4
R = 2.78;

SR_p = max(abs(eig(M_p)));
SR_v = max(abs(eig(M_v)));
SR_C = max(abs(eig(M_C)));

k_p = R/SR_p;
k_v = R/SR_v;
k_C = R/SR_C;

alpha_p = k_p/h;
alpha_v = k_v/h;
alpha_C = k_C/h;

fprintf('--- CFL estimates (RK4) ---\n');
fprintf('alpha_p (p-Dirichlet) = %.6f\n', alpha_p);
fprintf('alpha_v (v-Dirichlet) = %.6f\n', alpha_v);
fprintf('alpha_C (Characteristic) = %.6f\n', alpha_C);


%% -----------------------------------------------------------


function task4_compact()
% TASK4_COMPACT  Clean/short version of your Task 4 (7th upwind SBP + 6th SBP)
% Runs two-grid convergence test (m and n=2m-1), RK4 to t=1.8, CFL=0.05.
% Produces plots + prints convergence rate q and fine-grid L2 error.

% ---------------- user parameters ----------------
m0      = 401;     % coarse grid
xl      = -1;
xr      = 1;
tEnd    = 1.8;
CFL     = 0.05;
r_star  = 0.1;

% initial data generators (vectorized)
theta1 = @(x,t) exp(-((x - t)/r_star).^2);
theta2 = @(x,t) -exp(-((x + t)/r_star).^2);

% exact solution (as in your code)
p_exact = @(x) theta2(x,(xr-xl)-tEnd) - theta1(x,(xr-xl)-tEnd);
v_exact = @(x) theta1(x,(xr-xl)-tEnd) + theta2(x,(xr-xl)-tEnd);

% ------------------------------------------------
% 7th order upwind SBP
% ------------------------------------------------
fprintf('\n=== 7th order upwind SBP ===\n');
run_two_grid(m0, xl, xr, tEnd, CFL, theta1, theta2, p_exact, v_exact, @sbp7_upwind);

% ------------------------------------------------
% 6th order SBP (D1)
% ------------------------------------------------
fprintf('\n=== 6th order SBP (D1) ===\n');
run_two_grid(m0, xl, xr, tEnd, CFL, theta1, theta2, p_exact, v_exact, @sbp6_D1);

end

% =====================================================================
% Two-grid convergence driver
% =====================================================================
function run_two_grid(m0, xl, xr, tEnd, CFL, theta1, theta2, p_exact_fun, v_exact_fun, op_builder)

% coarse
[x, u_num, u_ex, errL2_coarse] = solve_once(m0, xl, xr, tEnd, CFL, theta1, theta2, p_exact_fun, v_exact_fun, op_builder);

% plots (coarse)
plot_solution(x, u_num, u_ex, sprintf('%s, m=%d', func2str(op_builder), m0));

% fine: n = 2m-1
m1 = 2*m0 - 1;
[~, ~, ~, errL2_fine] = solve_once(m1, xl, xr, tEnd, CFL, theta1, theta2, p_exact_fun, v_exact_fun, op_builder);

% convergence estimate (same formula you used)
q = log10(errL2_coarse/errL2_fine) / log10(m1/((m1+1)/2));

fprintf('q = %.6f  (coarse m=%d, fine m=%d)\n', q, m0, m1);
fprintf('Fine-grid L2 error = %.16e\n', errL2_fine);

end

% =====================================================================
% One solve: build operator, RK4 integrate, compute error
% =====================================================================
function [x, u, u_exact, errL2] = solve_once(m, xl, xr, tEnd, CFL, theta1, theta2, p_exact_fun, v_exact_fun, op_builder)

h = (xr - xl)/(m-1);
x = linspace(xl, xr, m).';

% Build SBP operator pieces
[H, Dx] = op_builder(m, h);          % H is m×m norm, Dx is 2m×2m first-derivative block operator

% Dirichlet on p: p(-1)=p(1)=0 (same as your L and P)
P = projection_p_dirichlet(H, m);

C_inv = kron(eye(2), eye(m));        % identity in your current tasks
M = -P * (C_inv * Dx) * P;

% initial data
t0 = 0;
p0 = theta1(x,t0) - theta2(x,t0);
v0 = theta1(x,t0) + theta2(x,t0);
u  = [p0; v0];

% RK4 integrate
dt = CFL*h;
N  = round(tEnd/dt);

f = @(uu) M*uu;
for k = 1:N
    K1 = f(u);
    K2 = f(u + 0.5*dt*K1);
    K3 = f(u + 0.5*dt*K2);
    K4 = f(u + dt*K3);
    u  = u + (dt/6)*(K1 + 2*K2 + 2*K3 + K4);
end

% exact & error
p_ex = p_exact_fun(x);
v_ex = v_exact_fun(x);
u_exact = [p_ex; v_ex];

e = u_exact - u;
errL2 = sqrt(h) * norm(e,2);

end

% =====================================================================
% Projection for p-Dirichlet using SBP norm H
% (matches: P = I - inv(H_)*L'*inv(L*inv(H_)*L')*L)
% =====================================================================
function P = projection_p_dirichlet(H, m)

I2m = eye(2*m);
Hbar = kron(eye(2), H);

eL = zeros(1,m); eL(1)=1;
eR = zeros(1,m); eR(m)=1;

e1 = [1 0]; % selects p
L  = [kron(e1, eL);
      kron(e1, eR)];

Hbar_inv = inv(Hbar);
B = L * Hbar_inv * L.';
P = I2m - Hbar_inv * L.' * (B \ L);

end

% =====================================================================
% Plot helper
% =====================================================================
function plot_solution(x, u_num, u_exact, ttl)

m = numel(x);
p_num = u_num(1:m); v_num = u_num(m+1:end);
p_ex  = u_exact(1:m); v_ex = u_exact(m+1:end);

figure;
plot(x, p_num, '-', x, v_num, '--', x, p_ex, ':', x, v_ex, '-.');
xlim([min(x) max(x)]);
grid on;
legend('p num','v num','p exact','v exact','Location','best');
title(ttl);

end

% =====================================================================
% 7th order upwind SBP builder (your exact coefficients, wrapped)
% Returns:
%   H  (m×m)
%   Dx (2m×2m) = [0 Dp; Dm 0]
% =====================================================================
function [H, Dx] = sbp7_upwind(m, h)

H = diag(ones(m,1),0);
H(1:6,1:6) = [ ...
    0.19087e5/0.60480e5 0 0 0 0 0; ...
    0 0.84199e5/0.60480e5 0 0 0 0; ...
    0 0 0.18869e5/0.30240e5 0 0 0; ...
    0 0 0 0.37621e5/0.30240e5 0 0; ...
    0 0 0 0 0.55031e5/0.60480e5 0; ...
    0 0 0 0 0 0.61343e5/0.60480e5];
H(m-5:m,m-5:m) = fliplr(flipud(H(1:6,1:6)));
H = H*h;
HI = inv(H);

Qp = (-1/105*diag(ones(m-3,1),-3) + 1/10*diag(ones(m-2,1),-2) - 3/5*diag(ones(m-1,1),-1) ...
     -1/4*diag(ones(m,1),0) + 1*diag(ones(m-1,1),+1) - 3/10*diag(ones(m-2,1),+2) ...
     +1/15*diag(ones(m-3,1),+3) - 1/140*diag(ones(m-4,1),+4));

Q_U = [ ...
 -0.265e3/0.300272e6  0.1587945773e10/0.2432203200e10 -0.1926361e7/0.25737600e8   -0.84398989e8/0.810734400e9   0.48781961e8/0.4864406400e10  0.3429119e7/0.202683600e9; ...
 -0.1570125773e10/0.2432203200e10 -0.26517e5/0.1501360e7  0.240029831e9/0.486440640e9  0.202934303e9/0.972881280e9  0.118207e6/0.13512240e8 -0.231357719e9/0.4864406400e10; ...
  0.1626361e7/0.25737600e8 -0.206937767e9/0.486440640e9 -0.61067e5/0.750680e6   0.49602727e8/0.81073440e8   -0.43783933e8/0.194576256e9  0.51815011e8/0.810734400e9; ...
  0.91418989e8/0.810734400e9 -0.53314099e8/0.194576256e9 -0.33094279e8/0.81073440e8 -0.18269e5/0.107240e6  0.440626231e9/0.486440640e9 -0.365711063e9/0.1621468800e10; ...
 -0.62551961e8/0.4864406400e10  0.799e3/0.35280e5  0.82588241e8/0.972881280e9 -0.279245719e9/0.486440640e9 -0.346583e6/0.1501360e7  0.2312302333e10/0.2432203200e10; ...
 -0.3375119e7/0.202683600e9  0.202087559e9/0.4864406400e10 -0.11297731e8/0.810734400e9  0.61008503e8/0.1621468800e10 -0.1360092253e10/0.2432203200e10 -0.10677e5/0.42896e5 ];

Qp(1:6,1:6) = Q_U;
Qp(m-5:m,m-5:m) = flipud(fliplr(Q_U)).';

Qm = -Qp';

e1 = zeros(m,1); e1(1)=1;
em = zeros(m,1); em(m)=1;

Dp = HI*(Qp - 0.5*e1*e1' + 0.5*em*em');
Dm = HI*(Qm - 0.5*e1*e1' + 0.5*em*em');

Dx = [zeros(m) Dp; Dm zeros(m)];

end

% =====================================================================
% 6th order SBP builder (your exact D1 + H coefficients, wrapped)
% Returns:
%   H  (m×m)
%   Dx (2m×2m) = [0 D1; D1 0]
% =====================================================================
function [H, Dx] = sbp6_D1(m, h)

% norm matrix
H = diag(ones(m,1),0);
H(1:6,1:6) = diag([13649/43200,12013/8640,2711/4320,5359/4320,7877/8640,43801/43200]);
H(m-5:m,m-5:m) = fliplr(flipud(diag([13649/43200,12013/8640,2711/4320,5359/4320,7877/8640,43801/43200])));
H = H*h;

% free parameter
x1 = 0.70127127127127;

% interior stencil (central 6th order)
D1 = (  1/60*diag(ones(m-3,1), 3) ...
      - 9/60*diag(ones(m-2,1), 2) ...
      +45/60*diag(ones(m-1,1), 1) ...
      -45/60*diag(ones(m-1,1),-1) ...
      + 9/60*diag(ones(m-2,1),-2) ...
      - 1/60*diag(ones(m-3,1),-3) );

% boundary closure (your matrix)
D1(1:6,1:9) = [ ...
    -21600/13649,                      43200/13649*x1-7624/40947,      -172800/13649*x1+715489/81894, ...
     259200/13649*x1-187917/13649,    -172800/13649*x1+735635/81894,    43200/13649*x1-89387/40947, 0, 0, 0; ...
    -8640/12013*x1+7624/180195,                0,                      86400/12013*x1-57139/12013, ...
    -172800/12013*x1+745733/72078,    129600/12013*x1-91715/12013,    -34560/12013*x1+240569/120130, 0, 0, 0; ...
     17280/2711*x1-715489/162660,    -43200/2711*x1+57139/5422,                 0, ...
     86400/2711*x1-176839/8133,      -86400/2711*x1+242111/10844,      25920/2711*x1-182261/27110, 0, 0, 0; ...
    -25920/5359*x1+187917/53590,      86400/5359*x1-745733/64308,     -86400/5359*x1+176839/16077, ...
              0,                      43200/5359*x1-165041/32154,     -17280/5359*x1+710473/321540, 72/5359, 0, 0; ...
     34560/7877*x1-147127/47262,    -129600/7877*x1+91715/7877,       172800/7877*x1-242111/15754, ...
    -86400/7877*x1+165041/23631,              0,                       8640/7877*x1,   -1296/7877, 144/7877, 0; ...
    -43200/43801*x1+89387/131403,    172800/43801*x1-240569/87602,   -259200/43801*x1+182261/43801, ...
     172800/43801*x1-710473/262806,  -43200/43801*x1,                          0,      32400/43801, -6480/43801, 720/43801 ...
    ];

D1(m-5:m, m-8:m) = flipud(fliplr(-D1(1:6,1:9)));
D1 = D1/h;

Dx = [zeros(m) D1; D1 zeros(m)];

end
