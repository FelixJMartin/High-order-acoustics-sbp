``` matlab 


%% --- def parameters ---
m = 200;
h = 2/(m-1);
rho = 1;
c = 1;


%% --- construct matrices ---
c1 = rho * c^2 * speye(m^2);    % m^2×m^2
c2 = (1/rho) * speye(m^2);      % m^2×m^2
C_inv = blkdiag(c1, c2, c2);    % 3m^2×3m^2
%C = inv(C_inv);                 % 3m^2×3m^2


%   7th order upwind SBP operator         --------------


H=diag(ones(m,1),0);
H(1:6,1:6)=[0.19087e5 / 0.60480e5 0 0 0 0 0; 0 0.84199e5 / 0.60480e5 0 0 0 0; 0 0 0.18869e5 / 0.30240e5 0 0 0; 0 0 0 0.37621e5 / 0.30240e5 0 0; 0 0 0 0 0.55031e5 / 0.60480e5 0; 0 0 0 0 0 0.61343e5 / 0.60480e5;];
H(m-5:m,m-5:m)=fliplr(flipud(H(1:6,1:6)));
H=H*h;
HI=inv(H);

   
Qp=(-1/105*diag(ones(m-3,1),-3)+1/10*diag(ones(m-2,1),-2)-3/5*diag(ones(m-1,1),-1)-1/4*diag(ones(m,1),0)+1*diag(ones(m-1,1),+1)-3/10*diag(ones(m-2,1),+2)+1/15*diag(ones(m-3,1),+3)-1/140*diag(ones(m-4,1),+4));
Q_U =[-0.265e3 / 0.300272e6 0.1587945773e10 / 0.2432203200e10 -0.1926361e7 / 0.25737600e8 -0.84398989e8 / 0.810734400e9 0.48781961e8 / 0.4864406400e10 0.3429119e7 / 0.202683600e9; -0.1570125773e10 / 0.2432203200e10 -0.26517e5 / 0.1501360e7 0.240029831e9 / 0.486440640e9 0.202934303e9 / 0.972881280e9 0.118207e6 / 0.13512240e8 -0.231357719e9 / 0.4864406400e10; 0.1626361e7 / 0.25737600e8 -0.206937767e9 / 0.486440640e9 -0.61067e5 / 0.750680e6 0.49602727e8 / 0.81073440e8 -0.43783933e8 / 0.194576256e9 0.51815011e8 / 0.810734400e9; 0.91418989e8 / 0.810734400e9 -0.53314099e8 / 0.194576256e9 -0.33094279e8 / 0.81073440e8 -0.18269e5 / 0.107240e6 0.440626231e9 / 0.486440640e9 -0.365711063e9 / 0.1621468800e10; -0.62551961e8 / 0.4864406400e10 0.799e3 / 0.35280e5 0.82588241e8 / 0.972881280e9 -0.279245719e9 / 0.486440640e9 -0.346583e6 / 0.1501360e7 0.2312302333e10 / 0.2432203200e10; -0.3375119e7 / 0.202683600e9 0.202087559e9 / 0.4864406400e10 -0.11297731e8 / 0.810734400e9 0.61008503e8 / 0.1621468800e10 -0.1360092253e10 / 0.2432203200e10 -0.10677e5 / 0.42896e5;];

Qp(1:6,1:6)=Q_U;
Qp(m-5:m,m-5:m)=flipud( fliplr(Q_U(1:6,1:6) ) )'; %%% This is different from standard SBP

Qm=-Qp';

e_1=zeros(m,1);e_1(1)=1;
e_m=zeros(m,1);e_m(m)=1;

Dp=HI*(Qp-1/2*e_1*e_1'+1/2*e_m*e_m') ;

Dm=HI*(Qm-1/2*e_1*e_1'+1/2*e_m*e_m') ;

%%
Im = speye(m);          % mxm
dx1 = kron(Dp, Im);     % m^2xm^2
dx2 = kron(Dm, Im);     % m^2xm^2
z = sparse(m^2, m^2);   % m^2xm^2
Dx = [z dx1 z
      dx2 z z
      z  z  z];         % 3m^2×3m^2

%%
dy1 = kron(Im, Dp);     % m^2xm^2
dy2 = kron(Im, Dm);     % m^2xm^2
Dy = [z z dy1
      z  z  z
      dy2 z z];         % 3m^2×3m^2

%%
e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];

eW = kron(e_1, Im);
eE = kron(e_m, Im);
eS = kron(Im, e_1);
eN = kron(Im, e_m);

LW = kron(e2, eW.');   % v på W
LE = kron(e2, eE.');   % v på E
LS = kron(e3, eS.');   % w på S  
LN = kron(e3, eN.');   % w på N  

L = [LW
    LE
    LS
    LN];

%%
H_bar = kron(sparse(H), sparse(H));
H_tilde = kron(speye(3), H_bar)/C_inv;

%%
P = speye(3*m^2) - H_tilde\(L.'/(L*(H_tilde\L.'))*L);

M = L*(H_tilde\L.');

fprintf('size(L) = %d x %d\n', size(L,1), size(L,2));
fprintf('rank(L) = %d\n', rank(full(L)));
fprintf('size(M) = %d x %d\n', size(M,1), size(M,2));
fprintf('rank(M) = %d\n', rank(full(M)));
fprintf('condest(M) = %.3e\n', condest(M));

%% --- time-integrate with RK4
A = sparse(-P*C_inv*(Dx + Dy)*P);
CFL = 0.05;
x = linspace(-1, 1, m);     
y = linspace(-1, 1, m);
[X,Y] = meshgrid(x,y);

p0 = exp(-100*(X.^2 + Y.^2));   % m x m
v0 = zeros(m,m);               % m x m
w0 = zeros(m,m);               % m x m

u0 = [p0(:); v0(:); w0(:)];
u  = P*u0;

f = @(u) A*u; 

dt = CFL * h;
% t = 0 --> IV ----------------------------------
p_num = u(1:m^2);
v_num = u(m^2 + 1:2*m^2);
w_num = u(2*m^2 + 1:end);

Pnum = reshape(p_num, m, m);
Vnum = reshape(v_num, m, m);
Wnum = reshape(w_num, m, m);


% Visualize the results
figure;
surf(X, Y, Pnum);
xlabel('X-axis');
ylabel('Y-axis');
title('Pressure at t = 1')
zlabel('Pressure');
title('Pressure Distribution at t = 0');
colorbar;
% t = 1 --------------------------------------------
Nstep = round(1 / dt);

for n = 1:Nstep % stegar igenom varje tidssteg
    K1 = f(u);
    K2 = f(u + 0.5*dt*K1);
    K3 = f(u + 0.5*dt*K2);
    K4 = f(u + dt*K3);
    u = u + (dt/6)*(K1 + 2*K2 + 2*K3 + K4);
end

p_num = u(1:m^2);
v_num = u(m^2 + 1:2*m^2);
w_num = u(2*m^2 + 1:end);

Pnum = reshape(p_num, m, m);
Vnum = reshape(v_num, m, m);
Wnum = reshape(w_num, m, m);


% Visualize the results
figure;
surf(X, Y, Pnum);
xlabel('X-axis');
ylabel('Y-axis');
title('Pressure at t = 1')
zlabel('Pressure');
title('Pressure Distribution at t = 0.5');
colorbar;


% t = 2 ----------------------------------------------
Nstep = round(1 / dt);

for n = 1:Nstep % stegar igenom varje tidssteg
    K1 = f(u);
    K2 = f(u + 0.5*dt*K1);
    K3 = f(u + 0.5*dt*K2);
    K4 = f(u + dt*K3);
    u = u + (dt/6)*(K1 + 2*K2 + 2*K3 + K4);
end

p_num = u(1:m^2);
v_num = u(m^2 + 1:2*m^2);
w_num = u(2*m^2 + 1:end);

Pnum = reshape(p_num, m, m);
Vnum = reshape(v_num, m, m);
Wnum = reshape(w_num, m, m);


% Visualize the results
figure;
surf(X, Y, Pnum);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Pressure');
title('Pressure Distribution at t = 2');
colorbar;


%% TASK 2 ------------------------------------------------------------ 
% 1) Parameters + grid ---
m = 200;
h = 2/(m-1);

x = linspace(-1, 1, m);
y = linspace(-1, 1, m);
[X,Y] = meshgrid(x,y);

CFL  = 0.05;

rho  = ones(m,m);
c    = ones(m,m);
beta = zeros(m,m);

right   = (X > 0);      % right half-plane
rho_new = 0.5;
rho(right) = rho_new;   % material change

% --- 2) Material matrix C_inv (depends on rho,c) ---
rho_vec = rho(:);
c_vec   = c(:);

c1 = spdiags(rho_vec .* (c_vec.^2), 0, m^2, m^2);   % rho*c^2
c2 = spdiags(1./rho_vec,            0, m^2, m^2);   % 1/rho

C_inv = blkdiag(c1, c2, c2);                         % 3m^2 x 3m^2

%  7th order upwind SBP operator   


H=diag(ones(m,1),0);
H(1:6,1:6)=[0.19087e5 / 0.60480e5 0 0 0 0 0; 0 0.84199e5 / 0.60480e5 0 0 0 0; 0 0 0.18869e5 / 0.30240e5 0 0 0; 0 0 0 0.37621e5 / 0.30240e5 0 0; 0 0 0 0 0.55031e5 / 0.60480e5 0; 0 0 0 0 0 0.61343e5 / 0.60480e5;];
H(m-5:m,m-5:m)=fliplr(flipud(H(1:6,1:6)));
H=H*h;
HI=inv(H);

   
Qp=(-1/105*diag(ones(m-3,1),-3)+1/10*diag(ones(m-2,1),-2)-3/5*diag(ones(m-1,1),-1)-1/4*diag(ones(m,1),0)+1*diag(ones(m-1,1),+1)-3/10*diag(ones(m-2,1),+2)+1/15*diag(ones(m-3,1),+3)-1/140*diag(ones(m-4,1),+4));
Q_U =[-0.265e3 / 0.300272e6 0.1587945773e10 / 0.2432203200e10 -0.1926361e7 / 0.25737600e8 -0.84398989e8 / 0.810734400e9 0.48781961e8 / 0.4864406400e10 0.3429119e7 / 0.202683600e9; -0.1570125773e10 / 0.2432203200e10 -0.26517e5 / 0.1501360e7 0.240029831e9 / 0.486440640e9 0.202934303e9 / 0.972881280e9 0.118207e6 / 0.13512240e8 -0.231357719e9 / 0.4864406400e10; 0.1626361e7 / 0.25737600e8 -0.206937767e9 / 0.486440640e9 -0.61067e5 / 0.750680e6 0.49602727e8 / 0.81073440e8 -0.43783933e8 / 0.194576256e9 0.51815011e8 / 0.810734400e9; 0.91418989e8 / 0.810734400e9 -0.53314099e8 / 0.194576256e9 -0.33094279e8 / 0.81073440e8 -0.18269e5 / 0.107240e6 0.440626231e9 / 0.486440640e9 -0.365711063e9 / 0.1621468800e10; -0.62551961e8 / 0.4864406400e10 0.799e3 / 0.35280e5 0.82588241e8 / 0.972881280e9 -0.279245719e9 / 0.486440640e9 -0.346583e6 / 0.1501360e7 0.2312302333e10 / 0.2432203200e10; -0.3375119e7 / 0.202683600e9 0.202087559e9 / 0.4864406400e10 -0.11297731e8 / 0.810734400e9 0.61008503e8 / 0.1621468800e10 -0.1360092253e10 / 0.2432203200e10 -0.10677e5 / 0.42896e5;];

Qp(1:6,1:6)=Q_U;
Qp(m-5:m,m-5:m)=flipud( fliplr(Q_U(1:6,1:6) ) )'; %%% This is different from standard SBP

Qm=-Qp';

e_1=zeros(m,1);e_1(1)=1;
e_m=zeros(m,1);e_m(m)=1;

Dp=HI*(Qp-1/2*e_1*e_1'+1/2*e_m*e_m') ;

Dm=HI*(Qm-1/2*e_1*e_1'+1/2*e_m*e_m') ;


Im = speye(m);          % mxm
dx1 = kron(Dp, Im);     % m^2xm^2
dx2 = kron(Dm, Im);     % m^2xm^2
z = sparse(m^2, m^2);   % m^2xm^2
Dx = [z dx1 z
      dx2 z z
      z  z  z];         % 3m^2×3m^2


dy1 = kron(Im, Dp);     % m^2xm^2
dy2 = kron(Im, Dm);     % m^2xm^2
Dy = [z z dy1
      z  z  z
      dy2 z z];         % 3m^2×3m^2


e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];

eW = kron(e_1, Im);
eE = kron(e_m, Im);
eS = kron(Im, e_1);
eN = kron(Im, e_m);

LW = kron(e2, eW.');   % v på W
LE = kron(e2, eE.');   % v på E
LS = kron(e3, eS.');   % w på S   <-- ÄNDRING
LN = kron(e3, eN.');   % w på N   <-- ÄNDRING

L = [LW
    LE
    LS
    LN];


H_bar = kron(sparse(H), sparse(H));
H_tilde = kron(speye(3), H_bar)/C_inv;


P = speye(3*m^2) - H_tilde\(L.'/(L*(H_tilde\L.'))*L);

M = L*(H_tilde\L.');

fprintf('size(L) = %d x %d\n', size(L,1), size(L,2));
fprintf('rank(L) = %d\n', rank(full(L)));
fprintf('size(M) = %d x %d\n', size(M,1), size(M,2));
fprintf('rank(M) = %d\n', rank(full(M)));
fprintf('condest(M) = %.3e\n', condest(M));

% --- time-integrate with RK4
beta_vec = beta(:);   % m^2 x 1
D = blkdiag(spdiags(beta_vec,0,m^2,m^2), sparse(m^2,m^2), sparse(m^2,m^2));
A = sparse(-P*C_inv*(Dx + Dy + D)*P);
CFL = 0.05;

p0 = exp(-100*(X.^2 + Y.^2));   % m x m
v0 = zeros(m,m);               % m x m
w0 = zeros(m,m);               % m x m

u0 = [p0(:); v0(:); w0(:)];
u  = P*u0;

f = @(u) A*u;

dt = CFL * h;
% t = 0 --> IV ----------------------------------
p_num = u(1:m^2);
v_num = u(m^2 + 1:2*m^2);
w_num = u(2*m^2 + 1:end);

Pnum = reshape(p_num, m, m);
Vnum = reshape(v_num, m, m);
Wnum = reshape(w_num, m, m);


% Visualize the results
figure;
surf(X, Y, Pnum);
xlabel('X-axis');
ylabel('Y-axis');
title('Pressure at t = 1')
zlabel('Pressure');
title('Pressure Distribution at t = 0');
colorbar;
% t = 1 --------------------------------------------
Nstep = round(1 / dt);

for n = 1:Nstep % stegar igenom varje tidssteg
    K1 = f(u);
    K2 = f(u + 0.5*dt*K1);
    K3 = f(u + 0.5*dt*K2);
    K4 = f(u + dt*K3);
    u = u + (dt/6)*(K1 + 2*K2 + 2*K3 + K4);
end

p_num = u(1:m^2);
v_num = u(m^2 + 1:2*m^2);
w_num = u(2*m^2 + 1:end);

Pnum = reshape(p_num, m, m);
Vnum = reshape(v_num, m, m);
Wnum = reshape(w_num, m, m);


% Visualize the results
figure;
surf(X, Y, Pnum);
xlabel('X-axis');
ylabel('Y-axis');
title('Pressure at t = 1')
zlabel('Pressure');
title('Pressure Distribution at t = 1');
colorbar;


% t = 2 ----------------------------------------------
Nstep = round(1 / dt);

for n = 1:Nstep % stegar igenom varje tidssteg
    K1 = f(u);
    K2 = f(u + 0.5*dt*K1);
    K3 = f(u + 0.5*dt*K2);
    K4 = f(u + dt*K3);
    u = u + (dt/6)*(K1 + 2*K2 + 2*K3 + K4);
end

p_num = u(1:m^2);
v_num = u(m^2 + 1:2*m^2);
w_num = u(2*m^2 + 1:end);

Pnum = reshape(p_num, m, m);
Vnum = reshape(v_num, m, m);
Wnum = reshape(w_num, m, m);


% Visualize the results
figure;
surf(X, Y, Vnum);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Velocity'); 
title(sprintf('Velocity distribution at t = 2, rho_{new} = %.4f', rho_new))
colorbar;

```
