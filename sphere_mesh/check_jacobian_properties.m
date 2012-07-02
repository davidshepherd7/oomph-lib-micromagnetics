
% Get jacobian
load matrices/cr_jacobian;
A = ascii2matlab(cr_jacobian);

%Get iLU preconditioner
[L,U] = ilu(A,struct('type','ilutp','droptol',1e-5));

%Set B to be something...
load matrices/residual;
B = doublevec2matlab(residual);

% gmres solve
tol = 1e-6;
max_it = size(A,1); % max allowed by matlab

% Try solving without preconditioner:
display('Without preconditioner:')
[x_np,flag_np,relres_np,iter_np] = gmres(A,B,[],tol,max_it);
flag_np
relres_np
iter_np



% and with an iLU preconditioner
display('With precondtioner:')
[x,flag,relres,iter] = gmres(A,B,[],tol,max_it,L,U);
flag
relres
iter