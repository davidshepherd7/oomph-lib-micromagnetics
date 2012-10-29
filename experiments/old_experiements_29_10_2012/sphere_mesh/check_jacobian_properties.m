
% Get jacobian
load matrices/cr_jacobian;
A = ascii2matlab(cr_jacobian);

%Get condition number - not accurate...
%condition_number = cond(full(A))

%Set B to be something...
load matrices/residual;
B = doublevec2matlab(residual);

%Get iLU preconditioner
[L,U] = ilu(A,struct('type','ilutp','droptol',1e-5));

% gmres options
tol = 1e-10;
max_it = size(A,1); % max allowed by matlab

% Try solving without preconditioner:
display('Without preconditioner:')
[x_np,flag_np,relres_np,iter_np] = gmres(A,B,[],tol,max_it);
flag_np
relres_np
iter_np
max_res_np = max(x_np)

% and with an iLU preconditioner
display('With precondtioner:')
[x_pre,flag_pre,relres_pre,iter_pre] = gmres(A,B,[],tol,max_it,L,U);
flag_pre
relres_pre
iter_pre
max_res_pre = max(x_pre)

% just LU solve:
display('Pure LU solve:')
[L_complete,U_complete,P] = ilu(A,struct('type','ilutp','droptol',0)); %P should be identity
y_lu = mldivide(L_complete,P*B); %solving system using LU decomp: we first invert wrt L
x_lu = mldivide(U_complete,y_lu);  %then wrt U, both steps are easy because triangular.
max_res_lu = max(x_lu)

