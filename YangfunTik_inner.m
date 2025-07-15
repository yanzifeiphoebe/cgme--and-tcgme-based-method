
function [x_04,x_05,x_06,err_04,err_05,err_06] = YangfunTik_inner(A,xtrue,X,V,L,k)

 xnorm = norm(L*xtrue);
 N=size(A,2);
% [x_k,err_k,res_k,sol_k,flag,relres,iter] = YangfunTik(A,b,xtrue,X,V,L,k)
%%%%%% Compute the lsqrTik and the relative errors %%%%%%%%%%%%%%%%%%
for i=1:k
    Q = V(:,1:i);
    d=L*X(:,i);
    z_04 = lsqr(@(z,tflag)afun(z,L,Q,tflag),d,1e-4,N);
    z_05 = lsqr(@(z,tflag)afun(z,L,Q,tflag),d,1e-5,N);
    z_06 = lsqr(@(z,tflag)afun(z,L,Q,tflag),d,1e-6,N);
    x_04(:,i) = X(:,i) - z_04;
    x_05(:,i) = X(:,i) - z_05;
    x_06(:,i) = X(:,i) - z_06;
    err_04(i) = norm(L*(x_04(:,i)-xtrue))/xnorm;
    err_05(i) = norm(L*(x_05(:,i)-xtrue))/xnorm;
    err_06(i) = norm(L*(x_06(:,i)-xtrue))/xnorm;
%     err_freeL(i) = norm((x_k(:,i)-xtrue))/norm(xtrue);
%     rho(i,1) = norm(A*x_k(:,i)-b);
%     eta(i,1) = norm(L*x_k(:,i));
    
    
end
end
function y =afun(z,A,B,transp_flag)
if strcmp(transp_flag,'transp')   % y = (A(I_n-BB^T))' * z;
    s = A'*z;
    t=B'*s;
    t=B*t;
    y = s-t;
elseif strcmp(transp_flag,'notransp') % y = (A(I_n-BB^T)) * z;
    s = A*z;
    t=B'*z;
    t=B*t;
    t = A*t;
    y = s-t;
    
end
end