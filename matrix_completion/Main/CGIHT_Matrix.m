function [Mout, Out] = CGIHT_Matrix(m,n,r,Omega,data,start,opts)

reltol = opts.rel_res_tol; 
maxiter = opts.maxit;
rate_limit = 1-opts.rel_res_change_tol;
relres = reltol*norm(data);
itres = zeros(maxiter,1);
conv_rate = 0;

U = start.U;
S = diag(start.sigma);
V = start.V;

clear opts;
clear start;

p = length(data);
[I,J] = ind2sub([m,n],Omega);
res_on_omega_matrix = sparse(I,J,data,m,n,p);

M = U*S*V';
res_on_omega = data - M(Omega); % p*1 vector
updateSval(res_on_omega_matrix,res_on_omega,p);
res = norm(res_on_omega);


iter = 1;
itres(iter) = res;

% note that sparse(full) is necessary
% or else d will point to the same area as res_on_omega_matrix
d = sparse(full(res_on_omega_matrix));
Ad_proj = partXY(U',U'*d, I, J, p);

% iteration
while  ((res >= relres) && (iter <= maxiter)) % && (conv_rate < rate_limit)) 
    % compute alpha
    % <UU'*r, UU'*d> = <r, UU'*d> = <r, A.*(UU'*d)>
    % since res_on_omega is sparse
    % alpha = d_proj*res_on_omega/norm(d_proj)^2;
    alpha = Ad_proj*res_on_omega;
    tmp = norm(Ad_proj)^2;
    if abs(alpha) < 1000*tmp
      alpha = alpha/tmp;
    else
      alpha = 1;
    end
    
    % compute M+alpha*p_new
    M = M + alpha * d;
    
    [U,S,V] = svds(M,r);
    M = U*S*V'; % can be further accelerated here

    res_on_omega = data - M(Omega); % p*1 vector
    updateSval(res_on_omega_matrix,res_on_omega,p);    
    res = norm(res_on_omega);

    % res_proj = A.*(UU'*res_on_omega) = A.*(U(U'*res_on_omega))
    % d_proj = A.*(UU'*d) = A.*(U(U'*d))
    Ares_proj = partXY(U', U'*res_on_omega_matrix, I, J, p);
    Ad_proj = partXY(U', U'*d, I, J, p);
    
    %beta = res_proj*d_proj'/norm(d_proj)^2; 
    beta = Ares_proj*Ad_proj';
    tmp = norm(Ad_proj)^2;

    if abs(beta) < 1000 * tmp
      beta = beta/tmp;
    else
      beta = 0;
    end
   
    d = res_on_omega_matrix - beta*d; 

    % update d_proj
    Ad_proj = Ares_proj - beta * Ad_proj;

    iter = iter + 1;
    
    itres(iter) = res;
    conv_rate = (itres(iter)/itres(max(1,iter-15)))^(1/min(15,iter-1));
    
    if iter >= 500 && conv_rate >= rate_limit
        break;
    end
end

Mout = M;

Out.itrelres = itres(1:iter)/norm(data);
Out.iter = iter;
Out.reschg = abs(1-conv_rate);

