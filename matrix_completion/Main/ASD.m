function [Mout, Out] = ASD(m,n,r,Omega,data,start,opts)

reltol = opts.rel_res_tol;
maxiter = opts.maxit;
verbosity = opts.verbosity;
rate_limit = 1-opts.rel_res_change_tol;
relres = reltol * norm(data);

p = length(data);

% create a sparse matrix
[I, J] = ind2sub([m,n],Omega);
diff_on_omega_matrix = sparse(I,J,data,m,n,p);

if ~isempty(start)
  X = start.L;
  Y = start.R;
  if verbosity
    fprintf('initial point provided!\n')
  end
else
  X = randn(m,r);
  Y = randn(r,n);
  if verbosity
    fprintf('no initial point!\n')
  end
end

clear opts;
clear start;

Xt = X';
diff_on_omega = data'-partXY(Xt,Y,I,J,p);
res = norm(diff_on_omega);

iter = 1;

itres = zeros(maxiter,1);
itres(iter) = res;

conv_rate = 0;

% iteration
while iter <= maxiter &&  res >=  relres && conv_rate <= rate_limit
    % gradient for X
    updateSval(diff_on_omega_matrix,diff_on_omega,p);
    grad_X = diff_on_omega_matrix*Y';
       
    % stepsize for X
    grad_Xt = grad_X';
    delta_XY = partXY(grad_Xt,Y,I,J,p);
    tx = norm(grad_X,'fro')^2/norm(delta_XY)^2;
    
    % update X
    X = X + tx*grad_X;
    
    diff_on_omega = diff_on_omega-tx*delta_XY;
    
   % gradient for Y
    updateSval(diff_on_omega_matrix,diff_on_omega,p);
    Xt = X';
    grad_Y = Xt*diff_on_omega_matrix;
    
    % stepsize for Y
    delta_XY = partXY(Xt,grad_Y,I,J,p);
    ty = norm(grad_Y,'fro')^2/norm(delta_XY)^2;
    
    % update Y
    Y = Y + ty*grad_Y;
   
    diff_on_omega = diff_on_omega-ty*delta_XY;
    res = norm(diff_on_omega);
    
    iter = iter + 1;
    
    itres(iter) = res;
    
    conv_rate = (itres(iter)/itres(max(1,iter-15)))^(1/min(15,iter-1));
end

Mout = X*Y;

Out.itrelres = itres(1:iter)/norm(data);
Out.iter = iter;
Out.reschg = abs(1-conv_rate);
