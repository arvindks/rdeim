function [u,s,v] = randsvd(A,Q,k)
  
    B = Q'*A;
    [ub,s,v] = svd(B,0);
    u = Q*ub;       u = u(:,1:k);
    s = diag(s);    s = s(1:k);
    v = v(:,1:k);
end