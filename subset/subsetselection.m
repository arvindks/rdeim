function [p,err] = subsetselection(V, method)
    % Inputs
    % V       :  n x k  matrix with orthonormal columns
    % method  :  {'pqr', 'rrqr'} subset selection method
    % 
    % Outputs
    % p       :  k x 1 (subset of indices {1,..,n})
    % err     :  \| V(p,:)^+ \|_2  DEIM error constant
   
    k = size(V,2);
    if strcmp(method, 'pqr')
        [~,~,p] = qr(V',0);
        p = p(1:k); 
    elseif strcmp(method, 'rrqr')
        [p,~] = srrqr(V',k);    
    end
    
    err = norm(pinv(V(p,:)));
end
