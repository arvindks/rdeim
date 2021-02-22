function p = deim(V)
%
% This function determines the indices for a "good" set of n rows
% from an input matrix V using the DEIM algorithm 
% by (Chaturantabut and Sorensen) and V is assumed to be orthonormal 
% % Input:
%   V, an m by n real orthonormal matrix 
%   
% Output:
%   p, the indices of the n selected rows of V

    k = size(V,2);
    p = zeros(k,1);
    [~,p1] = max(abs(V(:,1)));  p(1) = p1;
    for l = 2:k
        vl = V(:,l);
        plm1 = p(1:l-1);
        c = V(plm1,1:l-1)\vl(plm1);
        r = vl - V(:,1:l-1)*c;

        [~,pl] = max(abs(r));
        p(l) = pl;
    end

end