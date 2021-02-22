function [selected_indices, selected_columns] = srrqr(M, k)
%
% This function determines the indices for a "good" set of k columns
% from an input matrix M by 
% using the Gu-Eisenstat strong rank revealing QR factorization.
%
% Input:
%   M, an m by n real matrix of rank >= k
%   k, the number of columns to choose
%
% Output:
%   selected_indices, the indices of the k selected columns of M
%   selected_columns, the k selected columns of M
%
        
% Initializations
  [m, n] = size(M);
  P = 1:n;                 % P keeps track of column permutations
  [Q, R, P] = qr(M, 0);  % Golub QR with column pivoting
  % this guarantees that leading k columns of R are lin. indep.
  increasefound = true;
  counter_perm=0;

    while (increasefound)
        A = R(1:k,1:k); 
        AinvB = A\R(1:k,k+1:n); % Form A^{-1}B
        if m<=n C = R(k+1:m, k+1:n);
        else C = R(k+1:n, k+1:n);
        end % if
        
        %compute the column norms of C
        gamma = zeros(n-k, 1);
        for ccol = 1:n-k
            gamma(ccol) = norm(C(:,ccol), 2);
        end
        
        %find row norms of A^-1
        [U, S, V] = svd(A);
        Ainv = V*diag(1./diag(S))*U';
        omega = zeros(k, 1);
        for arow = 1:k
            omega(arow) = norm(Ainv(arow,:), 2);
        end
        
        %find indices i and j that maximize 
        %ainv(i,j)^2 + (w(i)*gamma(j))^2
        tmp = omega*gamma';
        F = AinvB.^2 + tmp.^2;
        [i, j] = find(F>1, 1);
        if (isempty(i))           % finished
            increasefound = false;
        else  %we can increase |det(A)|
         counter_perm = counter_perm +1;
            R(:,[i j+k]) = R(:, [j+k i]);  % permute columns i and j
            P([i j+k]) = P([j+k i]);
            [Q, R] = qr(R, 0); % retriangularize R
        end   %if
    end   %while

norm_of_inverse = norm(inv(R(1:k,1:k)));
if m<=n, residual_norm = norm(R(k+1:m,k+1:n));
else residual_norm = norm(R(k+1:n,k+1:n));
end   %if
selected_indices = sort(P(1:k),'ascend');
selected_columns = M(:,selected_indices);