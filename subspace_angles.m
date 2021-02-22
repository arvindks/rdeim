function theta = subspace_angles(u,uh)
    % Compute the angles b/w subspaces when matrices have 
    % orthonormal columns
    % 
    % Inputs
    % U      :  n x k  matrix with orthonormal columns
    % Uh     :  n x k' matrix with orthonormal columns
    % 
    % Outputs
    % theta  : min{k,k'} x 1 canonical angles b/w subspacea
    % 
    A = u'*uh;
    s = svd(A,0);
    theta = real(acos(s));
end