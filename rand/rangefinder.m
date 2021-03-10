function Q = rangefinder(A, k, p, q)

    % Inputs
    % A       : (m x n) matrix of interest
    % k       : (int => 1) target rank
    % p       : (int => 0) oversampling parameter
    % q       : (int => 0) number of subspace iterations
    % 
    % Outputs
    % Q       : (m x (k+p)) randomized basis 
    
    [~,n] = size(A);
    Omega = randn(n,k+p);       % Draw random matrix
    Y     = A*Omega;            % Sketch
    [Q,~] = qr(Y, 0);           % Ortho basis
    
    
    % subspace iterations
    for i = 1:q
       Y = A'*Q; 
       [Q, ~] = qr(Y, 0);
       
       Y = A*Q;
       [Q, ~] = qr(Y, 0);
    end
    
    B = Q'*A;
    [u,~,~] = svd(B, 0);
    Q = Q*u(:,1:k);
end