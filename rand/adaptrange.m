function [Q,errest,iter] = adaptrange(A, maxiter, nb, tol, varargin)
	% Adaptively estimate the matrix range
    % based on Algorithm 4.2 HMT Siam Review 2011. 
    % 
	% Inputs	
	% A         :   Matrix A of size m x n
	% maxiter   : 	max number of iterations
	% nb        : 	block size for incrementing
	% tol       :	relative error norm that is allowable	
	% 
	% Outputs
	% Q         : 	approximation of matrix range
	% errest    : 	randomized estimate of ||A - QQ'A||
    % iter      :	number of iterations, app. rank = iter*nb
    % normtype  :   {'two' (default), 'fro'}
    % 
    % Written by A.K. Saibaba 6/2017, Modified 7/2017
    
	[m,n] = size(A);
    
    if nargin < 5
        normtype = 'two';
    else 
        normtype = varargin{1};
    end
    
	% Do one step of subspace iteration to estimate range
	Omega = randn(n,nb);
	Y = A*Omega;
    [Q,~] = qr(Y,0);

    Anrm = norm(A,'fro');
    
    
    % Factor to test the error norm (see Lemma 4.1 HMT SIAM Review, 2011)
	alpha = exp(log(1./0.01)/nb);		% determine alpha based on 
	fact = alpha*sqrt(2/pi);            % 1% failure probability 

	iter = maxiter;
	for i = 1:maxiter
        
		% Estimate the error
		Omega = randn(n,nb);            % generate new samples 
		Y = A*Omega;
	
        Ynrm = sqrt(sum(Y.^2,1));
		res = Y - Q*(Q'*Y);             % compute residual
		resnrm = sqrt(sum(res.^2,1));	% residual columnwise 2-norm
		
        switch normtype% estimate of the error
            case 'two'
                errest = max(resnrm)*fact/max(Ynrm);
            case 'fro'
                errest = sqrt(min(m,n))*max(resnrm)*fact/Anrm;        
        end
        if errest < tol                 % break if error is below tolerance
			iter = i;
			break;
        end	
        
        % Extend the basis otherwise and keep going
        Q = orth([Q, Y]);			
	end


end
