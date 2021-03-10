function [p,d,err] = randls(V,varargin)
    % Code from the kappaSQ package:https://arxiv.org/pdf/1402.0642.pdf

    params = inputParser;
    params.addParameter('epsilon', 0.99, @(x) x > 0); %Accuracy
    params.addParameter('delta', 0.1, @(x) x > 0); %Failure probability
    params.addParameter('method','default'); %
    params.addParameter('factor',3); %factor
    params.KeepUnmatched = true;
    params.parse(varargin{:});
    
    
    epsilon  = params.Results.epsilon;
    delta    = params.Results.delta;
    method   = params.Results.method;
    fact     = params.Results.factor;
    
    [m,n] = size(V);
    
    li = leverageScores(V);
    if strcmp(method,'default')
        %Determining number of samples required
        
        mu = max(li);	%Compute coherence
        t = floor(1./mu);

        [li,~] = sort(li, 'descend');
        tau = mu*sum(li(1:t)) + (1-t*mu)*li(t+1);	

        c = (2/3)*m*(3*tau+epsilon*mu)*log(2*n/delta)/(epsilon^2);
        c = ceil(c);
    
    elseif strcmp(method, 'const')
        c = ceil(fact*n*log(n));
    end
    
    %Extracting the relevant rows of B
    [p, SV] = sample_leverage(V,c,li); 
    
    d   = 1./sqrt(c*li(p));
    err = norm(pinv(SV))*max(d);
end 

function li = leverageScores(Q)
    li = sum(Q.^2, 2);
end

function [p,SQ] = sample_leverage(Q,c,li)
    [m,n] = size(Q);
    p = randsample(m,c,true,li/n); % We are setting beta = 1 here.
    SQ = (1/sqrt(c/m))*Q(p,:);
end

