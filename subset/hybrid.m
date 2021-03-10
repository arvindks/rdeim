function [p,pls,err,d] = hybrid(V,varargin)

    % Reasonable selection of parameters.
    params = inputParser;
    params.addParameter('epsilon', 0.99, @(x) x > 0); %Accuracy
    params.addParameter('delta', 0.1, @(x) x > 0); %Failure probability
    params.addParameter('det', 'rrqr'); 
    params.addParameter('method','default'); %
    params.addParameter('factor',3); %factor
    params.KeepUnmatched = true;
    params.parse(varargin{:});
    
    
    epsilon  = params.Results.epsilon;
    delta    = params.Results.delta;
    det      = params.Results.det;
    method   = params.Results.method;
    factor   = params.Results.factor;
    
    % Stage 1: Leverage Score approach
    [pls,d,~] = randls(V,'epsilon',epsilon,'delta',delta,'method',method,'factor',factor);
    Vls = diag(d)*V(pls,:);

    % Stage 2: Deterministic subset selection
    [p,~] = subsetselection(Vls,det);
    d = d(p);
    p = pls(p);
    
    % Compute error in DEIM approximation
    err = norm(pinv(diag(d)*V(p,:))*diag(d));
    
end
