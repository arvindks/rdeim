function [p,pls,err] = hybrid(V,varargin)

    
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
    
    [pls,d,~] = randls(V,'epsilon',epsilon,'delta',delta,'method',method,'factor',factor);
    Vls = diag(d)*V(pls,:);

    [p,~] = subsetselection(Vls,det);
    p = pls(p);
   
    err = norm(pinv(V(p,:)));
    
end