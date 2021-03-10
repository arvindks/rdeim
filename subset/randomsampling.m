function [ind, SQ] = randomsampling(Q, c, varargin)
    % Code from the kappaSQ package:https://arxiv.org/pdf/1402.0642.pdf
    
    params = inputParser;
    params.addParameter('sampling','Bernoulli');
    params.KeepUnmatched = true;
    params.parse(varargin{:});

    
    sampling = params.Results.sampling;
    m = size(Q,1);
    
    if strcmp(sampling,'Bernoulli')
        ind = find(rand(1,m)<c/m);
    elseif strcmp(sampling,'randperm')
        v = randperm(m);    ind = v(1:ceil(c));
    elseif strcmp(sampling,'exactly')
        ind = randsample(m,ceil(c));    
    end
   
    SQ = sqrt(m/c)*Q(ind,:);
    
end