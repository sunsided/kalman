function [my, sigma, Y, wm, wc] = unscented(fun, my, sigma, varargin)
    % scaled unscented transform with 2*n+1 sigma points

    % enforce the canonical row vector
    if isrow(my)
        my = my';
    end

    % default values
    % http://ais.informatik.uni-freiburg.de/teaching/ws12/mapping/pdf/slam05-ukf.pdf
    % http://ei.uni-paderborn.de/fileadmin/Elektrotechnik/FG-NTH/http/download/estimation.pdf
    defaultKappa = 1;
    defaultAlpha = 1; % 1e-3 is typical, 1 results in regular unscented transform
    defaultBeta = 2;
    
    % number of output variables
    defaultNout = 0;
    
    % parse inputs
    p = inputParser;
    addRequired(p, 'fun', @(fh) isa(fh,'function_handle'));
    addRequired(p, 'my', @isvector);
    addRequired(p, 'sigma', @ismatrix);
    addOptional(p, 'kappa', defaultKappa, @isnumeric);
    addOptional(p, 'alpha', defaultAlpha, @isnumeric);
    addOptional(p, 'beta', defaultBeta, @isnumeric);
    addOptional(p, 'n_out', defaultNout, @isnumeric);
    parse(p, fun, my, sigma, varargin{:});   
    
    % obtain the number of states
    n = numel(my);

    % determine free parameters
    kappa = p.Results.kappa;
    alpha = p.Results.alpha;
    beta = p.Results.beta;
    n_out = p.Results.n_out;
    lambda = alpha^2*(n+kappa) - n;

    % calculate matrix square root of adjusted covariance matrix
    gamma = sqrt(n+lambda);
% TODO: compute via diagonalization
    S_adjusted = gamma * sqrtm(sigma);
        
    % calculate sigma points
    X = nan(n,2*n+1);      % prepare array of vectors
    X(:,1) = my;            % bootstrap first sigma point
    for i=1:n
        X(:,i+1)   = my + S_adjusted(:,i);
        X(:,i+n+1) = my - S_adjusted(:,i);
    end
    
    % prepare the array for the sigma points. 
    % Since we can't determine it from the state vector alone (the 
    % evaluation function might return a vector of any dimension) we use
    % the hint given in the function call or fall back to a dynamically 
    % growing array if no hint was given.
    if n_out ~= defaultNout
        Y = nan(n_out, size(X,2));
    else
        Y = [];
    end
    
    % transform the calculated sigma points
    for i=1:size(X,2)
        Y(:,i) = fun(X(:,i));
    end
    
    % calculate weights for expectation determination
    wm0 = lambda/(n+lambda); % = 1 - n/(alpha^2*(n+kappa))
    wmi = 0.5/(n+lambda);
    
    % calculate weights for covariance determination
    wc0 = wm0 + (1 - alpha^2 + beta);
    wci = wmi;

    % bundle for output
    wm = [wm0, wci*ones(1,2*n)];
    wc = [wc0, wci*ones(1,2*n)];
    
    % calculate the new expectation
    my = 0;
    for i=1:2*n+1
        my = my + wm(i)*Y(:,i);
    end
    
    % calculate the new covariance
    sigma = zeros(size(sigma));
    for i=1:2*n+1
        sigma = sigma + wc(i)*(Y(:,i)-my)*(Y(:,i)-my)';
    end
end