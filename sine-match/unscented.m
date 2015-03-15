function [my, sigma, Y, wm, wc] = unscented(fun, my, sigma, varargin)
    % scaled unscented transform with 2*n+1 sigma points

    % enforce the canonical row vector
    if isrow(my)
        my = my';
    end

    % default values
    % http://ais.informatik.uni-freiburg.de/teaching/ws12/mapping/pdf/slam05-ukf.pdf
    % http://ei.uni-paderborn.de/fileadmin/Elektrotechnik/FG-NTH/http/download/estimation.pdf
    defaultKappa = 0;
    defaultAlpha = 1e-3; % is given as typical, 1 results in regular unscented transform
    defaultBeta = 2;
    
    % number of output variables
    defaultNout = 0;
    
    % default constraint function
    defaultConstraintFun = @(x) x;
    
    % parse inputs
    p = inputParser;
    addRequired(p, 'fun', @(fh) isa(fh,'function_handle'));
    addRequired(p, 'my', @isvector);
    addRequired(p, 'sigma', @ismatrix);
    addOptional(p, 'kappa', defaultKappa, @(x) isnumeric(x) && (x >= 0));
    addOptional(p, 'alpha', defaultAlpha, @(x) isnumeric(x) && (x > 0) && (x <= 1));
    addOptional(p, 'beta', defaultBeta, @isnumeric);
    addOptional(p, 'n_out', defaultNout, @isnumeric);
    addOptional(p, 'constraint', defaultConstraintFun, @(fh) isa(fh,'function_handle'));
    parse(p, fun, my, sigma, varargin{:});   
    
    % obtain the number of states
    n = numel(my);

    % determine constraint function
    constraint = p.Results.constraint;
    
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
    X(:,1) = my;           % bootstrap first sigma point
    for i=1:n
        X(:,i+1)   = my + S_adjusted(:,i);
        X(:,i+n+1) = my - S_adjusted(:,i);
    end
    
    % prepare the array for the transformed sigma points. 
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
        
        % apply constraints to the transformed sigma points
        % as described in "Constrained State Estimation Using the 
        % Unscented Kalman Filter" by Kandepu, Imsland and Foss
        Y(:,i) = constraint(Y(:,i));
    end
    
    % make sure we're not doing anything stupid below
    clearvars X;
    
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
    sigma = zeros(numel(my), numel(my));
    for i=1:2*n+1
        sigma = sigma + wc(i)*(Y(:,i)-my)*(Y(:,i)-my)';
    end
end