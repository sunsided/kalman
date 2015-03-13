function [my, sigma] = unscented(fun, my, sigma)
    % scaled unscented transform with 2*n+1 sigma points

    % enforce the canonical row vector
    if isrow(my)
        my = my';
    end
    
    % obtain the number of states
    n = numel(my);

    % determine free parameters
    % http://ais.informatik.uni-freiburg.de/teaching/ws12/mapping/pdf/slam05-ukf.pdf
    kappa = 1;
    alpha = 1; % 1e-3;
    beta = 2;
    lambda = alpha^2*(n+kappa) - n;

    % calculate matrix square root of adjusted covariance matrix
% TODO: compute via diagonalization
    S_adjusted = alpha*sqrt(n+kappa) * sqrtm(sigma);
        
    % calculate sigma points
    X = nan(n,2*n+1);      % prepare array of vectors
    X(:,1) = my;            % bootstrap first sigma point
    for i=1:n
        X(:,i+1)   = my + S_adjusted(:,i);
        X(:,i+n+1) = my - S_adjusted(:,i);
    end
    
    % transform the calculated sigma points
    Y = nan(size(X));
    for i=1:size(X,2)
        Y(:,i) = fun(X(:,i));
    end
    
    % calculate weights for expectation determination
    wm0 = lambda/(n+lambda); % = 1 - n/(alpha^2*(n+kappa))
    wmi = 0.5/(n+lambda);
    
    % calculate weights for covariance determination
    wc0 = wm0 + (1 - alpha^2 + beta);
    wci = wmi;
    
    % calculate the new expectation
    my = wm0*Y(:,1);
    for i=2:2*n+1
        my = my + wmi*Y(:,i);
    end
    
    % calculate the new covariance
    sigma = wc0*(Y(:,1)-my)*(Y(:,1)-my)';
    for i=2:2*n+1
        sigma = sigma + wci*(Y(:,i)-my)*(Y(:,i)-my)';
    end
end