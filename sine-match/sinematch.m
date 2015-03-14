real_frequency = .42;         % Hz
real_phase     = deg2rad(15); % radians
real_amplitude = 2;           % <scalar>
real_offset    = -5;          % <scalar>

T_start        = 0;           % seconds
T_end          = 10;          % seconds
N_samples      = 1000;        % <number of samples>

% generate the time vector
time_vector = linspace(T_start, T_end, N_samples);

% generate the data vector
real_omega = 2*pi*real_frequency;
real_data = real_amplitude*sin(real_omega*time_vector+real_phase) + real_offset;

% plot the real data
close all;
plot(time_vector, real_data);
xlabel('t [s]');
ylabel('a*sin(\omegat+\phi)+b');

% "it it is easier to approximate
%  a probability distribution than it is to approximate
%  an arbitrary nonlinear function or transformation"
% J. K. Uhlmann, �Simultaneous map building and localization for
% real time applications,� transfer thesis, Univ. Oxford, Oxford, U.K.,
% 1994.

% set initial state estimate
x = [1;  % frequency [Hz]
     0;  % phase [rad]
     1;  % amplitude [<scalar>]
     0]; % offset [<scalar>]
 
% set initial state covariance
P = 1000*diag(ones(size(x)));
 
% simulate
for i=1:numel(time_vector);
    % obtain current time and time delta to last step
    t = time_vector(i);
    if i > 1
        T = t-time_vector(i-1);
    else
        T = 0;
    end
    
    % define the nonlinear state transition function
    state_transition_fun = @(x) x;

    % define the nonlinear observation function
    % note this is time dependent
    observation_fun = @(x) x(3)*sin(2*pi*x(1)*T+x(2))+x(4);
    
    % time update - propagate state
    [x_prior, P_prior, X, Xwm, Xwc] = unscented(state_transition_fun, ...
                                                x, P, ...
                                                'n_out', 4);
    
% TODO: add prediction noise R_t to P_prior
    
    % predict observations using the a-priori state
    % Note that the weights calculated by this function are the very
    % same as calculated above since we're still operating on 
    % the state vector (i.e. dimensionality didn't change).
    [z_estimate, S_estimate, Z] = unscented(observation_fun, ...
                                            x_prior, P_prior, ...
                                            'n_out', 1);
    
% TODO: add measurement noise Q_t to Sy
    
    % calculate state-observation cross-covariance
    Pxy = zeros(numel(x), numel(z_estimate));
    for j=1:numel(Xwc)
        Pxy = Pxy + Xwc(j)*(X(:,j)-x_prior)*(Z(:,j)-z_estimate)';
    end
    
    % calculate Kalman gain
    K = Pxy/S_estimate; % note the inversion of S!
    
    % obtain observation
    z = real_data(i); % ... easy way out, should use more fake data here
    
    % measurement update
    x_posterior = x_prior + K*(z - z_estimate);
    P_posterior = P_prior - K*S_estimate*K';
    
    % pass variables around
    x = x_posterior;
    P = P_posterior;
    
end