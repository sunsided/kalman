real_frequency = .42;         % Hz
real_phase     = deg2rad(15); % radians
real_amplitude = 2;           % <scalar>
real_offset    = -5;          % <scalar>

T_start        = 0;           % seconds
T_end          = 20;          % seconds
N_samples      = 2000;        % <number of samples>

% generate the time vector
time_vector = linspace(T_start, T_end, N_samples);
 
% generate the data vector
real_omega = 2*pi*real_frequency;
real_data = real_amplitude*sin(real_omega*time_vector+real_phase) + real_offset;

% generate the output buffers
observed_data = nan(size(real_data));
estimated_data = nan(size(real_data));

% plot the real data
close all;
plot(time_vector, real_data);
xlabel('t [s]');
ylabel('a*sin(\omegat+\phi)+b');

% "it it is easier to approximate
%  a probability distribution than it is to approximate
%  an arbitrary nonlinear function or transformation"
% J. K. Uhlmann, “Simultaneous map building and localization for
% real time applications,” transfer thesis, Univ. Oxford, Oxford, U.K.,
% 1994.

% set initial state estimate
x = [1;  % frequency [Hz]
     0;  % phase [rad]
     1;  % amplitude [<scalar>]
     0]; % offset [<scalar>]
 
% set initial state covariance
P = 1000*diag(ones(size(x)));
 
% define additive state covariance prediction noise
R = 0 * ...
    [1 0 0   0;     % frequency may change over time
     0 1 0   0;     % phase may change over time
     0 0 0.01 0;     % amplitude does not change
     0 0 0   0.01];  % offset does not change

% define additive measurement covariance prediction noise
z_sigma = .2;
Q = ones(size(1))*z_sigma;

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
    observation_fun = @(x) x(3)*sin(2*pi*x(1)*t+x(2))+x(4);
   
    % time update - propagate state
    [x_prior, P_prior, X, Xwm, Xwc] = unscented(state_transition_fun, ...
                                                x, P, ...
                                                'n_out', 4);
                                            
    % add prediction noise
    P_prior = P_prior + R;
    
    if mod(i,5) ~= 0
        % pass variables around
        x = x_prior;
        P = P_prior;
    else
        % predict observations using the a-priori state
        % Note that the weights calculated by this function are the very
        % same as calculated above since we're still operating on 
        % the state vector (i.e. dimensionality didn't change).
        [z_estimate, S_estimate, Z] = unscented(observation_fun, ...
                                                x_prior, P_prior, ...
                                                'n_out', 1);

        % add measurement noise Q to S_estimate
        S_estimate = S_estimate + Q;

        % calculate state-observation cross-covariance
        Pxy = zeros(numel(x), numel(z_estimate));
        for j=1:numel(Xwc)
            Pxy = Pxy + Xwc(j)*(X(:,j)-x_prior)*(Z(:,j)-z_estimate)';
        end

        % calculate Kalman gain
        K = Pxy/S_estimate; % note the inversion of S!

        % obtain observation
        z_error = z_sigma*randn(1);
        z = real_data(i) + z_error;

        % measurement update
        x_posterior = x_prior + K*(z - z_estimate);
        P_posterior = P_prior - K*S_estimate*K';

        % pass variables around
        estimated_data(i) = z_estimate;
        observed_data(i) = z;
        x = x_posterior;
        P = P_posterior;
    end
 end

% plot the estimated data
hold all;
valid = ~isnan(estimated_data);
plot(time_vector(valid), estimated_data(valid), 'r');
plot(time_vector(valid), observed_data(valid), 'm+');
