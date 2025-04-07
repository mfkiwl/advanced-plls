function [x_est_ds, u_est_ds, cycleSlipCount] = get_cycle_slips(y, lambda, dsFactor)
% get_cycle_slips Detect cycle slips using TVD with downsampling.
%
%   [x_est_ds, u_est_ds, cycleSlipCount] = get_cycle_slips(y, lambda, dsFactor)
%   takes the observed signal y, regularization parameter lambda, and an
%   optional downsampling factor dsFactor. It returns:
%
%       - x_est_ds: the denoised (piecewise constant) version of the downsampled signal.
%       - u_est_ds: the full solution [x; z], where z approximates |diff(x)|.
%       - cycleSlipCount: the number of cycle slips detected by thresholding z.
%
%   The optimization problem is:
%
%       min_{x,z} 0.5*||x - y_ds||^2 + lambda * sum(z)
%       s.t.   x(i+1) - x(i) - z(i) <= 0,
%              -[x(i+1) - x(i)] - z(i) <= 0,  for i = 1,...,n-1,
%              z >= 0.
%
%   y_ds is the downsampled version of y (if dsFactor > 1, otherwise y_ds = y).
%   Cycle slips are detected by counting the number of z values above a threshold.
%
%   If dsFactor is not provided, it defaults to 1 (no downsampling).

    % Default downsampling factor if not provided
    if nargin < 3
        dsFactor = 1;
    end

    % Downsample the signal if needed
    if dsFactor > 1
        y_ds = downsample(y, dsFactor);
    else
        y_ds = y;
    end

    % Define sizes based on downsampled signal
    n = length(y_ds);
    m = n - 1;          % number of differences
    N = n + m;          % total variables: [x; z]

    % Initial guess: x = y_ds, z = abs(diff(y_ds))
    x0 = y_ds;
    z0 = abs(diff(y_ds));
    u0 = [x0; z0];

    % Lower bounds: no bounds for x, z >= 0.
    lb = [-inf(n,1); zeros(m,1)];
    ub = [];

    % Build the difference operator for x as an m x n sparse matrix:
    E = spdiags([-ones(m,1), ones(m,1)], [0, 1], m, n);
    
    % Build inequality constraints in sparse form:
    % Constraint 1: [E, -I] * u <= 0
    % Constraint 2: [-E, -I] * u <= 0
    A = [E, -speye(m); -E, -speye(m)];
    b = zeros(2*m, 1);

    % Define the objective function
    function f = obj(u)
        x = u(1:n);
        z = u(n+1:end);
        f = 0.5 * sum((x - y_ds).^2) + lambda * sum(z);
    end

    % Set up fmincon options with Hessian approximation
    options = optimoptions('fmincon',...
        'Display', 'iter',...
        'Algorithm', 'interior-point',...
        'HessianApproximation', 'lbfgs');

    % Solve the optimization problem
    u_est_ds = fmincon(@obj, u0, A, b, [], [], lb, ub, [], options);
    x_est_ds = u_est_ds(1:n);
    
    % Extract the auxiliary variable z
    z_est = u_est_ds(n+1:end);
    
    % Cycle-slip detection: count indices where z exceeds a threshold.
    % (Adjust the threshold value based on your application's noise level.)
    threshold = 1.5;  
    cycleSlipCount = sum(z_est > threshold);
end
