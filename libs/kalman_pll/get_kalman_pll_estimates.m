function [state_estimates, error_covariance_estimates] = ...
    get_kalman_pll_estimates(received_signal,kalman_pll_config,initial_estimates,training_scint_model)
    % Update Step
    steps = 1:1:length(received_signal);
    state_estimates = zeros( ...
        length(received_signal), ...
        length(initial_estimates.x_hat_init) ...
        );
    state_estimates(1,:) = initial_estimates.x_hat_init.';
    error_covariance_estimates = zeros( ...
        length(received_signal), ...
        size(initial_estimates.P_hat_init,1), ...
        size(initial_estimates.P_hat_init,2) ...
    );
    x_hat_project_ahead = initial_estimates.x_hat_init;
    P_hat_project_ahead = initial_estimates.P_hat_init;
    F = kalman_pll_config.(training_scint_model).F;
    Q = kalman_pll_config.(training_scint_model).Q;
    H = kalman_pll_config.(training_scint_model).H;
    R = kalman_pll_config.(training_scint_model).R;
    W = kalman_pll_config.(training_scint_model).W;
    for step = steps
        if step > 1
            K = P_hat_project_ahead * H.' * (H * P_hat_project_ahead * H.' + R);
            x_hat_update = x_hat_project_ahead + K * angle(received_signal(step,1) .* exp(-1j * H * x_hat_project_ahead));
            P_hat_update = P_hat_project_ahead - K * H * P_hat_project_ahead;
        else
            x_hat_update = x_hat_project_ahead;
            P_hat_update = P_hat_project_ahead;
        end
        x_hat_project_ahead = F * x_hat_update + W;
        P_hat_project_ahead = F * P_hat_update * F.' + Q;
    end
end