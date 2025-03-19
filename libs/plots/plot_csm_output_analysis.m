function plot_csm_output_analysis(S4, psi_csm_samples, S4_differences)
    % Generate visualizations for S4 analysis
    figure;
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % 1. Histogram of S4 differences
    nexttile;
    histogram(S4_differences, 'Normalization', 'pdf', NumBins=50);
    title('Distribution of S4 Differences');
    xlabel('Absolute Difference |S4_{calculated} - S4_{input}|');
    ylabel('Probability Density');
    grid on;

    % 2. Plot of S4 differences over iterations
    nexttile;
    plot(1:length(S4_differences), S4_differences, '-o');
    title('Absolute Differences Over Iterations');
    xlabel('Iteration');
    ylabel('|S4_{calculated} - S4_{input}|');
    grid on;

    % 3. Histogram of |psi_csm| with Rice PDF
    nexttile;
    abs_psi_csm = abs(psi_csm_samples);
    histogram(abs_psi_csm, 'Normalization', 'pdf', 'DisplayName', '|psi_{csm}| Samples', NumBins=50);
    hold on;

    K = sqrt(S4^2 / (1 - S4^2)); % Rice K-factor
    x_values = linspace(0, 2.5, 200);
    rice_pdf = get_rice(x_values,K,1);
    plot(x_values, rice_pdf, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Rice PDF');
    hold off;

    title('Histogram of |psi_{csm}| with Rice PDF');
    xlabel('|psi_{csm}|');
    ylabel('Probability Density');
    legend show;
    grid on;

    % 4. PDF of the angle of psi_csm
    nexttile;
    angles = angle(psi_csm_samples);
    histogram(angles, 'Normalization', 'pdf', 'DisplayName', 'Angle of psi_{csm}', NumBins=50);
    title('PDF of the Angle of psi_{csm}');
    xlabel('Angle (radians)');
    ylabel('Probability Density');
    grid on;
end

function rice_pdf = get_rice(alpha, K, Omega)
% get_rice
% Computes the theoretical Rician PDF as per the article's formula.
%
% Syntax:
%   rice_pdf = get_rice(alpha, K, Omega)
%
% Description:
%   This function calculates the Rician probability density function (PDF)
%   based on the provided parameters, using the formula from the referenced
%   article.
%
% Inputs:
%   alpha - Values at which to compute the PDF (vector or scalar).
%   K     - Rician K-factor (scalar, must be >= 0).
%   Omega - Total power (scalar, must be > 0).
%
% Outputs:
%   rice_pdf - Theoretical Rician PDF values corresponding to input alpha.

% Input validation
validateattributes(alpha, {'numeric'}, {'real', 'nonnegative'}, 'get_rice', 'alpha');
validateattributes(K, {'numeric'}, {'scalar', 'real', 'nonnegative'}, 'get_rice', 'K');
validateattributes(Omega, {'numeric'}, {'scalar', 'real', 'positive'}, 'get_rice', 'Omega');

% Compute the Rician PDF
rice_pdf = (2 * alpha * (1 + K) / Omega) .* ...
           besseli(0, 2 * alpha .* sqrt(K + K^2) / Omega) .* ...
           exp(-K - alpha.^2 * (1 + K) / Omega);
end
