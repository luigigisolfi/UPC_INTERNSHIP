function STM = new_numericSTMvfield(ti, tf, xi, eps, vfield, hmin, hmax, tol)
%------------------------------------------------------------------------
% Computes the State Transition Matrix (STM) from ti to tf of an initial
% condition xi in the vector field vfield using finite differences
%------------------------------------------------------------------------

n = length(xi); % Determine the dimension of the system
STM = zeros(n, n); % Preallocate the STM matrix

parfor k = 1:n % Use parfor for parallel computing
    % Positive perturbation
    xa_plus = xi;
    xa_plus(k) = xi(k) + eps;
    xfp = new_propTITF_vfield(ti, xa_plus, tf, vfield, hmin, hmax, tol);

    % Negative perturbation
    xa_minus = xi;
    xa_minus(k) = xi(k) - eps;
    xfm = new_propTITF_vfield(ti, xa_minus, tf, vfield, hmin, hmax, tol);

    % Compute the k-th column of the STM
    STM(:, k) = (xfp - xfm) / (2 * eps);
end
end
