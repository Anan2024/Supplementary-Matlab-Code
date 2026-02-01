function UCL = bisection_method_to_refine_UCL(MRL_0, r, T_0, m, n, lambda, d,p)

% Additional parameters
max_iterations = 5000;  % Maximum iterations for bisection
iteration = 0;
tolerance = 0.01;     % Relative tolerance for convergence
converged = false;

 % Initial guesses
    H1 = 1;  % Small initial value
    H2 = 5;   % Large initial value

    MRL_H1 = calculate_IC_MRL(H1, lambda, p, r, m, n, T_0, d); 
    MRL_H2 = calculate_IC_MRL(H2, lambda, p, r, m, n, T_0, d);

    while iteration < max_iterations && ~converged
        iteration = iteration + 1;
        H = (H1 + H2) / 2;
        MRL = calculate_IC_MRL(H, lambda, p, r, m, n, T_0, d)

        rel_error = (MRL - MRL_0) / MRL_0;

        if abs(rel_error) <= tolerance
            UCL = H;
            converged = true;
        else
            % Update bounds based on relative error
            if rel_error > 0
                % MRL too large, need to decrease UCL
                H2 = H;
            else
                % MRL too small, need to increase UCL
                H1 = H;
            end
        end
    end
    fprintf('\nFinal result: UCL = %.6f\n', UCL);
    