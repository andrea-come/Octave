% Exercise 2.18: Fit sines to straight line

classdef fit_sines
  methods (Static)

    function S = sinesum(t, b)
    %SINESUM Computes the partial sum of a Fourier sine series.
    %  S = fit_sines.sinesum(t, b) evaluates the sum:
    %      S_N(t) = sum_{n=1}^{N} b_n * sin(n * t)
    %
    %  Inputs:
    %      t - Vector of time coordinates (scalar or vector).
    %      b - Vector of coefficients [b_1, b_2, ..., b_N].
    %
    %  Output:
    %      S - Vector containing the evaluated sum. It preserves the
    %          shape (row or column) of the input vector 't'.

      % Input Shape Handling
      % Check if input 't' is a row vector before reshaping it.
      % This flag is used later to ensure the output matches the input shape.
      input_is_row = (isvector(t) && size(t, 2) > 1);

      % Force inputs to be column vectors for matrix operations
      t = t(:); b = b(:);

      % Matrix Construction
      N = length(b); ns = 1 : N;

      % Compute the sine matrix using the Outer Product.
      % Dimensions: (M x 1) * (1 x N) -> (M x N) matrix.
      % Columns represent harmonics (n), Rows represent time steps (t).
      sinM = sin(t * ns);

      % Computation via Linear Combination
      % Multiply the matrix by the coefficient vector.
      % Dimensions: (M x N) * (N x 1) -> (M x 1) column vector.
      S = sinM * b;

      % Output Formatting
      % If the original input 't' was a row vector, transpose the result.
      if input_is_row
        S = S.';
      endif
    endfunction


    function test_sinesum()
    %TEST_SINESUM Verifies the correctness of the sinesum function.
    %  The function compares the computed result against the exact
    %  analytical solution for a specific test case:
    %    t = [-pi/2, pi/4]
    %    b = [4, -3] (corresponding to 4*sin(t) - 3*sin(2t))

      % 1. Define test parameters
      t = [-pi/2, pi/4]; b = [4, -3];

      % 2. Compute the result using the function under test
      computed_S = fit_sines.sinesum(t, b);

      % 3. Define the exact expected result (analytical solution)
      %    using a SLOW but SAFE method (loop)
      expected_S = zeros(size(t));
      for k = 1:length(b)
          expected_S = expected_S + b(k) * sin(k * t)
      endfor

      % 4. Verify the result with a tolerance
      tol = 1e-12;
      max_diff = max(abs(computed_S(:) - expected_S(:)));

      if max_diff < tol
        disp('Test passed: The calculated values match the expected ones.');
      else
        error('Test failed: Max difference is %e.', max_diff);
      endif
    endfunction


    function plot_compare(f, b, M)
    %PLOT_COMPARE Plots the original function f(t) and the sine approximation
    % S_N(t).
    %
    % Inputs:
    %   f - Handle to the function to approximate (e.g., @(t) t).
    %   b - Vector of coefficients [b_1, b_2, ..., b_N]..
    %   M - Number of plotting points (coordinates).

      % Standard interval for Fourier series is usually [-pi, pi]
      t = linspace(-pi, pi, M)';
      y_exact = f(t);
      y_approx = fit_sines.sinesum(t, b);
      n = length(b);

      % `b`  can be calculated using Least Squares
      %S_mat = sin(t * (1 : n));
      %b = S_mat \ y_exact;

      figure(1); % Open a new figure window
      plot(t, y_exact, 'b-', 'LineWidth', 2, 'DisplayName', 'Original f(t)');
      hold on;
      plot(t, y_approx, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('S_N(t) with N=%d', n));

      xlabel('Time t');
      ylabel('Function value');
      title(sprintf('Approximation quality check', n));
      legend('Location', 'northeast');
      grid on;
      hold off;
    endfunction


    function approx_err = compute_error(b, f, M)
      % COMPUTE_ERROR Calculates the root sum of squared errors.
      %   E = fit_sines.compute_error(b, f, M) computes the error measure:
      %   E = sqrt( sum( (f(t_i) - S_N(t_i))^2 ) )
      %   where t_i are M uniformly distributed points in [-pi, pi].
      %
      %   Inputs:
      %       b - Vector of coefficients for the sine sum.
      %       f - Function handle for the target function f(t).
      %       M - Number of sample points.
      t = linspace(-pi, pi, M)';  % Force column vector for safety
      y_exact = f(t);
      y_approx = fit_sines.sinesum(t, b);
      approx_err = sqrt(sum((y_exact - y_approx).^2));
    endfunction


    function trial(f, N)
    %TRIAL Interactive tool to test Fourier coefficients manually.
    %   Loop asks for coefficients one by one, plots the result,
    %   and computes the error.
      M = 500;
      bn = zeros(1, N);
      t = linspace(-pi, pi, M)';

      new_trial = true;
      while new_trial
        printf("Enter %d values:\n", N)
        for k = 1 : N
          msg = sprintf("Insert value b_%d: ", k);
          val = input(msg);
          if isempty(val)
            val = 0;
          endif
          bn(k) = val;
        endfor

        fit_sines.plot_compare(f, bn, M);
        approx_err = fit_sines.compute_error(bn, f, M);
        printf("The error is: %f\n", approx_err);

        answer = input("Another trial (y/n)? ", "s");
        if strcmpi(answer, 'n')
          new_trial = false;
          fprintf('Exiting trial.\n');
        endif
      endwhile
    endfunction


    function [best_b, min_err] = automatic_trial(f)
    %AUTOMATIC_TRIAL Brute-force search for the best 3 coefficients.
    %   Iterates through b1, b2, b3 in the range [-1, 1] with step 0.1
    %   to find the combination that minimizes the approximation error.
      M = 500;
      t = linspace(-pi, pi, M)';
      % Pre-compute the sine basis matrix to optimize performance.
      % This avoids re-evaluating the sine function ~9000 times inside the loop.
      ns = 1:3;
      S_mat = sin(t * ns);

      y_exact = f(t);
      best_b = [0, 0, 0];
      min_err = inf;

      % Iterate using integers [-10, 10] to prevent floating-point precision
      range_vals = -10:1:10;
      for b1 = range_vals
        for b2 = range_vals
          for b3 = range_vals
            % Normalize integer indices to get coefficients in range [-1.0, 1.0]
            b = [b1 b2 b3] / 10.0;
            % Perform direct matrix-vector multiplication instead of calling
            % external functions to minimize overhead.
            % Note: b' turns the row vector into a column for the product.
            y_approx = S_mat * b';
            approx_err = sqrt(sum((y_exact - y_approx).^2));
            if approx_err < min_err
              best_b = b;
              min_err = approx_err;
            endif
          endfor
        endfor
      endfor

      % Print results to screen
      fprintf('\n--- Optimization Complete ---\n');
      fprintf('Smallest error found: %f\n', min_err);
      fprintf('Best coefficients: b = [%.1f, %.1f, %.1f]\n', best_b);

      % Plot the best result found
      fit_sines.plot_compare(f, best_b, M);
    endfunction


    function [best_b, min_err] = fast_automatic_trial(f)
    %AUTOMATIC_TRIAL_VECTORIZED Fully vectorized brute-force search.
    %   Finds the best coefficients without using any for-loops.
    %   This implementation leverages matrix algebra for maximum performance.

      M = 500;
      t = linspace(-pi, pi, M)';

      % Pre-compute the Sine Basis Matrix (Mx3)
      ns = 1:3;
      S_mat = sin(t * ns);

      % Generate the Coefficients Grid, creating a grid of integers [-10, 10]
      vals = -10:1:10;
      [G1, G2, G3] = ndgrid(vals, vals, vals);

      % Flatten the 3D grids into a single matrix of coefficients (Kx3)
      % where K = 21*21*21 = 9261.
      % Divide by 10 to map integers back to [-1.0, 1.0].
      B_all = [G1(:), G2(:), G3(:)] / 10.0;

      % Vectorized Approximation Calculation
      % Multiply Sine Matrix (Mx3) by Transposed Coeffs (3xK).
      % Result Y_all is (MxK): each column is one trial approximation.
      Y_all = S_mat * B_all';

      % Vectorized Error Calculation
      y_true = f(t);

      % Compute difference (using broadcasting: Mx1 vs MxK)
      Diffs = y_true - Y_all;

      % Calculate Euclidean norm for each column (Result is 1xK)
      % sum(..., 1) sums along the rows (collapsing them)
      Errors = sqrt(sum(Diffs.^2, 1));

      % 5. Find the Minimum
      [min_err, min_idx] = min(Errors);

      % Retrieve the best coefficient set using the index
      best_b = B_all(min_idx, :);

      % Output results
      fprintf('\n--- Vectorized Optimization Complete ---\n');
      fprintf('Checked %d combinations.\n', length(Errors));
      fprintf('Smallest error: %f\n', min_err);
      fprintf('Best b: [%.1f, %.1f, %.1f]\n', best_b);

      fit_sines.plot_compare(f, best_b, M);
    endfunction

  endmethods
endclassdef
