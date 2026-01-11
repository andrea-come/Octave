% Exercise 2.15: Linear interpolation (Optimized Vector Version)
% This program uses linear interpolation to compute y between two consecutive
% measurements using efficient built-in functions for index search.

% 1. Data Definition
ys = [4.4, 2.0, 11.0, 21.5, 7.5];
N = length(ys);
ts = 0:(N - 1);

% Value to interpolate
t = 2.5;

% 2. Index Search using histc (Efficient)
% [~, idx] returns the bin index, `~` ignore the counts output.
% histc returns 0 if out of range, or N if t matches the very last element.
[~, idx] = histc(t, ts);

% 3. Boundary Checks & Edge Case Handling
if idx == 0
  if t < ts(1)
    error("Value %f is out of range (too low). Min value is %f.", t, ts(1));
  else
    error("Value %f is out of range (too high). Max value is %f.", t, ts(N));
  end
elseif idx == N
  % Special case: t is exactly the last stored time, use the last valid
  % interval [N-1 N] for the interpolation formula.
  idx = N - 1;
end

% 4. Linear Interpolation Calculation
% Formula: y = y_i + (t - t_i) * slope
dt = ts(idx + 1) - ts(idx);       % Width of the time interval
dy = ys(idx + 1) - ys(idx);       % Variation of y
slope = dy / dt;                  % Gradient

y = ys(idx) + (t - ts(idx)) * slope;

% 5. Output Result
fprintf("The interpolated value is y = %.4f\n", y);

% 6. Visualization
figure(1); clf;

% Known data points connected by lines
plot(ts, ys, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Known Data');
hold on;

% Onterpolated point
plot(t, y, 'r*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Interpolated Value');

% Graph formatting
grid on;
xlabel('Time (t)');
ylabel('Value (y)');
title(sprintf('Linear Interpolation at t = %.2f', t));
legend('Location', 'best');
hold off;
