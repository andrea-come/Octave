function res = my_polyarea(xs, ys)
% Computes the area of a 2D polygon given its vertices.
%  implementing the Gauss Area Formula (Shoelace Formula).
%
% The input vectors must satisfy the following conditions:
%       1. They must have the same length.
%       2. The vertices must be ordered consecutively (clockwise or 
%          counter-clockwise).
%
%   Inputs:
%       xs - Vector of x-coordinates (double).
%       ys - Vector of y-coordinates (double).
%
%   Output:
%       res - The calculated area (positive scalar).
%
  if nargin < 2
    error("Two input arguments are required: xs and ys.");
  endif
  x_len = length(xs); y_len = length(ys);
  if x_len ~= y_len
    error("The input vectors must have same length")
  endif
  xs = xs(:); ys = ys(:);
  n = x_len;
  res = 0.5 * abs( sum(xs .* circshift(ys, -1) - ys .* circshift(xs, -1)));
endfunction
