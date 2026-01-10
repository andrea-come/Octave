% Exercise 2.9: Find crossing points of two graphs
% Consider two functions f(x) = x and g(x) = x^2 on the interval [-4, 4]
%  Write a program that, by trial and error, finds approximately for which
%  values of x the two graphs cross, i.e., f(x) = g(x).
% Do this by considering N equally distributed points on the interval, at each
%  point checking whether |f(x) - g(x)| < ε, where ε is some small number.
% Let N and ε be user input to the program and let the result be printed
%  to screen. Run your program with N = 400 and ε = 0.01.
% Explain the output from the program. Finally, try also other values of N,
%  keeping the value of ε fixed. Explain your observations.

function crossing_points = crossing_2_graphs(N, epsilon)
  xs = linspace(-4, 4, N);
  crossing_points = xs(abs(xs - xs.^2 ) < epsilon);
endfunction
