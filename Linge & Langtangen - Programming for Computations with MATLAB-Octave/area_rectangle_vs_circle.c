% Consider one circle and one rectangle. The circle has a radius r = 10.6.
% The rectangle has sides a and b, but only a is known from the outset.
% Let a = 1.3 and write a program that uses a while loop to find the largest
%  possible integer b that gives a rectangle area smaller than, but as close
%  as possible to, the area of the circle.
% Run the program and confirm that it gives the right answer (which is b = 272).
% Filename: area_rectangle_vs_circle.m.

function [b, iter] = area_rectangle_vs_circle
  r = 10.6;
  area = pi * r^2;

  a = 1.3;
  b_min = 0; b_max = 1;

  iter = 0;
  while a * b_max < area
    b_min = b_max;
    b_max = b_max * 2;
    iter = iter + 1;
  endwhile
  b = (b_min + b_max) / 2;

  while (b_max - b_min) > 0.5
    iter = iter + 1;
    if a * b > area
      b_max = b;
    else
      b_min = b;
    endif
    b = (b_min + b_max) / 2;
    endwhile
  b = round(b);
endfunction
