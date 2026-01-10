% Exercise 2.10: Sort array with numbers
% The built-in function rand may be used to draw pseudo-random numbers for the
%  standard uniform distribution between 0 and 1 (exclusive at both ends).
%  See help rand.
% Write a script that generates an array of 6 random numbers between 0 and 10.
% The program should then sort the array so that numbers appear in increasing
%  order.
% Let the program make a formatted print of the array to the screen both before
%  and after sorting. The printouts should appear on the screen so that
%  comparison is made easy. Confirm that the array has been sorted correctly.
% Filename: sort_numbers.m



nums = floor(rand(1, 6) * 11);
nums_sorted = sort(nums);
nums_len = length(nums);

printf("Before sorting:\n")
disp(nums)

% Insertion Sort
for k = 2:(nums_len)
  val = nums(k);
	j = k - 1;
  while (j >= 1) && (val < nums(j))
    nums(j + 1) = nums(j);
    j = j - 1;
  endwhile
	nums(j+1) = val;
endfor

printf("After sorting:\n")
disp(nums)
printf("Using built-in function sort():\n")
disp(nums_sorted)
