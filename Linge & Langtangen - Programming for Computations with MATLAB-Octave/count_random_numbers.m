% Exercise 2.13 (modified): Frequency of random numbers
% Write a program that takes a positive integer N as input and draws
% N random integers in a configurable interval [n_min, n_max].
%
% The program must:
% 1. Count the occurrences of each integer within the interval.
% 2. Calculate and display the frequency (fraction) for each integer.
% 3. Calculate the mean and standard deviation of these frequencies to
%    evaluate the uniformity of the distribution.
% 4. Display a histogram of the generated numbers.



n_min = 1; n_max = 6;
msg = "Enter how many numbers you want to generate: ";

N = str2num(input(msg, 's'));
nums = randi([n_min, n_max], 1, N);    % array that stores the random numbers
counts = zeros(1, n_max - n_min + 1);  % array that stores each number's count

for k = 1:N
  j = nums(k) - n_min + 1;
  counts(j) = counts(j) + 1;
endfor
frequencies = counts / N;
#sort(nums) % for DEBUG

disp("")
for k = 1:(n_max - n_min + 1)
  printf("The number %d was generated %d times, frequency = %f\n", ...
          k + n_min - 1, counts(k), frequencies(k))
endfor
disp("")

printf("Mean of frequencies = %f (Should be constant: %f)\n", ...
        mean(frequencies), 1/(n_max - n_min + 1));
printf("Standard deviation = %f (Closer to 0 means fairer die)\n", ...
        std(frequencies));

% Histogram
hist(nums, n_min:n_max);
hold on;

ideal_count = N / (n_max - n_min + 1);

line([n_min-0.5, n_max+0.5], [ideal_count, ideal_count], ...
     "linestyle", "--", "color", "r", "linewidth", 2);

legend("Actual counts", "Theoretical ideal");
hold off;
