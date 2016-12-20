function y = schw(x)
% 
% Schwefel function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 2;
s = sum(-x.*sin(sqrt(abs(x))));
y = 418.9829*n+s;

r Search domain: −500 ≤ xi ≤ 500, i = 1, 2, . . . , n.

r Number of local minima: several local minima.

r The global minima: x* =  (1, …, 1), f(x*) = 0.