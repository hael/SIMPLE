function y = rast(x)
% 
% Rastrigin function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 2; 
s = 0;
for j = 1:n
    s = s+(x(j)^2-10*cos(2*pi*x(j))); 
end
y = 10*n+s;

r Search domain: −5.12 ≤ xi ≤ 5.12, i = 1, 2, . . . , n.

r Number of local minima: several local minima.

r The global minima: x* =  (0, …, 0), f(x*) = 0.

r MATLAB Code: download rast.m