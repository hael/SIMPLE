function y = griewank(x)
% 
% Griewank function
% Matlab Code by A. Hedar (Sep. 29, 2005).
% The number of variables n should be adjusted below.
% The default value of n =2.
% 
n = 2;
fr = 4000;
s = 0;
p = 1;
for j = 1:n; s = s+x(j)^2; end
for j = 1:n; p = p*cos(x(j)/sqrt(j)); end
y = s/fr-p+1;

r Search domain: −600 ≤ xi ≤ 600, i = 1, 2, . . . , n.

r Number of local minima: several local minima.

r The global minima: x* =  (0, …, 0), f(x*) = 0.