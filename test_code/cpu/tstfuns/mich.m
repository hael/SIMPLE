function y = mich(x)
% 
% Michalewicz function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n =2.
% 
n = 2; 
m = 10;
s = 0;
for i = 1:n;
    s = s+sin(x(i))*(sin(i*x(i)^2/pi))^(2*m);
end
y = -s;

r Search domain: 0 ≤ xi ≤ π, i = 1, 2, . . . , n.

r Number of local minima: several local minima.

r The global minima:          at n=2, f(x*) = -1.8013.

                              at n=5, f(x*) = --4.687658.

                              at n=10, f(x*) = -9.66015.