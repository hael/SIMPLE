function y = sum2(x)
% 
% Sum Squares function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 20.
% 
n = 20;
s = 0;
for j = 1:n  
    s=s+j*x(j)^2; 
end
y = s;

r Search domain: −10 ≤ xi ≤ 10, i = 1, 2, . . . , n.

r Number of local minima: no local minimum except the global one.

r The global minima: x* =  (0, …, 0), f(x*) = 0.