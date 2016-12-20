function y = zakh(x)
% 
% Zakharov function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = 2;
s1 = 0;
s2 = 0;
for j = 1:n;
    s1 = s1+x(j)^2;
    s2 = s2+0.5*j*x(j);
end
y = s1+s2^2+s2^4;

r Search domain: −5 ≤ xi ≤ 10, i = 1, 2, . . . , n.

r Number of local minima: no local minimum except the global one.

r The global minima: x* =  (0, …, 0), Zn(x*) = 0.
