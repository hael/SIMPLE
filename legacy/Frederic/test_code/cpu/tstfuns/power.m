function y = power(x)
% 
% Power Sum function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 4.
% 
n = 4;
b = [8,18,44,114];
s_out = 0;
for k = 1:n;
    s_in = 0;
    for j = 1:n
        s_in = s_in+x(j)^k;
    end
    s_out = s_out+(s_in-b(k))^2;
end
y = s_out;

r Search domain: −4 ≤ xi ≤ 5, i = 1, 2, . . . , n.

r The global minima: x* =  (3,-1,0,1, …, 3,-1,0,1), f(x*) = 0.