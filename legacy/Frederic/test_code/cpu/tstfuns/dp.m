function y = dp(x)
% 
% Dixon and Price function.
% Matlab Code by A. Hedar (Sep. 29, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 25.
% 
n = 25;
s1 = 0;
for j = 2:n;
    s1 = s1+j*(2*x(j)^2-x(j-1))^2;    
end
y = s1+(x(1)-1)^2;

r Search domain: −10 ≤ xi ≤ 10, i = 1, 2, . . . , n.

r The global minima: f(x*) = 0.
