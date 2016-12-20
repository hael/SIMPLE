function  y = trid(x)
% 
% Trid function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 10.
% 
n = 10;
s1 = 0;
s2 = 0;
for j = 1:n;
    s1 = s1+(x(j)-1)^2;    
end
for j = 2:n;
    s2 = s2+x(j)*x(j-1);    
end
y = s1-s2;

r Search domain: −n2 ≤ xi ≤ n2, i = 1, 2, . . . , n.

r Number of local minima: no local minimum except the global one.

r The global minima:          f(x*) = -50    for n=6,

                              f(x*) = -200 for n=10.