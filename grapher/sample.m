clear all
close all

epsilon = 10^(-5);
c1 = -1;
c2 = 0;
range = [-1 1];

g = @(x) airy(x/(epsilon)^(1/3));
h = @(x) airy(2,x/(epsilon)^(1/3));

f = @(x) c1*g(x) + c2*h(x);

fplot(f, range)
