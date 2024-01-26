
clc;
x = [1 2 3];
g1 = @(x) (x+1)/x.^2;
g2 = @(x) (x+1)./x.^2;
display(g1(x))
display(g2(x))