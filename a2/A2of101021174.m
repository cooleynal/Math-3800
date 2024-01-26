%{
{%
Math 3800 - A2
Nick Cooley
101021174
Date: 25 Jan 2024
%}

%% Q1

format long;
clc;

A = [0 1 -1 ; 2 -3 5 ; 1 -3 4];
[L, U, P] = lu(A);
L
U
P
LU = L*U
PA = P*A
isequal(PA, LU)



%% Q2

clc;
f = @(x) 3*x - x^(5/3);



%% Q3
clc;
f = @(x) x^3 - 3*x^2 -13*x;
f = @(x) x * (-13 + x * (-3 + x)); 

x_ = @(x) x - f(x)^2/(f(x +f(x) - f(x)))

%% Q4
clc;

f = @(x,y) [x^3 + 2^2 - 7 ; 3/(1 + x -2*y) - 4];