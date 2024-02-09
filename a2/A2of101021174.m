%{
Math 3800 - A2
Nick Cooley
101021174
Date: 29 Jan 2024
%}

%% Q1

clearvars; clc; format long;

A = [0 1 -1 ; 2 -3 5 ; 1 -3 4];
[L, U, P] = lu(A);
L
U
P
LU = L*U
PA = P*A
isequal(PA, LU)
rref(LU)


%% Q2

% place cobwebmodded.m in the same directory as this.
clearvars; clc; format long; clf;

% f = @(x) x * ( 3 - x^(2/3) );
f = @(x) 3*x - x^(5/3);
cobwebmodded(f, [0 6], 2, 1, 8);



%% Q3
clearvars; clc; format long; clf;
% polynomial of degree 3 has 1, 2,or  3 roots.

f = @(x) x .* (-13 + x .* (-3 + x)); 
x_ = @(x) x - f(x)^2/( f(x +f(x)) - f(x) );
df = @(x) -13 + x*(-6 + 3*x);

ca = 3;
cb = -6;
cc = -13;

cx1 = (- cb - (cb^2 - 4*ca*cc)^(1/2) )/(2*ca);
cx2 = (- cb + (cb^2 - 4*ca*cc)^(1/2) ) /(2*ca);

xs = linspace(-6, 8, 1e2);
ys = f(xs);
x = 6.; % implicit float

slots = 2024; % i=0 to i=2024  
for k = 1 : slots % -1 so we dont step 2025 times
    x = x_(x);
end
x
f(x)
% k
% x
% cx1
% cx2
% x = (3 + 61^(1/2) / 2)
% round(x, 8)
hold on;
% crits
plot(cx1, f(cx1), 'marker', 'o', 'MarkerEdgeColor', 'green');
plot(cx2, f(cx2), 'marker', 'o', 'MarkerEdgeColor', 'green');
% roots
plot(x, f(x), 'marker', 'o', 'MarkerEdgeColor', 'red');
plot(0, f(0), 'marker', 'o', 'MarkerEdgeColor', 'red');
plot( (3- 61^(1/2))/2, f( (3- 61^(1/2))/2 ),'marker', 'o', 'MarkerEdgeColor', 'red');
plot(xs, ys, 'color', 'blue');
% extra
yline(0,'k--');
xline(0,'k--');
axis([-6, 8, -44, +15])
hold off;


%% Q4
clearvars; clc; format long; 

f = @(v) [ v(1)^3+v(2)^2-7; 3/(1 + v(1) - 2*v(2))-4 ];
Jf = @(v) [ 3*v(1)^2 2*v(2) ; -3*(1+v(1)-2*v(2))^(-2) 6*(1+v(1)-2*v(2))^(-2) ];

a = 2;
b = 1;
v0 = [a; b];
% iteration_count = 1 + 51; 
iteration_count = 1 + 3;
arr = zeros(2, 1, iteration_count);
arr(:, :, 1) = v0;

for k = 1 : iteration_count - 1 % last indix is set.
    arr(:, :, k+1) = arr(:, :, k) - Jf(arr(:, :, k)) \f(arr(:, :, k));
end

% arr(:, :, 52) 
arr(:, :, 1)
arr(:, :, 2)
arr(:, :, 3)
arr(:, :, 4)

f(arr(:, :, 1))
f(arr(:, :, 4))
norm(f(arr(:, :, 1)) - [0; 0], 2)
norm(f(arr(:, :, 4)) - [0; 0], 2)

