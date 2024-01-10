% MATLAB Task.
% Nick Cooley.
% 101021174.
% Date: 09 Jan 2024.

%% Task 1

clc;

format long;
target = 3800;   % goal
step = 0.001;    % precision, I assume the domain is in the reals
m = 2000;        % initialize arbitrary less than expr(m) == goal

expr = @(m) ceil(m * sqrt( abs(sin(4) ) * log(20) ));

while abs(target - expr(m)) > step
    m = m + step;
end

% not displayed
% 3800 == expr(m);
% expr(m);

disp(m)
disp(expr(m));
disp(string(3800 == expr(m)));
fprintf("The starting value of m : %f satisfies the condition 3800 == expr(%f) which evaluates to: %s.\n", m, m, string(3800 == expr(m)));

% There exists a range of solutions that could be provided. Question only
% asks for a single solution, so im afraid to oversolve the problem and
% possible loose marks. Here we take the smaller value in the solution set.


%% Task 2

xmin = -5; xmax = 5;

f = @(x) 6*x.^3 - 2*x;
g = @(x) (8*x + 6) / (4*x.^2 - x + 3);
% g = @(x) (4*x.^2 - x + 3) \ (8*x + 6); % this works to

x = linspace(xmin, xmax, 111);
v = arrayfun(f, x);
w = arrayfun(g, x);

plot(x, v); hold on; plot(x, w);
% plot(x, v, 'b', x, w, 'r'); % this works to

clc; hold off;

xline(0); yline(0);
xlim([xmin, xmax]); ylim([xmin, xmax]);


%% Task 3

% To find a single critical point of p, we will take the derivative to attain the slope of p over all x.
% dp = (3x^2 +16x +b) dx
% We now will look for all zeros.
% zeros = ( -16 pm sqrt(16^2 - 4 * 3 * b) ) / (2 * 6)
% in particular: a single zero in dp requires: 16^2 -4 * 3 * b = 0
% 16^2 = 12*b, and b = 16^2/12.
% with dp, we find the zero is x = -16/6.

b = 16^2/12;
p = @(x) x.^3 + 8*x.^2 + b*x + 3;

clc; clf; hold off;
fplot(p); hold on;

critx = -16/6;
plot(critx, p(critx), 'O')

xline(0); yline(0);
xmin = -6; xmax = 1;
ymin = -25; ymax = 1;
xlim([xmin, xmax]); ylim([ymin, ymax]);
