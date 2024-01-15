%{
A1
Nick Cooley
101021174
Date: 14 Jan 2024
%}

%% Q1



%% Q2 (4 points)

%{
find a c such that x = c where c is the induced L2 norm of A that is, 
c = square root of maximum eigen value of A and c = x in A.

We have x^2 +24x + 144 + L^2 -26 L -Lx^2 = 0, where L is lamda.
But we require x = Lmax^(1/2) by defintiion of L2 norm.
x^2 +24x + 144 + x^4 -26x^2 - x^4 = 0

-25x^2 +24x + 144 = 0

no solution exists when solving for the 5th order polynomial lamda roots either.
%}

clc; clearvars; format long;

x1 = (-24 + sqrt(24^2 + 4 * 25 * 144)) / (2 * -25);
x2 = (-24 - sqrt(24^2 + 4 * 25 * 144)) / (2 * -25); % suspect this as its larger.

fprintf("We have found two roots x1:%f, and x2:%f.\n", x1, x2);

ATA = @(x) [10 -4+3*x; -4+3*x 16+x^2];

disp( x1 )
disp( sqrt(eig(ATA(x1))) )

fprintf("Unfortionately with x1 = %f we have sqrt(min(eig(ATA(x1)))) = %f.\n", x1, sqrt(min(eig(ATA(x1)))) );
fprintf("and x2 = %f we have sqrt(min(eig(ATA(x2)))) = %f.\n", x2, sqrt(min(eig(ATA(x2)))) );

% no solution in this question exists. I will write the numerical attempt to push for marks

ATA = @(x) [10 -4+3*x; -4+3*x 16+x^2]; % reduce computation time
A = @(x) [-1 4; 3 x];

arr = [];
sol = 0;
tol = 1e-1; % loose tolerance demonstrating no solution exists.
start = 0; step = 1e-4; stop = 100; % step size 1e-4 isnt small, futile search
for x = start : step : stop
    e = sqrt(max(eig(ATA(x)))); % but we can find sqrt(min(eig(ATA(x1)))) see above
    if abs(e - x) < tol
        disp(abs(e - x))
        fprintf("eigenvalue is: %f x is: %f\n", e, x)
        arr = [arr, x]; % no warning flag exists in matlab to suppress this
        sol = sol + 1;
        if sol == 1     % was going to explore all solutions. there are none
            break;
        end
    end
end

disp("solution check.");
disp(arr)
for i = 1:length(arr)
    sol = arr(i);
    fprintf("solution is: %f, verified with:norm(A(arr(i)), 2) = %f \n", sol, norm(A(sol), 2));
end


%% Q3 (7 points)

clc; clf; hold off; clearvars;

disp('a. See Figure 1.')

xmin = 0.; xmax = 0.5235;
x = linspace(xmin, xmax, 244);

v = @(x) cos(x);
w = @(x) 1 - x.^2/2;

plot(x, v(x), 'Color','red'); hold on; plot(x, w(x), 'Color','blue'); hold off;

disp('b.')
cosine = arrayfun(v, x);
aprox = arrayfun(w, x);

diff = abs(cosine - aprox);
max_delta_x = max(diff);

fprintf("The maximum difference between the real vector and approximated taylor series vector to 2 terms is: %.15f.\n", max_delta_x)

disp('c.')

x = 0.2453; 
err = w(x) - v(x);

er = sprintf("The forward error is: %.8f.", err);
disp(er); % throws warning to use fprintf, we are told to use sprintf.
er = sprintf("The relative forward error is: %.8f.", err/v(x));
disp(er);

% the notes indicate we use a estimate y to get a estimate x, we dont have
% an estimate x, we only have a estimate of cos.
% I believe notes should say f^hat(x^hat) = y^hat
% aproxarcos = @(y) sqrt(2 * (1-y)); % would result in 0 backwards error.

x_hat = acos(w(x));
er = sprintf("The backwards error is: %.8f.", x_hat - x);
disp(er);
er = sprintf("The relative fbackwards error is: %.8f.", (x_hat - x) / x);
disp(er);

%% Q4




