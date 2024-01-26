%{
{%
Math 3800 - A1
Nick Cooley
101021174
Date: 14 Jan 2024
%}

%% Q1

% See write-up


%% Q2 (4 points)

%{
find a c such that x = c where c is the induced L2 norm of A that is, 
c = square root of maximum eigen value of A'*A and c = x in A.

We have x^2 +24x + 144 + L^2 -26 L -Lx^2 = 0, where L is lamda.
But we require x = Lmax^(1/2) by defintiion of L2 norm.

x^2 +24x + 144 + x^4 -26x^2 - x^4 = 0

-25x^2 +24x + 144 = 0
%}

clc; clearvars; format long;

x1 = (-24 + sqrt(24^2 + 4 * 25 * 144)) / (2 * -25);
x2 = (-24 - sqrt(24^2 + 4 * 25 * 144)) / (2 * -25); % suspect this as its larger.

fprintf("We have found two roots x1:%f, and x2:%f.\n", x1, x2);

A = @(x) [-1 4; 3 x];
ATA = @(x) [10 -4+3*x; -4+3*x 16+x^2];

fprintf("Unfortionately with x1 = %f we have sqrt(min(eig(ATA(x1)))) = %f.\n", x1, sqrt(min(eig(ATA(x1)))) );
fprintf("and x2 = %f we have sqrt(min(eig(ATA(x2)))) = %f.\n", x2, sqrt(min(eig(ATA(x2)))) );

fprintf("\nWhat question 2 is really asking follows:\n")

%{
No solution exists where x is equal to the L2 norm of A(x)

After discussing with the Ms Nookala and Professor Cheung, it seems the above is not
desired. 
%}

tol = 1e-17;
arrn0 = []; arr0 = [];
diff = @(x) norm(A(x), 2) - sqrt(max(eig(ATA(x))));

fprintf("\nTolerance is %e.\n", tol)
for x = 1 : 1 : 100
    delta = diff(x);
    if abs(delta) < tol % I do not trust logical == 0
        arr0 = [arr0 ; x delta]; % no flag exists to suppress this warning.   equal
    else
        arrn0 = [arrn0 ; x delta]; % no flag exists to suppress this warning. not equal
        fprintf("x = %02d is such a value.\n", x); % we find lots.
    end
end


fprintf("\nThe values x that satisfy norm(A(x), 2) - sqrt(max(eig(ATA(x))) < %e are:\n", tol);
for x = 1 : length(arr0)
    fprintf("%03d with a difference of %.18e\n", arr0(x, 1), arr0(x, 2)  )
end

fprintf("\nThe values x that are outside the open interval %e are:\n", tol)
for x = 1 : length(arrn0)
   fprintf("%03d with a difference of %.18e\n", arrn0(x, 1), arrn0(x, 2) )
end


%% Q3

clc; clf; hold off; clearvars; format long;

disp('a. See Figure 1.')

xmin = 0.; xmax = 0.5235;
x = linspace(xmin, xmax, 244);

cos_ = @(x) 1 - x.^2/2;

plot(x, cos(x), 'Color','red'); hold on; plot(x, cos_(x), 'Color','blue'); hold off;

disp('b.')

cosine_ = arrayfun(cos_, x);
diff = abs(cos(x) - cosine_);
max_delta_x = max(diff);

fprintf("The maximum difference between the real vector and approximated \n taylor series vector to 2 terms is: %.15f.\n", max_delta_x)

disp('c.')

x = 0.2453; 
y_ = cos_(x);
y = cos(x);
err = y_ - y;

er = sprintf("The forward error is: %.8f", err);
disp(er); % throws warning to use fprintf, we are told to use sprintf.

er = sprintf("The relative forward error is: %.8f", err/y);
disp(er);

x_ = acos(y_);

er = sprintf("\nThe backwards error is: %.8f", x_ - x);
disp(er);

er = sprintf("The relative backwards error is: %.8f", (x_ - x) / x);
disp(er);


%% Q4

% See write-up

