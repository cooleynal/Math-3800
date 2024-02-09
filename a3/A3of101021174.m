%% q1a
clc; clearvars, close;
% format long;
format rat;
% syms x % not allowed

X = [0, 2, 4, 6];
Y = [0, 7 0, -2];

len = length(X);
for i = 1: len
    denom = 1;
    str = '';
    for j = 1: len
        if j == i
            continue;
        end 
    str = [str, '(x - ', num2str(X(j)), ')'];   % syms x % not allowed
    denom = denom * (X(i) - X(j));
    end
    fprintf('term: %d\n', i);
    coe = Y(i)/denom
    str                                         % syms x % not allowed
    
end

% p = @(x) (19/48) * x.^3 + (-33/8) * x.^2 + (61/6) * x + 0;
p = @(x) (-1/24) * (x - 0) .* (x - 2) .* (x - 4) + (7/16) * (x - 0) .* (x - 4) .* (x - 6);

hold on;
for i = 1: len
    plot(X(i), Y(i), 'marker', 'o', 'MarkerEdgeColor', 'green');
end

xs = linspace(-0.5, 8, 1e2);
plot(xs, p(xs), '-b');

yline(0,'k--');
xline(0,'k--');
axis([-1, 9, -4, +10])
hold off;


%% q2 b
clc; clearvars;
format rat;

X = [0, 2, 4, 6];
Y = [0, 7 0, -2];

len = length(X);

M = zeros(len, len + 1);
M(:, 1) = X'; % could have initalized as column vector, this is cooler.
M(:, 2) = Y';

for i = 1: len
    for j = 1: len
        if j + i > len
            continue;
        end 
        M(i+j, i+2) = (M(i+j, i+1) - M(i+j-1, i+1)) / (M(i+j, 1) - M(j, 1) ) ;

    end
end

M

for i = 1:len
    str = num2str(rat(M(i, i+1))); % format rat, and rat are not even the same... MATLAB...
    for j = 1:i-1
        str = [str, '(x - ', num2str(X(j)), ')'];  
    end
    str
end


% p = @(x) (19/48) * x.^3 + (-33/8) * x.^2 + (61/6) * x + 0;
p = @(x) (19/48) * (x - 0) .* (x - 2) .* (x - 4) + (-7/4) * (x - 0) .* (x - 2) + (7/2)  * (x - 0);


hold on;
for i = 1: len
    plot(X(i), Y(i), 'marker', 'o', 'MarkerEdgeColor', 'green');
end

xs = linspace(-0.5, 8, 1e2);
plot(xs, p(xs), '-b');

yline(0,'k--');
xline(0,'k--');
axis([-1, 9, -4, +10])
hold off;



%% q3
clc; clearvars;
format long;

%% q4
clc; clearvars;
format long;


%% q5
clc; clearvars;
format long;


%% q6
clc; clearvars;
format long;