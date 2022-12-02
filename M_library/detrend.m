function X = detrend(Y)

[T, M] = size(Y);
x = 1:T;
x = [ones(T, 1), x'];
xx = x'*x;
if M == 1
    xy = x'*Y;
    b = xx\xy;
    X = Y - x*b;
else
    X = Y;
    for i = 1:M
        xy = x'*Y(:, i);
        b = xx\xy;
        X(:, i) = Y(:, i) - x*b;
    end
    
end
end