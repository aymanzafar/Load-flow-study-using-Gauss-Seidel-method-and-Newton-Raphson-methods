function y_new_ndd = newton_divided_difference_interpolation(x, y, x_new)
    n = length(x);
    y_new_ndd = zeros(1, length(x_new));

    for k = 1:length(x_new)
        N = length(x) - 1;
        F = zeros(N+1, N+1);
        F(:,1) = y(:);
        for j = 2:N+1
            for i = 1:N+2-j
                F(i,j) = (F(i+1,j-1) - F(i,j-1))/(x(i+j-1) - x(i));
            end
        end
        y_new_ndd(k) = F(1,2);
        for i = 2:N
            temp = F(1,i+1);
            for j = 1:i-1
                temp = temp*(x_new(k) - x(j)) + F(1,j+1);
            end
            y_new_ndd(k) = y_new_ndd(k) + temp;
        end
    end
end

function y_new_lagrange = lagrange_interpolation(x, y, x_new)
    n = length(x);
    y_new_lagrange = zeros(1, length(x_new));

    for k = 1:length(x_new)
        L = ones(1,n);
        for i = 1:n
            for j = 1:n
                if j ~= i
                    L(i) = L(i) * (x_new(k) - x(j)) / (x(i) - x(j));
                end
            end
        end
        y_new_lagrange(k) = sum(y .* L);
    end
end

% Example dataset from Exercise 1
x = [1, 2, 3, 4, 5];
y = [2, 4, 6, 8, 10];
x_new = [2.5, 3.5, 4.5];

% Newton's divided difference interpolation
y_new_ndd = newton_divided_difference_interpolation(x, y, x_new)

% Lagrange interpolation
y_new_lagrange = lagrange_interpolation(x, y, x_new)

% Verify agreement
if isequal(y_new_ndd, y_new_lagrange)
    fprintf('Agreed\n');
else
    fprintf('Disagreed\n');
end
