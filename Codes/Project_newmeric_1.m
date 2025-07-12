clc; close all; clear all;
disp('Which mathode you want');
disp('For Gauss-Seidel Prees 1');
disp('For raphson method solution Prees 2');
z = input('Input Only 1 or 2...');
ybus = input('Input ybus in Matrix... ');
p = input('Input p in Matrix... ');
q = input('Input q in Matrix... ');
mv = input('Input mv in Matrix... ');
th = zeros(5, 1); % Initial angles are set to zero
v = mv.* exp(1i .* th);

             % The Gauss-Seidel iteration

if z == 1 
iteration = 0;
max_iterations = 1000;
del = 1; % Initialize a value larger than the tolerance
while del > 1e-4 && iteration < max_iterations
    iteration = iteration + 1;
    %P-Q buses
    for i = 2:4
        tmpl = (p(i) - 1i * q(i)) / conj(v(i));
        tmp2 = 0;
        for k = 1:5
            if i == k
                tmp2 = tmp2 + 0;
            else
                tmp2 = tmp2 + ybus(i, k) *v(k);
            end
        end
        v(i) = (tmpl - tmp2) / ybus(i, i);
    end
    %P-V bus (Bus 5)
    q5 = 0;
    for i = 1:5
        q5 = q5 + ybus(5, i) *v(i);
    end
    q5 = -imag(conj(v(5)) * q5);
    tmpl = (p(5) - 1i * q5) / conj(v(5));
    tmp2 = 0;
    for k = 1:4
        tmp2 = tmp2 + ybus(5, k) * v(k);
    end
    vt  = (tmpl  - tmp2) / ybus(5, 5);
    v (5) =  abs(vt) * exp(1i * angle(v(5)));
    % Calculate S for all buses
    S = zeros (5, 1);
    for i = 1:5
        sm = 0;
        for k = 1:5
            sm = sm + ybus(i, k) * v(k);
        end
        S(i) = conj(v(i)) * sm;
    end
    % Calculate the mismatch
    delp = p - real(S);
    delq = q + imag(S);
    delpq = [delp(2:5); delq(2:5)];
    del = max(abs(delpq));
    fprintf('Iteration: %d, Max Mismatch: %.6f\n', iteration, del);
end

% Display results
disp('Bus Voltages and Angles:')
for i = 1:5
    fprintf('Bus %d Voltage Magnitude: %.4f, Voltage Angle: %.4f degrees\n', ...
    i, abs(v(i)), rad2deg(angle(v(i))));
end
end

                   % The Newton-Raphson iteration

if z == 2
iteration = 0;
max_iterations = 1000;
tolerance = 1e-4;
delta = inf;

while delta > tolerance && iteration < max_iterations
    iteration = iteration + 1;
    
    % Construct Jacobian matrix
    J = zeros(10, 10); % Jacobian matrix size (5 buses, 2 unknowns per bus)
    F = zeros(10, 1); % Function evaluations
    
    for i = 1:5
        for j = 1:5
            if i ~= j
                % Off-diagonal elements
                J(i, j) = -imag(conj(v(i)) * ybus(i, j));
                J(i+5, j+5) = -real(conj(v(i)) * ybus(i, j));
            else
                % Diagonal elements
                sm = 0;
                for k = 1:5
                    if k ~= i
                        sm = sm + ybus(i, k) * v(k);
                    end
                end
                J(i, j) = real(conj(v(i)) * sm) + q(i);
                J(i+5, j+5) = -imag(conj(v(i)) * sm) - p(i);
            end
        end
        % Function evaluations
        sm = 0;
        for k = 1:5
            if k ~= i
                sm = sm + ybus(i, k) * v(k);
            end
        end
        F(i) = v(i) * conj(sm) - p(i) - 1i * q(i);
        F(i+5) = imag(v(i) * conj(sm)) - q(i) + 1i * p(i);
    end
    
    % Newton-Raphson step
    increment = J \ (-F);
    v = v + increment(1:5);
    th = th + increment(6:10);
    
    % Recalculate maximum mismatch
    S = zeros(5, 1);
    for i = 1:5
        sm = 0;
        for k = 1:5
            sm = sm + ybus(i, k) * v(k);
        end
        S(i) = conj(v(i)) * sm;
    end
    delp = p - real(S);
    delq = q + imag(S);
    delpq = [delp(2:5); delq(2:5)];
    delta = max(abs(delpq));
    
    fprintf('Iteration: %d, Max Mismatch: %.6f\n', iteration, delta);
end

% Display results
disp('Bus Voltages and Angles:')
for i = 1:5
    fprintf('Bus %d Voltage Magnitude: %.4f, Voltage Angle: %.4f degrees\n', ...
        i, abs(v(i)), rad2deg(angle(v(i))));
end
end