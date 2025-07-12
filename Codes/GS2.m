% The Y bus matrix is
ybus = [2.6923-13.41151 -1.9231+9.61541 0 0 -0.7692+3.8462i
-1.9231+9.6154i 3.6538-18.1942i -0.9615+4.8077i 0 -0.7692+3.8462i
0 -0.9615+4.8077i 2.2115-11.0027i -0.7692+3.8462i -0.4808+2.4030i
0 0 -0.7692+3.8462i 1.1538-5.6742i -0.3846+1.9231i
-0.7692+3.8462i -0.7692+3.8462i -0.4808+2.4030i -0.3846+1.9231i 2.4030-23.6942i];

% The given parameters and initial conditions are
p = [0; -0.96; -0.35; -0.16; 0.24];
q = [0; 0.62; 0.14; 0.08; -0.35];
mv = [1.05; 1; 1; 1; 1.02];
th = zeros(5, 1); % Initial angles are set to zero
v = mv.* exp(1i .* th);

% The Gauss-Seidel iteration
iteration = 0;
max_iterations = 100;
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