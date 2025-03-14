% Given parameters
Is = 0.01e-12;  % Forward bias saturation current (0.01 pA)
Ib = 0.1e-12;   % Breakdown saturation current (0.1 pA)
Vb = 1.3;       % Breakdown voltage (V)
Gp = 0.1;       % Parasitic parallel conductance (Ω- ¹)

% Voltage vector from -1.95V to 0.7V with 200 steps
V = linspace(-1.95, 0.7, 200);

% Compute ideal diode current
I_ideally = Is .* (exp(1.2 .* V / 0.025) - 1) ... % Ideal diode term  
                + Gp .* V ... % Parallel resistor term  
                - Ib .* (exp(1.2 .* (- (V + Vb)) / 0.025) - 1); % Breakdown term
                
% Add 20% random noise to simulate experimental variation
noise = I_ideally .* (0.2 * randn(size(I_ideally)));
I_noisyy = I_ideally + noise;

% Plot data
figure;
subplot(2,1,1)
plot(V, I_ideally, 'b', 'LineWidth', 2); hold on;
plot(V, I_noisyy, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
xlabel('Voltage (V)');
ylabel('Current (A)');
title('Diode I-V Characteristics');
legend('Ideal Data', 'Noisy Data');
grid on;

% Logarithmic plot
subplot(2,1,2)
semilogy(V, abs(I_ideally), 'b', 'LineWidth', 2); hold on;
semilogy(V, abs(I_noisyy), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
xlabel('Voltage (V)');
ylabel('Current (A) (Log Scale)');
title('Diode I-V Characteristics (Log Scale)');
legend('Ideal Data', 'Noisy Data');
grid on;

% Compute polynomial fits (4th and 8th order)
p4 = polyfit(V, I_noisyy, 4);   % 4th order polynomial
p8 = polyfit(V, I_noisyy, 8);   % 8th order polynomial

% Evaluate polynomials at V points
I_poly4 = polyval(p4, V);
I_poly8 = polyval(p8, V);

% Plot original data and polynomial fits
figure;
subplot(2,1,1);
plot(V, I_noisyy, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
plot(V, I_poly4, 'b', 'LineWidth', 2);
plot(V, I_poly8, 'g', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A)');
title('Polynomial Fitting to Noisy Data');
legend('Noisy Data', '4th Order Fit', '8th Order Fit');
grid on;

% Logarithmic plot
subplot(2,1,2);
semilogy(V, abs(I_noisyy), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
semilogy(V, abs(I_poly4), 'b', 'LineWidth', 2);
semilogy(V, abs(I_poly8), 'g', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A) (Log Scale)');
title('Polynomial Fitting (Log Scale)');
legend('Noisy Data', '4th Order Fit', '8th Order Fit');
grid on;

% Given parameters (used for data generation)
Is = 0.01e-12;  % Forward bias saturation current (0.01 pA)
Ib = 0.1e-12;   % Breakdown saturation current (0.1 pA)
Vb = 1.3;       % Breakdown voltage (V)
Gp = 0.1;       % Parasitic parallel conductance (Ω⁻¹)

% Voltage vector from -1.95V to 0.7V with 200 steps
V = linspace(-1.95, 0.7, 200);

% Compute ideal diode current
I_ideally = Is .* (exp(1.2 .* V / 0.025) - 1) ... % Ideal diode term
         + Gp .* V ... % Parallel resistor term
         - Ib .* (exp(1.2 .* (- (V + Vb)) / 0.025) - 1); % Breakdown term

% Add 20% random noise to simulate experimental variation
noise = I_ideally .* (0.2 * randn(size(I_ideally))); 
I_noisyy = I_ideally + noise;

% Define the custom nonlinear model
model = fittype('A*(exp(1.2*x/0.025)-1) + B*x - C*(exp(1.2*(-(x+D))/0.025)-1)', ...
                'independent', 'x', 'coefficients', {'A', 'B', 'C', 'D'});

% Initial guesses for the parameters
initial_guesses = [1e-12, 0.01, 1e-12, 1];

% Fit the model to the noisy data
[fit_result, gof] = fit(V.', I_noisyy.', model, 'StartPoint', initial_guesses);

% Extract the fitted parameters
A_fit = fit_result.A;
B_fit = fit_result.B;
C_fit = fit_result.C;
D_fit = fit_result.D;

% Generate fitted current values
I_fit = A_fit * (exp(1.2 * V / 0.025) - 1) + B_fit * V - C_fit * (exp(1.2 * (-(V + D_fit)) / 0.025) - 1);

% Given known values for B and D (from Equation 1)
B_known = 0.01;  % B (Gp) = 0.01 from the original data generation
D_known = 1.3;   % D (Vb) = 1.3 from the original data generation

% Define the custom nonlinear model with B and D fixed
model_fixed_B_D = fittype('A*(exp(1.2*x/0.025)-1) + 0.01*x - C*(exp(1.2*(-(x+1.3))/0.025)-1)', ...
                          'independent', 'x', 'coefficients', {'A', 'C'});

% Initial guesses for parameters A and C
initial_guesses_fixed_B_D = [1e-12, 1e-12];  % Initial guesses for A and C

% Fit the model to the noisy data with B and D fixed
[fit_result_fixed_B_D, gof_fixed_B_D] = fit(V.', I_noisyy.', model_fixed_B_D, 'StartPoint', initial_guesses_fixed_B_D);

% Extract the fitted parameters A and C
A_fit_fixed_B_D = fit_result_fixed_B_D.A;
C_fit_fixed_B_D = fit_result_fixed_B_D.C;

% Generate fitted current values with B and D fixed
I_fit_fixed_B_D = A_fit_fixed_B_D * (exp(1.2 * V / 0.025) - 1) + B_known * V - C_fit_fixed_B_D * (exp(1.2 * (-(V + D_known)) / 0.025) - 1);

D = 1.3;  % Known value of D

% Define the custom nonlinear model with fixed D
model_fixed_D = fittype('A*(exp(1.2*x/0.025)-1) + B*x - C*(exp(1.2*(-(x+1.3))/0.025)-1)', ...
                        'independent', 'x', 'coefficients', {'A', 'B', 'C'});

% Initial guesses for parameters A, B, and C
initial_guesses_fixed_D = [1e-12, 0.01, 1e-12];  % Initial guesses for A, B, and C

% Fit the model to the noisy data with D fixed
[fit_result_fixed_D, gof_fixed_D] = fit(V.', I_noisyy.', model_fixed_D, 'StartPoint', initial_guesses_fixed_D);

% Extract the fitted parameters
A_fit_fixed_D = fit_result_fixed_D.A;
B_fit_fixed_D = fit_result_fixed_D.B;
C_fit_fixed_D = fit_result_fixed_D.C;

% Generate fitted current values with fixed D
I_fit_fixed_D = A_fit_fixed_D * (exp(1.2 * V / 0.025) - 1) + B_fit_fixed_D * V - C_fit_fixed_D * (exp(1.2 * (-(V + D)) / 0.025) - 1);

% Plot the data and the fitted curve with D fixed
figure;
subplot(2,1,1);
plot(V, I_noisyy, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
plot(V, I_fit_fixed_D, 'b', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A)');
title('Nonlinear Curve Fitting with Fixed D to Noisy Data');
legend('Noisy Data', 'Fitted Curve (Fixed D)');
grid on;

% Logarithmic plot
subplot(2,1,2);
semilogy(V, abs(I_noisyy), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
semilogy(V, abs(I_fit_fixed_D), 'b', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A) (Log Scale)');
title('Nonlinear Curve Fitting with Fixed D (Log Scale)');
legend('Noisy Data', 'Fitted Curve (Fixed D)');
grid on;

% Plot the data and the fitted curve with B and D fixed
figure;
subplot(2,1,1);
plot(V, I_noisyy, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
plot(V, I_fit_fixed_B_D, 'b', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A)');
title('Nonlinear Curve Fitting with Fixed B and D to Noisy Data');
legend('Noisy Data', 'Fitted Curve (Fixed B, D)');
grid on;

% Logarithmic plot
subplot(2,1,2);
semilogy(V, abs(I_noisyy), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
semilogy(V, abs(I_fit_fixed_B_D), 'b', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A) (Log Scale)');
title('Nonlinear Curve Fitting with Fixed B and D (Log Scale)');
legend('Noisy Data', 'Fitted Curve (Fixed B, D)');
grid on;

% Display the fitted parameters A and C
disp('Fitted Parameters (with B and D fixed):');
fprintf('A (Is)  = %.3e A\n', A_fit_fixed_B_D);
fprintf('C (Ib)  = %.3e A\n', C_fit_fixed_B_D);


% Plot the data and the fitted curve
figure;
subplot(2,1,1);
plot(V, I_noisyy, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
plot(V, I_fit, 'b', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A)');
title('Nonlinear Curve Fitting to Noisy Data');
legend('Noisy Data', 'Fitted Curve');
grid on;

% Logarithmic plot
subplot(2,1,2);
semilogy(V, abs(I_noisyy), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
semilogy(V, abs(I_fit), 'b', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A) (Log Scale)');
title('Nonlinear Curve Fitting (Log Scale)');
legend('Noisy Data', 'Fitted Curve');
grid on;

% Display the fitted parameters
disp('Fitted Parameters:');
fprintf('A (Is)  = %.3e A\n', A_fit);
fprintf('B (Gp)  = %.3e Ω⁻¹\n', B_fit);
fprintf('C (Ib)  = %.3e A\n', C_fit);
fprintf('D (Vb)  = %.3f V\n', D_fit);

% Define inputs (voltage) and targets (current)
inputs = V.';    % Voltage as column vector
targets = I_noisyy.';  % Noisy current data as column vector

% Define neural network structure
hiddenLayerSize = 10;   % 10 neurons in hidden layer
net = fitnet(hiddenLayerSize);

% Configure training parameters
net.divideParam.trainRatio = 70/100;  % 70% training data
net.divideParam.valRatio = 15/100;    % 15% validation data
net.divideParam.testRatio = 15/100;   % 15% test data

% Train the neural network
[net, tr] = train(net, inputs, targets);

% Use trained network to generate output predictions
outputs = net(inputs);

% Compute errors and performance
errors = gsubtract(outputs, targets);
performance = perform(net, targets, outputs);

% Display neural network model
view(net);

% Plot the results
figure;
subplot(2,1,1);
plot(V, I_noisyy, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
plot(V, outputs, 'b', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A)');
title('Neural Network Fitting for Diode Data');
legend('Noisy Data', 'Neural Net Fit');
grid on;

% Logarithmic plot
subplot(2,1,2);
semilogy(V, abs(I_noisyy), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
semilogy(V, abs(outputs), 'b', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A) (Log Scale)');
title('Neural Network Fit (Log Scale)');
legend('Noisy Data', 'Neural Net Fit');
grid on;

% Display performance
fprintf('Neural Network Performance: %.6f\n', performance);

% Define inputs (voltage) and targets (current)
inputs = V.';      % Voltage as column vector
targets = I_noisyy.';  % Noisy current data as column vector

% Train GPR model
gprModel = fitrgp(inputs, targets);

% Predict current values using the trained GPR model
I_gpr = predict(gprModel, inputs);

% Plot original data and GPR fit
figure;
subplot(2,1,1);
plot(V, I_noisyy, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
plot(V, I_gpr, 'b', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A)');
title('Gaussian Process Regression (GPR) Fit for Diode Data');
legend('Noisy Data', 'GPR Fit');
grid on;

% Logarithmic plot
subplot(2,1,2);
semilogy(V, abs(I_noisyy), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); hold on;
semilogy(V, abs(I_gpr), 'b', 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A) (Log Scale)');
title('Gaussian Process Regression Fit (Log Scale)');
legend('Noisy Data', 'GPR Fit');
grid on;

% Display GPR model details
disp('Trained Gaussian Process Regression Model:');
disp(gprModel);
