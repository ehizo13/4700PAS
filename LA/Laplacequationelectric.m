% Define grid size
nx = 50; % Number of points in x-direction
ny = 50; % Number of points in y-direction
V = zeros(nx, ny); % Initialize potential matrix

% Iteration parameters
maxIter = 5000;    % Maximum number of iterations
tolerance = 1e-6;  % Convergence criterion
error = inf;       % Initialize error
iter = 0;          % Iteration counter

% Set up figure for visualization
figure;
colormap(jet);
frames = struct('cdata', [], 'colormap', []); % Preallocate movie frames

%% --- First Phase: Solve Laplace's Equation ---
fprintf('Solving with initial boundary conditions...\n');

% % Apply initial boundary conditions
% V(:,1) = 1;      % Left boundary (V = 1)
% V(:,end) = 0;    % Right boundary (V = 0)
% V(1,:) = V(2,:); % Top boundary (Neumann: ∂V/∂y = 0)
% V(end,:) = V(end-1,:); % Bottom boundary (Neumann: ∂V/∂y = 0)

% Apply initial boundary conditions
V(:,1) = 1;      % Left boundary (V = 1)
V(:,end) = 1;    % Right boundary (V = 0)
V(1,:) = 0; % Top boundary (Neumann: ∂V/∂y = 0)
V(end,:) = 0; % Bottom boundary (Neumann: ∂V/∂y = 0)

while error > tolerance && iter < maxIter
    V_old = V; % Store previous iteration
   
    % Update the solution for interior points
    for i = 2:nx-1
        for j = 2:ny-1
            V(i,j) = 0.25 * (V(i+1,j) + V(i-1,j) + V(i,j+1) + V(i,j-1));
        end
    end

    % Reapply boundary conditions
    V(:,1) = 1;      
    V(:,end) = 1;    
    V(1,:) = 0;
    V(end,:) = 0;

    % Compute error
    error = max(max(abs(V - V_old)));
    iter = iter + 1;
   
    % Capture frames for movie
    if mod(iter, 10) == 0
        frames(iter/10) = getframe(gcf);
    end
   
    % Display progress
    if mod(iter, 100) == 0
        fprintf('Iteration: %d, Max Error: %.6e\n', iter, error);
        surf(V)
        shading interp
        pause(0.1)
    end
end

fprintf('Phase 1 converged after %d iterations with max error %.6e\n', iter, error);

%% --- Compute the Electric Field ---
fprintf('Computing the electric field...\n');

[Ey, Ex] = gradient(-V); % Compute electric field as negative gradient of V

% Generate grid for plotting
[X, Y] = meshgrid(1:ny, 1:nx);

%% --- Plot Ex and Ey using surf() ---
figure;
subplot(1,2,1);
surf(X, Y, Ex, 'EdgeColor', 'none');
title('Electric Field Ex');
xlabel('X'); ylabel('Y'); zlabel('Ex');
colorbar;
view(3); % 3D view

subplot(1,2,2);
surf(X, Y, Ey, 'EdgeColor', 'none');
title('Electric Field Ey');
xlabel('X'); ylabel('Y'); zlabel('Ey');
colorbar;
view(3); % 3D view

%% --- Plot the Vector Field using quiver() ---
figure;
quiver(X, Y, Ex, Ey, 'k');
title('Electric Field Vector Field');
xlabel('X'); ylabel('Y');
axis equal;
xlim([1, ny]);
ylim([1, nx]);
