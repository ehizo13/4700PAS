% Grid dimensions
nx = 50;
ny = 50;
V = zeros(ny, nx);

%Set boundary conditions
V(:, 1) = 1 ;%Left boundary
V(:, end) = 0; %Right boundary
V(end, :) = V(end-1, :) ; %Top boundary (Insulating: ∂V/∂y = 0)
V(1, :) = V(2, :) ; %Bottom boundary (Insulating: ∂V/∂y = 0)

max_iter = 1000;  %Set the maximum number of iterations
tolerance = 1e-5;  %Convergence tolerance
error = inf;
iter = 0;

%Iterative Solution using Gaussian method
figure
while error > tolerance && iter < max_iter
    V_older = V;
    for i = 2:nx-1
        for j = 2:ny-1
            V(i,j) = 0.25 * (V_older(i+1,j) + V_older(i-1,j) + V_older(i, j+1) + V_older(i, j-1));
        end
    end

    %Apply Boundary Conditions
    V(:,1) = 1;
    V(:,end) = 0;
    V(1,:) = V(2,:);
    V(end,:) = V(end-1,:);

    %Compute the Error
    error = max(max(abs(V-V_older)));
    iter = iter +1;

    %Display the progress of 100 iters
    if mod(iter, 100) == 0
        fprintf('Iteration: %d, Max error: %.6f\n', iter, error);

        surf(V)
        shading interp
        pause(0.1)
    end
end

% plt.ion()
% fig, ax = plt.subplots()
% for it = 1:max_iter
%     V_new = V;
% 
%     % Update the potential field
%     for i=1:ny-1
%         for j = nx-1
%             V_new(i, j) = 0.25 * (V(i+1, j) + V(i-1, j) + V(i, j+1) + V(i, j-1));
% 
%         end
%     end
%     % Plot the potential
%     ax.clear()
%     c = ax.imshow(V, cmap='jet', origin='lower')
%     fig.colorbar(c)
%     ax.set_title(f"Iteration: {it}")
%     plt.pause(0.01)
%     plt.ioff()

% end



