% Grid Size
nx = 50; 
ny = 50; 
N = nx * ny;  % Total number of grid points
% Sparse matrix initialization
G = sparse(N, N);
% Constructing the finite-difference matrix with modified region
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;  % Convert 2D index to 1D index
        
        if i == 1 || i == nx || j == 1 || j == ny
            % Boundary Condition: u = 0
            G(n, :) = 0;    
            G(n, n) = 1;    
        else
            % Interior points
            % if i > 10 && i < 20 && j > 10 && j < 20
            %    G(n, n) = -2;  % Modified region (weaker Laplacian)
            % else
                G(n, n) = -4;  % Standard Laplacian coefficient
            % end
            G(n, n - 1) = 1;   % Left neighbor
            G(n, n + 1) = 1;   % Right neighbor
            G(n, n - ny) = 1;  % Bottom neighbor
            G(n, n + ny) = 1;  % Top neighbor
        end
    end
end
% Compute the 9 smallest eigenvalues and eigenvectors
[E, D] = eigs(G, 9, 'SM');
% Display the 9 smallest eigenvalues
disp('Smallest 9 Eigenvalues:');
disp(diag(D));
% Plot the first 9 eigenmodes
figure;
for k = 1:9
    subplot(3, 3, k);
    imagesc(reshape(E(:, k), ny, nx)); % Reshape eigenvector to 2D
    colormap jet; colorbar;
    title(['Mode ', num2str(k), ', \lambda = ', num2str(D(k, k), '%.4f')]);
end

