function AddEllipticAtomicArray(a, b, X0, Y0, VX0, VY0, InitDist, Temp, Type)
global C
global x y AtomSpacing
global nAtoms
global AtomType Vx Vy Mass0 Mass1 Mass2

% Choose mass based on Type (0, 1, or 2)
if Type == 0
    Mass = Mass0;
elseif Type == 1
    Mass = Mass1;
elseif Type == 2
    Mass = Mass2;
end

L = (2*a - 1) * AtomSpacing;
W = (2*b - 1) * AtomSpacing;

% xp(1, :) = linspace(-L/2, L/2, 2*a);
% yp(1, :) = linspace(-W/2, W/2, 2*b);

xp = linspace(-L/2, L/2, 2*a);  % 1D array
yp = linspace(-W/2, W/2, 2*b);  % 1D array

numAtoms = 0;

for i = 1:2*a
    for j = 1:2*b
        % Check if the point is inside the ellipse
        if (xp(i)^2 / (a*AtomSpacing)^2 + yp(j)^2 / (b*AtomSpacing)^2) <= 1
            numAtoms = numAtoms + 1;
            x(nAtoms + numAtoms) = xp(i);
            y(nAtoms + numAtoms) = yp(j);
            
            % Assign the atom type
            AtomType(nAtoms + numAtoms) = Type;
        end
    end
end

% Add randomness to positions
x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + X0;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + Y0;

% Set velocities based on temperature and mass
if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp / Mass);
    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

% Adjust velocities to match initial conditions
Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end
