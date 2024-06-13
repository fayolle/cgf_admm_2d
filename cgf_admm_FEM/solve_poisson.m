% Solve poisson equation by finite element method
function U = solve_poisson(triangles, coordinates, dirichlet, neumann, pre)
free_nodes = setdiff(unique(triangles), unique(dirichlet));

% Load the precomputed stiffness matrix
A = pre.A;

% Volume forces
b = sparse(size(coordinates, 1), 1);
for j = 1:size(triangles, 1)
    b(triangles(j, :)) = ...
        b(triangles(j, :)) + ...
        det([1, 1, 1; coordinates(triangles(j, :), :)']) * ...
        compute_volume_force(sum(coordinates(triangles(j, :), :))/3)/6;
end

% Neumann conditions
for j = 1 : size(neumann, 1)
    b(neumann(j, :)) = ...
        b(neumann(j, :)) + ...
        norm(coordinates(neumann(j, 1), :) - ...
        coordinates(neumann(j, 2), :))*...
        compute_neumann(sum(coordinates(neumann(j, :), :))/2)/2;
end

% Dirichlet conditions
u = sparse(size(coordinates, 1), 1);
u(unique(dirichlet)) = compute_dirichlet(coordinates(unique(dirichlet),:));
b = b - A * u;

% Computation of the solution
u(free_nodes) = A(free_nodes, free_nodes) \ b(free_nodes);

U = u;
end


% Helper functions to customize the volume force, the Dirichlet and
% Neumann boundary conditions:

% volume force (i.e. rhs of the poisson problem)
function volume_force = compute_volume_force(x)
volume_force = ones(size(x, 1), 1);
end


% Neumann boundary condition
function neumann_bc = compute_neumann(x)
neumann_bc = -ones(size(x, 1), 1);
end


% Dirichlet boundary condition
function dirichlet_boundary_value = compute_dirichlet(x)
dirichlet_boundary_value = zeros(size(x, 1), 1);
end

