function iterative_dist_cgf(input_basename, output_basename, k, init_sol, r, nupdate)
% Given a file basename (no extension) corresponding to a meshed domain
% compute a smoothed distance field to the mesh boundary computed with
% cgf_admm.
%

if (nargin == 3)
    init_sol = 'Poisson_normalized';
    r = 1.0; % relaxation parameter. See [Tai4] for recommendations
    nupdate = 1; % trigger a Poisson solve every iteration
end

if (nargin == 4)
    r = 1.0; % relaxation parameter. See [Tai4] for recommendations
    nupdate = 1; % trigger a Poisson solve every iteration
end


% Read the mesh
nodes_filename = strcat(input_basename, '.1.node');
nodes = ReadNodes(nodes_filename);

% We need only the x and y coordinates
nodes = nodes(:, 1:2);

elements_filename = strcat(input_basename, '.1.ele');
elements = ReadElements(elements_filename);

segments_filename = strcat(input_basename, '.1.poly');
segments = ReadPoly(segments_filename);


% Precompute shared structures
pre = precompute(nodes, elements);

% Compute initial solution
[d0, d0_normalized] = compute_poisson_dist(nodes, segments, elements, pre);


if (strcmpi(init_sol,'Poisson'))
    % Compute iteratively a refined solution for k iterations
    dk = cgf_admm(nodes, segments, elements, d0, k, r, nupdate, pre);
elseif (strcmpi(init_sol,'Poisson_normalized'))
    % Alternative: start from the normalized solution to delta u = -1
    dk = cgf_admm(nodes, segments, elements, d0_normalized, k, r, nupdate, pre);
end


% Save field
output_filename = strcat(output_basename, '.vtk');
save_2d_field_as_vtk(nodes, elements, dk, 'node', output_filename, 'Iterative normalization');

heightfield_name = strcat(output_basename, '_heightfield.vtk');
save_heightfield_as_vtk(nodes, elements, dk, heightfield_name, 'Iterative normalization');


% matlab plots
X = nodes(:, 1);
Y = nodes(:, 2);
Z = zeros(size(X));

figure,
trisurf(elements, X, Y, Z, full(dk), 'EdgeColor', 'none'),
title('CGF-ADMM');
colormap jet;
campos([0,0,10]);
axis off;
axis equal;

end
