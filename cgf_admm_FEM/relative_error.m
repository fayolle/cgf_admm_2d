% Relative error between k iteration starting from the Poisson dist
% and the exact dist.
function relative_error(input_basename, output_basename, k, r, nupdate)

% Read the mesh
nodes_filename = strcat(input_basename, '.1.node');
nodes = ReadNodes(nodes_filename);

% We need only the x and y coordinates
nodes = nodes(:, 1:2);

elements_filename = strcat(input_basename, '.1.ele');
elements = ReadElements(elements_filename);

segments_filename = strcat(input_basename, '.1.poly');
segments = ReadPoly(segments_filename);

% unsigned distance from the nodes to the boundary
dist = dist_to_polygon(nodes, segments);

% Precompute shared structures
pre = precompute(nodes, elements);

% Compute initial solution (should I use the normalized solution?)
[~, d0_normalized] = compute_poisson_dist(nodes, segments, elements, pre);

% Compute iteratively a refined solution for order k
dk = cgf_admm(nodes, segments, elements, d0_normalized, k, r, nupdate, pre);

% Compute the relative error
relative_error = compute_relative_error(dk, dist);

% Some stats:
fprintf("Min relative error: %f\n", min(relative_error(:)));
fprintf("Max relative error: %f\n", max(relative_error(:)));
fprintf("Mean relative error: %f\n", mean(relative_error(:)));

% Save field
output_filename = strcat(output_basename, '.vtk');
save_2d_field_as_vtk(nodes, elements, relative_error, 'node', output_filename, 'Relative error');

end


function rel_err = compute_relative_error(dk, dist, epsilon)
if nargin ~= 3
    epsilon = 1e-7;
end

rel_err = abs(dk - dist);
rel_err = rel_err ./ (abs(dist) + epsilon);

end

