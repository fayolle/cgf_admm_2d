% Relative error between k iteration starting from the Poisson dist
% and the exact dist.

function value_error(input_basename, output_basename, k, r, nupdate)

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

% Compute initial solution 
[~, d0_normalized] = compute_poisson_dist(nodes, segments, elements, pre);

% Compute iteratively a refined solution for order k
dk = cgf_admm(nodes, segments, elements, d0_normalized, k, r, nupdate, pre);

% Compute the (absolute) value error
value_error = abs(dk - dist);

% Some stats:
fprintf("Min value error: %f\n", min(value_error(:)));
fprintf("Max value error: %f\n", max(value_error(:)));
fprintf("Mean value error: %f\n", mean(value_error(:)));

% Save value error field for visualization
output_filename = strcat(output_basename, '.vtk');
save_2d_field_as_vtk(nodes, elements, value_error, 'node', output_filename, 'Value error');

end
