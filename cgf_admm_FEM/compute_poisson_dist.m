function [u, u_normalized] = compute_poisson_dist(nodes, segments, elements, pre)
% Assume Dirichlet boundary conditions
dirichlet_segments = segments;
neumann_segments = [];

t = elements;
p = nodes;

[u, u_normalized] = poisson_dist(t, p(:,1:2), dirichlet_segments, neumann_segments, pre);

end
