% X is needed to compute the volume force:
% compute_volume_force = div X
% I need to compute an approximation of div X on a triangle mesh defined
% by t and p
function phi = solve_poisson_div(triangles, coordinates, X, dirichlet, neumann, pre)

free_nodes=setdiff(unique(triangles), unique(dirichlet));

A = pre.A;

%b = sparse(size(coordinates, 1), 1);
divx = compute_integrated_divergence(triangles, coordinates, X);
b = -divx;

% Neumann conditions
for j = 1 : size(neumann, 1)
    b(neumann(j, :))=b(neumann(j, :)) + ...
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

phi = u;
end


function divx = compute_integrated_divergence(T, P, X)
% given a vector field X sampled at each face
% compute the integrated divergence of this vector field associated with a
% given vertex

num_vert = size(P, 1);
divx = zeros(num_vert, 1);

% pre-compute the face neighboring information for each vertex
vertex_face_ring = compute_vertex_face_ring(T);

for v = 1:num_vert
    % given a vertex index v, find the adjacent triangles
    neigh_faces = vertex_face_ring{v};
    num_tri = size(neigh_faces,2);

    for t = 1:num_tri
        tri_div = compute_integrated_divergence_triangle(v, T(neigh_faces(t),:), P, X(neigh_faces(t),:));
        divx(v) = divx(v) + tri_div;
    end
end

divx = 0.5 * divx;

end


function divx = compute_integrated_divergence_triangle(i, tri, P, Xi)
% order the vertices: i, j, k are given in ccw order
if i == tri(1)
    j = tri(2);
    k = tri(3);
elseif i == tri(2)
    j = tri(3);
    k = tri(1);
elseif i == tri(3)
    j = tri(1);
    k = tri(2);
end

% i is the vertex at which we compute the integrated divergence
e1 = P(j,:) - P(i,:);
e2 = P(k,:) - P(i,:);

% compute the angle
v1 = P(i,:) - P(k,:);
v2 = P(j,:) - P(k,:);
theta1 = compute_angle(v1, v2);
v1 = P(k,:) - P(j,:);
v2 = P(i,:) - P(j,:);
theta2 = compute_angle(v1, v2);

% integrated divergence associated to i for the current triangle
divx = cot(abs(theta1)) * dot(e1, Xi) + cot(abs(theta2)) * dot(e2, Xi);

end


function d = dot(v1, v2)
d = v1(:,1).*v2(:,1) + v1(:,2).*v2(:,2);
end


function c = cross(v1, v2)
cz = v1(:,1).*v2(:,2) - v2(:,1).*v1(:,2);
c = [0 0 cz];
end


function theta = compute_angle(v1, v2)
theta = atan2(norm(cross(v1,v2)), dot(v1,v2));
end



% Neumann boundary condition
function neumann_bc = compute_neumann(nodes)
neumann_bc = -ones(size(nodes, 1), 1);
end


% Dirichlet boundary condition
function dirichlet_boundary_value = compute_dirichlet(nodes)
dirichlet_boundary_value = zeros(size(nodes, 1), 1);
end
