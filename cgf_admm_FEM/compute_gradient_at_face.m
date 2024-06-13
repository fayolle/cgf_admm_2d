function dU_face = compute_gradient_at_face(U, T, P)
% pre-compute each triangle area
tri_areas = compute_triangle_areas(T, P);

% pre-compute gradient of U at each face
dU_face = compute_derivative_at_triangle(U, T, P, tri_areas);
end


% Compute the area of a triangle defined by three vertices: a, b, c.
% Non-vectorized version.
function s = compute_triangle_area(a, b, c)
acx = a(1) - c(1);
bcx = b(1) - c(1);
acy = a(2) - c(2);
bcy = b(2) - c(2);
two_area = acx * bcy - acy * bcx;
s = 0.5 * two_area;
end


function tri_areas = compute_triangle_areas(T, P)
num_tris = size(T, 1);
tri_areas = zeros(num_tris, 1);
for t = 1:num_tris
    curr_tri = T(t,:);

    % vert indices in ccw order
    i = curr_tri(1);
    j = curr_tri(2);
    k = curr_tri(3);

    % The triangle area:
    A = compute_triangle_area(P(i,:), P(j,:), P(k,:));

    tri_areas(t) = A;
end
end


% From E. Suli notes on FEM p. 37-38, eq. 2.16, 2.18
function dU_face = compute_derivative_at_triangle(U, T, P, tri_areas)
num_tri = size(T, 1);
dU_face = zeros(num_tri, 2);

for t=1:num_tri
    curr_tri = T(t,:);

    % use the same notation as in the linked document above
    i = curr_tri(1);
    j = curr_tri(2);
    k = curr_tri(3);

    % field value at each vertex
    V1 = U(i);
    V2 = U(j);
    V3 = U(k);

    % nodes coordinates
    x1 = P(i, 1);
    y1 = P(i, 2);
    x2 = P(j, 1);
    y2 = P(j, 2);
    x3 = P(k, 1);
    y3 = P(k, 2);

    % The triangle area:
    A = tri_areas(t);

    dudxsi = V2 - V1;
    dudeta = V3 - V1;

    dudx = 1.0 / (2*A) * ((y3-y1)*dudxsi - (y2-y1)*dudeta);
    dudy = 1.0 / (2*A) * (-(x3-x1)*dudxsi + (x2-x1)*dudeta);

    dU_face(t,:) = [dudx dudy];
end
end
