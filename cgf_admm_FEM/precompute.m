function pre = precompute(nodes, elements)
% Precompute quantities re-used in the program
% and store them in the structure "pre".
% For now, it contains only the stiffness matrix:
%   pre.A

A = sparse(size(nodes, 1), size(nodes, 1));

% Assembly of the stiffness matrix
for j = 1:size(elements, 1)
    A(elements(j, :), elements(j, :)) = ...
        A(elements(j, :), elements(j, :)) ...
        + compute_local_assembly_matrix(nodes(elements(j, :), :));
end

pre.A = A;

end


function M = compute_local_assembly_matrix(vertices)
d = size(vertices, 2);
G = [ones(1, d+1); vertices'] \ [zeros(1, d); eye(d)];
M = det([ones(1, d+1); vertices']) * G * G' / prod(1:d);
end

