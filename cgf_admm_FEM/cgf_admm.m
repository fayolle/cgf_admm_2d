function uk = cgf_admm(nodes, segments, elements, u0, niter, r, nupdate, pre)
% cgf_admm(nodes, segments, elements, u0, niter, r, nupdate)
% Compute an approximate distance to the boundary of the domain specified
% by the triangulation (nodes, segments, elements) by minimizing
% \int (|p|-1)^2 with constraint p=\nabla u and u=0 on the
% boundary by ADMM.
% Start the iterative process with the initial solution u0.
%
% (nodes, segments, elements): Triangulated domain
% u0: initial solution
% niter: number of iterations
% r: relaxation parameter
% nupdate: number of iterations until the next Poisson solve is triggered
%

% initial solution
ui = u0;
gradui = compute_gradient_at_face(ui, elements, nodes);

% last column contains a boundary marker that needs to be removed
nodes = nodes(:,1:2);

% Lagrange multipliers
num_ele = size(elements, 1);
lambda = zeros(num_ele, 2);

q = zeros(num_ele, 2);
p = zeros(num_ele, 2);
graduip1 = zeros(num_ele, 2);


for i = 1:niter
    % grad u and lambda fixed, optimize for p
    q(:,1) = gradui(:,1) - lambda(:,1) ./ r;
    q(:,2) = gradui(:,2) - lambda(:,2) ./ r;

    norm_q = compute_norm(q);

    p(:,1) = ((2.0 + r .* norm_q) ./ ((2.0 + r) .* norm_q)) .* q(:,1);
    p(:,2) = ((2.0 + r .* norm_q) ./ ((2.0 + r) .* norm_q)) .* q(:,2);

    zeroidx = (norm_q == 0);
    p(zeroidx,:) = 0;

    % Update grad(u_{i+1})
    graduip1(:,1) = p(:,1) + lambda(:,1) ./ r;
    graduip1(:,2) = p(:,2) + lambda(:,2) ./ r;


    % Update Lagrange multipliers
    lambda(:,1) = lambda(:,1) + r .* (p(:,1) - gradui(:,1));
    lambda(:,2) = lambda(:,2) + r .* (p(:,2) - gradui(:,2));

    gradui(:,1) = graduip1(:,1);
    gradui(:,2) = graduip1(:,2);


    % Update ui (and gradui) every 'nupdate' iterations
    if ((mod(i, nupdate) == 0) && i ~= niter)
        % Given grad(u_i), obtain u by solving:
        % delta u = div grad u_i
        dirichlet_segments = segments;
        neumann_segments = [];
        phi = solve_poisson_div(elements, nodes, gradui, dirichlet_segments, neumann_segments, pre);

        % Thresholding to guarantee ui=0 on the boundary (just in case)
        c = compute_threshold(segments, phi);
        ui = phi - c;

        % Update the gradient
        gradui = compute_gradient_at_face(ui, elements, nodes);
    end
end


% Given grad(u_{i+1}), obtain u by solving:
% delta u = div grad u_{i+1}
dirichlet_segments = segments;
neumann_segments = [];
phi = solve_poisson_div(elements, nodes, gradui, dirichlet_segments, neumann_segments, pre);


% Thresholding to guarantee ui=0 on the boundary (just in case)
c = compute_threshold(segments, phi);
uk = phi - c;

end


% Private functions

function c = compute_threshold(segments, phi)
% Compute c = mean(phi(vi)) where vi are points on the input surface
vi = segments(:);
vi = unique(vi(:));
c = mean(phi(vi(:)));
end

