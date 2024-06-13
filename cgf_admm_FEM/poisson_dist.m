% Compute the distance to a surface Gamma approximated by a polyline by:
% * solving the Poisson equation: -\Delta u = 1
%   with Dirichlet boundary condition: u = 0 on Gamma_D
% * normalize the solution u with Tucker normalization

function [U, U_normalized] = poisson_dist(T, P, dirichlet_segments, neumann_segments, pre)
% solve Poisson equation with Dirichlet b.c.: u=0
U = solve_poisson(T, P, dirichlet_segments, neumann_segments, pre);

% compute derivative of U at each node
dU1 = compute_derivative_norm(U, T, P);

% apply Tucker normalization
U_normalized = compute_tucker_normalization(U, dU1);
end


function U_normalized = compute_tucker_normalization(U, dU1)
dU2 = dU1 .* dU1;
U_normalized = 2*U ./ (dU1 + sqrt(dU2 + 2*U));

% if dU and U are 0, U_normalized can be NaN
nan_idx = isnan(U_normalized);
U_normalized(nan_idx) = 0.0;
end


% Compute the L2 norm of the gradient of u at each grid node
function dU1 = compute_derivative_norm(U, T, P)
dU = compute_derivative_at_node(U, T, P);
dU1 = compute_L2_norm(dU);
end


function dU1 = compute_L2_norm(dU)
dU1 = dU(:,1) .* dU(:,1) + dU(:,2) .* dU(:,2);
dU1 = sqrt(dU1);
end
