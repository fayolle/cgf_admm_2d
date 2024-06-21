function Idist = cgf_admm(I, niter, r, nupdate)
% cgf_admm(I, niter, r, nupdate)
% Compute an approximate distance to the boundary of the input image I by
% minimizing \int (|p|-1)^2 with constraint p=\nabla u and u=0 on the
% boundary by ADMM. 
%
% I: input (black and white) image
% niter: number of iterations
% r: relaxation parameter
% nupdate: number of iterations until the next Poisson solve is triggered 
% 

[row,col] = find(I > 0);
v = (1:length(row))';
G = zeros(size(I));

for i = 1:length(row)
    G(row(i),col(i)) = v(i);
end

% delsq(G) is the 5-point discrete (negative) Laplacian 
D = delsq(G);

% init with Poisson dist
U = Poisson_dist(I);

% Tucker normalization 
[gx, gy] = gradient(U);
dU2 = gx.*gx + gy.*gy; 
dU1 = sqrt(dU2);
U_normalized = 2.0.*U ./ (dU1 + sqrt(dU2 + 2.0.*U));
nan_idx = isnan(U_normalized);
U_normalized(nan_idx) = 0.0;

%[gx, gy] = gradient(U);
[gx, gy] = gradient(U_normalized);

lambdax = zeros(size(I));
lambday = zeros(size(I));

% admm loop
for i=1:niter
    % step 1: compute p
    qx = gx - lambdax./r;
    qy = gy - lambday./r;
    q = sqrt(qx.*qx+qy.*qy);

    c = (2+r.*q)./((2+r).*q);
    px = c.*qx;
    py = c.*qy;

    px(abs(q)<eps) = 0.0;
    py(abs(q)<eps) = 0.0;

    % step 2: update u
    if mod(i, nupdate)==0
        px2 = px + (1.0/r).*lambdax;
        py2 = py + (1.0/r).*lambday;
        div = -divergence(px2, py2);

        u = D\div(G>0);
        U = G;
        U(G>0) = full(u(G(G>0)));

        % step 3: update lambda
        lambdax = lambdax + r.*(px - gx);
        lambday = lambday + r.*(py - gy);

        [gx, gy] = gradient(U);
    else
        oldgx = gx;
        oldgy = gy;

        % step 2: update grad u
        gx = px + lambdax ./ r;
        gy = py + lambday ./ r;

        % step 3: update lambda
        lambdax = lambdax + r.*(px - oldgx);
        lambday = lambday + r.*(py - oldgy);

    end
end

% last Poisson solve
px2 = px + (1.0/r).*lambdax;
py2 = py + (1.0/r).*lambday;
div = -divergence(px2, py2);

u = D\div(G>0);
U = G;
U(G>0) = full(u(G(G>0)));

Idist = real(U);

end


function Idist = Poisson_dist(I)
[row,col] = find(I > 0);
v = (1:length(row))';
G = zeros(size(I));

for i = 1:length(row)
    G(row(i),col(i)) = v(i);
end

% delsq(G) is the 5-point discrete (negative) Laplacian
D = delsq(G);

% Number of interior points
N = sum(G(:)>0);

% right-hand side = 1
f = ones(N,1);

% solve Poisson dist
u = D\f;

U = G;
U(G>0) = full(u(G(G>0)));

% Poisson-dist normalization
[Ux, Uy]=gradient(U);
dU = sqrt(Ux.^2 + Uy.^2);
U_normalized = -dU + sqrt(2.*U + dU.^2);
U = U_normalized;

Idist = real(U);
end
