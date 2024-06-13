% Parameters
niter = 50; 
init_sol = 'Poisson_normalized'; % default is 'Poisson_normalized'
r = 1.0; % relaxation parameter. Default: 1
nupdate = 5; % trigger a Poisson solver every nupdate updates. Default: 1 

iterative_dist_cgf('../data/disk', 'disk', 50);
iterative_dist_cgf('../data/square', 'square', 50, init_sol);
iterative_dist_cgf('../data/riderr', 'riderr', 50, init_sol, r, nupdate);
