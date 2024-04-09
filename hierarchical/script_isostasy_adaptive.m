%% Adaptive simulation of the lithosphere in Europe

% Topography data
prb.file = fopen('../data/earth/Earth2014.RET2014.1min.geod.bin');
prb.data = fread(prb.file, [21600, 10800], 'int16', 'ieee-be');
fclose(prb.file);

prb.geod = geo_location('europe'); % geographic coordinates in decimal degrees
prb.data = 1e-3 * rot90(geo_locate(prb.data, prb.geod)); % rock-equivalent topography in km
prb.mpsz = fliplr(size(prb.data)); % topography map size in px
prb.grsz = [min(size(prb.data)), min(size(prb.data))]; % grid size in px
prb.dmsz = prb.mpsz ./ prb.grsz; % domain size in ul

% Problem parameters
prb.er = 6371; % Earth radius in km
prb.ul = 2 * pi * prb.er * min(prb.geod.sz) / 360; % km per unit length
prb.km = 1 / prb.ul; % unit length per km
prb.em = 65 / prb.km; % Young's modulus in 1e12kg/ul/s^2
prb.nu = 0.25; % Poisson's ratio -1 <= nu <= 0.5
prb.th = 16 * prb.km; % effective elastic thickness in ul
prb.cr = 30 * prb.km; % standard crustal thickness in ul
prb.gr = 9.81e-3 * prb.km; % gravitational acceleration in ul/s^2
prb.wr = 2.67 / prb.km^3 * prb.gr; % specific reference weight of rock in 1e12kg/ul^2/s^2
prb.wm = 3.33 / prb.km^3 * prb.gr; % specific weight of the upper mantle in 1e12kg/ul^2/s^2
prb.rt = (prb.wm - prb.wr) / prb.wm; % crustal depth-to-thickness ratio
prb.dp = -prb.cr * prb.rt; % initial depth for the mid-surface of the lithosphere in ul

prb.em = @(x, y) prb.em * ones(size(x));
prb.nu = @(x, y) prb.nu * ones(size(x));
prb.th = @(x, y) prb.th * ones(size(x));

prb.tp = @(x, y) eval_uv(prb.data, x, y, prb.grsz) * prb.km; % rock-equivalent topography in ul
prb.lf = @(x, y) -prb.wr * (prb.tp(x, y) - prb.dp); % topographic weight in 1e12kg/ul/s^2
prb.bf = @(x, y) (prb.wm - prb.wr) * ones(size(x)); % specific weight of displaced mass in 1e12kg/ul^2/s^2

% Geometry parameters
mdl.nrb = nrbtform(nrb_unit_square, vecscale([prb.dmsz, 1]));

% Simulation parameters
deg = 4;
sub = 1;
ref = 2;

sim.deg = [deg, deg]; % spline degree
sim.reg = [deg - 2, deg - 2]; % spline regularity
sim.sub = [sub, sub] .* prb.dmsz; % number of subdivisions
sim.ngp = [deg + 1, deg + 1]; % number of Gaussian quadrature points

sim.typ = 'simplified'; % type of hierarchical B-splines
sim.thb = 0; % truncated hierarchical B-splines
sim.ref = [ref, ref]; % number of subdivisions for the refinement
sim.num = 12; % number of refinement steps

% Visualization parameters
vis.res = [501, 501]; % resolution per patch

% Adaptive finite element simulation
sim = sim_plate_adaptive(prb, mdl, sim);

% Visualization of the results
vis_mesh_hier(sim.msh_hrc);
vis_plate_hier(prb, sim, vis);
