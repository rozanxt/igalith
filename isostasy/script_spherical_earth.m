%% Global simulation of Earth's lithosphere as a spherical shell

% Topography data
prb.file = fopen('../data/earth/Earth2014.RET2014.1min.geod.bin');
prb.data = fread(prb.file, [21600, 10800], 'int16', 'ieee-be');
fclose(prb.file);

prb.mpsz = 10 * [360, 180]; % topography map size in px
prb.data = 1e-3 * rot90(imresize(prb.data, prb.mpsz)); % rock-equivalent topography in km

% Problem parameters
prb.er = 6371; % Earth radius in km
prb.ul = prb.er; % km per unit length
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

prb.em = @(x, y, z) prb.em * ones(size(x));
prb.nu = @(x, y, z) prb.nu * ones(size(x));

prb.rd = @(x, y, z) sqrt(x.^2 + y.^2 + z.^2); % radial distance in ul
prb.sr = @(x, y, z) prb.ul * (prb.rd(x, y, z) - 1); % radial displacement in km
prb.tp = @(x, y, z) eval_uv(prb.data, (1 + atan2(y, x) / pi) / 2, 0.5 + asin(min(max(z, -1), 1)) / pi, prb.mpsz) * prb.km; % rock-equivalent topography in ul
prb.fm = @(x, y, z) -prb.wr * (prb.tp(x, y, z) - prb.dp); % topographic load in 1e12kg/ul/s^2
prb.lf = @(x, y, z) permute(repmat(prb.fm(x, y, z) ./ prb.rd(x, y, z), [1, 1, 3]) .* cat(3, x, y, z), [3, 1, 2]); % radial gravitational force in 1e12kg/ul/s^2
prb.bf = @(x, y, z) (prb.wm - prb.wr) * ones(size(x)); % specific weight of displaced mass in 1e12kg/ul^2/s^2

% Geometry parameters
mdl = mdl_quad_sphere((1 + prb.dp) * [1, 1, 1]);

% Simulation parameters
deg = 4;
sub = 8;

sim.deg = [deg, deg]; % spline degree
sim.reg = [deg - 2, deg - 2]; % spline regularity
sim.sub = [sub, sub]; % number of subdivisions
sim.ngp = [deg + 1, deg + 1]; % number of Gaussian quadrature points

% Visualization parameters
vis.res = [301, 301]; % resolution per patch

% Finite element simulation
sim = sim_shell_isostasy(prb, mdl, sim);

% Visualization of the results
vis_earth_map(prb, mdl, vis);
vis_shell_map(prb, mdl, sim, vis);
