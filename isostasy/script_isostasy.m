%% Airy-Heiskanen and Vening-Meinesz model of isostasy

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
prb.lf = @(x, y) -prb.wr * (prb.tp(x, y) - prb.dp); % topographic load in 1e12kg/ul/s^2
prb.bf = @(x, y) (prb.wm - prb.wr) * ones(size(x)); % specific weight of displaced mass in 1e12kg/ul^2/s^2

% Geometry parameters
mdl.nrb = {nrbtform(nrb_unit_square, vecscale([prb.dmsz, 1]))}; % geometry of the domain
mdl.fcs = {}; % patch interfaces
mdl.flp = false(1, numel(mdl.fcs)); % orientation flips
mdl.dbd = []; % Dirichlet boundary sides

% Simulation parameters
deg = 4;
sub = 16;

sim.deg = [deg, deg]; % spline degree
sim.reg = [deg - 2, deg - 2]; % spline regularity
sim.sub = [sub, sub] .* prb.dmsz; % number of subdivisions
sim.ngp = [deg + 1, deg + 1]; % number of Gaussian quadrature points

% Visualization parameters
vis.res = [501, 501]; % resolution per patch

% Airy-Heiskanen model of local isostasy
iso.tp = @(x, y) eval_uv(prb.data, x, y, prb.grsz); % rock-equivalent topography in km
iso.ct = @(x, y) prb.cr * prb.ul + prb.wm / (prb.wm - prb.wr) * iso.tp(x, y); % crustal thickness in km
iso.mh = @(x, y) -prb.cr * prb.ul - prb.wr / (prb.wm - prb.wr) * iso.tp(x, y); % Mohorovicic depth in km

% Vening-Meinesz model of flexural isostasy
sim = sim_plate_isostasy(prb, mdl, sim);

% Visualization of the results
vis_airy_heiskanen(prb, iso, 'moho');
vis_plate_map(prb, mdl, sim, vis);
