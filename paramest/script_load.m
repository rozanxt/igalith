%% Parameter identification of the effective elastic thickness

% Topography data
prb.file = fopen('../data/earth/Earth2014.RET2014.1min.geod.bin');
prb.topo = fread(prb.file, [21600, 10800], 'int16', 'ieee-be');
fclose(prb.file);

prb.geod = geo_location('europe'); % geographic coordinates in decimal degrees
prb.topo = 1e-3 * rot90(geo_locate(prb.topo, prb.geod)); % rock-equivalent topography in km
prb.mpsz = fliplr(size(prb.topo)); % topography map size in px
prb.grsz = [min(size(prb.topo)), min(size(prb.topo))]; % grid size in px
prb.dmsz = prb.mpsz ./ prb.grsz; % domain size in ul

% Mohorovicic depth data
[I, J, prb.moho] = grdread2('../data/europe/Europe_moho_depth_2007.grd');
K = find(I == prb.geod.lon(1)) : (find(I == prb.geod.lon(2)) - 1);
L = find(J == prb.geod.lat(1)) : (find(J == prb.geod.lat(2)) - 1);
prb.moho = flipud(prb.moho(L, K)); % Mohorovicic depth in km
prb.mhsz = fliplr(size(prb.moho)); % Mohorovicic depth map size in px

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

prb.mh = @(x, y) -eval_uv(prb.moho, x, y, prb.mhsz) * prb.km; % Mohorovicic depth in ul

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

% Optimization parameters
deg = 4;
sub = 16;

opt.deg = [deg, deg]; % spline degree
opt.reg = [deg - 2, deg - 2]; % spline regularity
opt.sub = [sub, sub] .* prb.dmsz; % number of subdivisions
opt.ngp = [deg + 1, deg + 1]; % number of Gaussian quadrature points

% Visualization parameters
vis.res = [501, 501]; % resolution per patch

%% Optimization
options.func_val_tol = 1e-4;
options.grad_nrm_tol = 1e-4;
options.max_iter = 1e2;
options.min_step_size = 1e-4;
options.max_step_size = 1e2;
options.init_step_size = 1e-2;
options.sigma = 1e-4;
options.tau = 0.5;
options.output_fcn = @output_fcn;

opt.nrb = nrbrefine(mdl.nrb{1}, opt.deg, opt.reg, opt.sub);
opt.msh = setup_msh(opt.nrb, opt.ngp);
opt.spc = sp_nurbs(opt.nrb, opt.msh);

opt.init_dsn = prb.km * ones(opt.spc.ndof, 1);
opt.final_dsn = descent_method(@(dsn) problem_load(dsn, prb, mdl, sim, opt), opt.init_dsn, options);

[obj, dsc, nrm, prb, mdl, sim] = problem_load(opt.final_dsn, prb, mdl, sim, opt);

%% Visualization

% Topographic load
[val, grd] = sp_eval(opt.final_dsn, opt.spc, geo_load(opt.nrb), vis.res);
[X, Y] = deal(squeeze(grd(1, :, :)), squeeze(grd(2, :, :)));
[lon, lat] = geo_coords(X, Y, prb.geod, prb.dmsz);
vis_contour_map(lon, lat, val * prb.ul, 'fname', 'Estimated topographic load from observed data', 'cname', 'Rock-equivalent topography (km)');

% Lithospheric depression
vis_plate_map(prb, mdl, sim, vis);
