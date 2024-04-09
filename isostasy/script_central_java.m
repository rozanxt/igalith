%% Multi-patch simulation of the lithosphere in Central Java

% Topography data
prb.file = fopen('../data/earth/Earth2014.RET2014.1min.geod.bin');
prb.data = fread(prb.file, [21600, 10800], 'int16', 'ieee-be');
fclose(prb.file);

prb.geod = geo_location('central_java'); % geographic coordinates in decimal degrees
prb.data = 1e-3 * rot90(geo_locate(prb.data, prb.geod)); % rock-equivalent topography in km
prb.mpsz = fliplr(size(prb.data)); % topography map size in px
prb.grsz = [min(size(prb.data)), min(size(prb.data))] / 9; % grid size in px
prb.dmsz = prb.mpsz ./ prb.grsz; % domain size in ul

% Problem parameters
prb.er = 6371; % Earth radius in km
prb.ul = 2 * pi * prb.er * min(prb.geod.sz) / 360 / 9; % km per unit length
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
part = 'single';

switch part
    case 'single'
        mdl.nrb = {nrb_quad_tile([0, 0; 9, 0; 0, 9; 9, 9])};
        mdl.fcs = {};
    case 'partial'
        mdl.nrb = {nrb_quad_tile([0, 2; 2, 2; 0, 7; 4, 7]), ...
                   nrb_quad_tile([2, 2; 5, 1; 4, 7; 7, 8]), ...
                   nrb_quad_tile([5, 1; 9, 0; 7, 8; 9, 8]), ...
                   nrb_quad_tile([4, 7; 7, 8; 5, 9; 6, 9])};
        mdl.fcs = {[1, 2; 2, 1], [2, 2; 3, 1], [2, 4; 4, 3]};
    case 'full'
        mdl.nrb = {nrb_quad_tile([0, 2; 2, 2; 0, 7; 4, 7]), ...
                   nrb_quad_tile([2, 2; 5, 1; 4, 7; 7, 8]), ...
                   nrb_quad_tile([5, 1; 9, 0; 7, 8; 9, 8]), ...
                   nrb_quad_tile([4, 7; 7, 8; 5, 9; 6, 9]), ...
                   nrb_quad_tile([0, 7; 4, 7; 0, 9; 5, 9]), ...
                   nrb_quad_tile([7, 8; 9, 8; 6, 9; 9, 9]), ...
                   nrb_quad_tile([0, 0; 2, 0; 0, 2; 2, 2]), ...
                   nrb_quad_tile([2, 0; 5, 0; 2, 2; 5, 1]), ...
                   nrb_quad_tile([5, 0; 9, 0; 5, 1; 9, 0])};
        mdl.fcs = {[1, 2; 2, 1], ...
                   [2, 2; 3, 1], ...
                   [2, 4; 4, 3], ...
                   [1, 4; 5, 3], ...
                   [4, 1; 5, 2], ...
                   [3, 4; 6, 3], ...
                   [4, 2; 6, 1], ...
                   [1, 3; 7, 4], ...
                   [7, 2; 8, 1], ...
                   [2, 3; 8, 4], ...
                   [3, 3; 9, 4], ...
                   [8, 2; 9, 1]};
end

mdl.flp = false(1, numel(mdl.fcs));
mdl.dbd = [];

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
sim = sim_plate_isostasy(prb, mdl, sim);

% Visualization of the results
vis_central_java(prb, mdl, sim, vis);
vis_mesh_map(prb, mdl, sim, 'grid', true);
