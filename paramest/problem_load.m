% Compute descent direction for the topographic load as design variable
function [obj, dsc, nrm, prb, mdl, sim] = problem_load(dsn, prb, mdl, sim, opt)
    prb.fr = @(x, y) (prb.th(x, y).^3 / 12) .* prb.em(x, y) ./ (1 - prb.nu(x, y).^2); % flexural rigidity
    prb.ld = @(x, y) prb.fr(x, y) .* prb.nu(x, y); % modified Lame parameter
    prb.m2 = @(x, y) prb.fr(x, y) .* (1 - prb.nu(x, y)); % two times shear modulus
    
    prb.tp = @(x, y) reshape(eval_sp(dsn, opt.spc, geo_load(opt.nrb), [x(:), y(:)]), size(x)); % topographic elevation in ul
    prb.lf = @(x, y) -prb.wr * (prb.tp(x, y) - prb.dp); % topographic weight in kg/ul/s^2
    prb.bf = @(x, y) (prb.wm - prb.wr) * ones(size(x)); % specific weight of displaced mass in kg/ul^2/s^2
    
    sim = compute_state(prb, mdl, sim);
    
    prb.wd = @(x, y) reshape(eval_sp(sim.sol, sim.spc{1}, geo_load(sim.nrb{1}), [x(:), y(:)]), size(x)) + prb.dp;
    prb.tt = @(x, y) prb.wr * (prb.wd(x, y) - prb.mh(x, y));
    sim = compute_adjoint_state(prb, mdl, sim);
    
    [sim.obj, sim.drv] = objective_load(prb, sim);
    sim = compute_gradient(mdl, sim);
    
    obj = sim.obj;
    dsc = -sim.grd;
    nrm = sim.nrm;
end
