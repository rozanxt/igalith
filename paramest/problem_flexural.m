% Compute descent direction for the flexural parameter as design variable
function [obj, dsc, nrm, prb, mdl, sim] = problem_flexural(dsn, prb, mdl, sim, opt)
    prb.fp = @(x, y) reshape(eval_sp(dsn, opt.spc, geo_load(opt.nrb), [x(:), y(:)]), size(x)); % flexural parameter
    prb.ld = @(x, y) prb.fp(x, y) .* prb.nu(x, y); % modified Lame parameter
    prb.m2 = @(x, y) prb.fp(x, y) .* (1 - prb.nu(x, y)); % two times shear modulus
    
    sim = compute_state(prb, mdl, sim);
    
    prb.wd = @(x, y) reshape(eval_sp(sim.sol, sim.spc{1}, geo_load(sim.nrb{1}), [x(:), y(:)]), size(x)) + prb.dp;
    prb.tt = @(x, y) prb.wr * (prb.wd(x, y) - prb.mh(x, y));
    sim = compute_adjoint_state(prb, mdl, sim);
    
    [sim.obj, sim.drv] = objective_flexural(prb, sim);
    sim = compute_gradient(mdl, sim);
    
    obj = sim.obj;
    dsc = -sim.grd;
    nrm = sim.nrm;
end
