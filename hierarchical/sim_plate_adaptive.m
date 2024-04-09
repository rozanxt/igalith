% Adaptive finite element simulation of a linearized Kirchhoff-Love plate
function sim = sim_plate_adaptive(prb, mdl, sim)
    % Material parameters
    prb.fr = @(x, y) (prb.th(x, y).^3 / 12) .* prb.em(x, y) ./ (1 - prb.nu(x, y).^2); % flexural rigidity
    prb.ld = @(x, y) prb.fr(x, y) .* prb.nu(x, y); % modified Lame parameter
    prb.m2 = @(x, y) prb.fr(x, y) .* (1 - prb.nu(x, y)); % two times shear modulus
    
    % Refinement procedure
    sim.nrb = nrbrefine(mdl.nrb, sim.deg, sim.reg, sim.sub);
    sim.msh = setup_msh(sim.nrb, sim.ngp);
    sim.spc = sp_nurbs(sim.nrb, sim.msh);
    
    sim.msh_hrc = hierarchical_mesh(sim.msh, sim.ref);
    sim.spc_hrc = hierarchical_space(sim.msh_hrc, sim.spc, sim.typ, sim.thb, sim.reg);
    
    sim.mrk_elm = {[]};
    
    for k = 1 : sim.num
        [sim.msh_hrc, sim.new_cel] = hmsh_refine(sim.msh_hrc, sim.mrk_elm);
        sim.mrk_fcn = compute_functions_to_deactivate(sim.msh_hrc, sim.spc_hrc, sim.mrk_elm, 'elements');
        sim.spc_hrc = hspace_refine(sim.spc_hrc, sim.msh_hrc, sim.mrk_fcn, sim.new_cel);
        
        sim.sol = solve_plate_hier(sim.spc_hrc, sim.msh_hrc, prb);
        sim.est = adaptivity_estimate_multilevel(sim.sol, sim.spc_hrc, sim.msh_hrc, prb);
        
        adt.flag = 'elements';
        adt.mark_strategy = 'MS';
        adt.mark_param = 0.5;
        [sim.mrk_elm, sim.mrk_num] = adaptivity_mark(sim.est, sim.msh_hrc, sim.spc_hrc, adt);
    end
end
