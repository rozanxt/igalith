% Compute solution to the state equation
function sim = compute_state(prb, mdl, sim)
    % Finite element assembly
    sim.nrb = cell(1, numel(mdl.nrb));
    sim.msh = cell(1, numel(mdl.nrb));
    sim.spc = cell(1, numel(mdl.nrb));
    sim.mat = sparse([]);
    sim.rhs = sparse([]);
    sim.dfs = cell(1, numel(mdl.nrb));
    sim.ndf = 0;
    for pc = 1 : numel(mdl.nrb)
        nrb = nrbrefine(mdl.nrb{pc}, sim.deg, sim.reg, sim.sub);
        msh = setup_msh(nrb, sim.ngp);
        spc = sp_nurbs(nrb, msh);
        
        lpm = op_laplaceu_laplacev_tp(spc, spc, msh, prb.ld);
        hsm = op_gradgradu_gradgradv_tp(spc, spc, msh, prb.m2);
        bym = op_u_v_tp(spc, spc, msh, prb.bf);
        mat = lpm + hsm + bym;
        rhs = op_f_v_tp(spc, msh, prb.lf);
        
        sim.nrb{pc} = nrb;
        sim.msh{pc} = msh;
        sim.spc{pc} = spc;
        
        sim.mat = blkdiag(sim.mat, mat);
        sim.rhs = [sim.rhs; rhs];
        
        sim.dfs{pc} = (sim.ndf + 1) : (sim.ndf + spc.ndof);
        sim.ndf = sim.ndf + spc.ndof;
    end
    clear pc nrb msh spc lpm hsm bym mat rhs;
    
    % Constraints and boundary conditions
    [sim.csm, sim.lnk] = setup_c1c(sim, mdl);
    sim.ddf = setup_dbc(sim.spc, sim.dfs, mdl.dbd);
    
    % Unconstrained degrees of freedom
    sim.dof = unique(sim.lnk);
    sim.dof = setdiff(sim.dof, sim.ddf);
    
    % Computation of minimal coordinates
    idm = speye(sim.ndf);
    sim.ns0 = idm(sim.lnk, :);
    sim.ns0 = sim.ns0(:, sim.dof);
    sim.ns1 = null(sim.ns0' * sim.csm * sim.ns0);
    sim.nsm = sim.ns0 * sim.ns1;
    
    % Projection onto constraint space
    mat = sim.nsm' * sim.mat * sim.nsm;
    rhs = sim.nsm' * sim.rhs;
    vec = mat \ rhs;
    
    % Solution vector
    sim.sol = sim.nsm * vec;
end
