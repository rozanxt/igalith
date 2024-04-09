% Flexural isostasy for a linearized Koiter shell
function sim = sim_shell_isostasy(prb, mdl, sim)
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
        spc = setup_spc(nrb, msh);
        
        ksm = op_KL_shells_tp(spc, spc, msh, prb.em, prb.nu, prb.th);
        bym = op_u_v_tp(spc, spc, msh, prb.bf);
        mat = ksm + bym;
        rhs = op_f_v_tp(spc, msh, prb.lf);
        
        sim.nrb{pc} = nrb;
        sim.msh{pc} = msh;
        sim.spc{pc} = spc;
        
        sim.mat = blkdiag(sim.mat, mat);
        sim.rhs = [sim.rhs; rhs];
        
        sim.dfs{pc} = (sim.ndf + 1) : (sim.ndf + spc.ndof);
        sim.ndf = sim.ndf + spc.ndof;
    end
    clear pc nrb msh spc ksm bym mat rhs zmc;
    
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
