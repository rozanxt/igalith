% Compute solution to the adjoint state equation
function sim = compute_adjoint_state(prb, mdl, sim)
    % Assemble right-hand side
    sim.ttr = sparse([]);
    for pc = 1 : numel(mdl.nrb)
        nrb = nrbrefine(mdl.nrb{pc}, sim.deg, sim.reg, sim.sub);
        msh = setup_msh(nrb, sim.ngp);
        spc = sp_nurbs(nrb, msh);
        
        ttr = op_f_v_tp(spc, msh, prb.tt);
        sim.ttr = [sim.ttr; ttr];
    end
    clear pc nrb msh spc ttr;
    
    % Projection onto constraint space
    mat = sim.nsm' * sim.mat * sim.nsm;
    rhs = sim.nsm' * sim.ttr;
    vec = mat \ rhs;
    
    % Solution vector
    sim.adj = sim.nsm * vec;
end
