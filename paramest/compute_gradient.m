% Compute gradient of the objective function
function sim = compute_gradient(mdl, sim)
    % Finite element assembly
    sim.shm = sparse([]);
    for pc = 1 : numel(mdl.nrb)
        nrb = nrbrefine(mdl.nrb{pc}, sim.deg, sim.reg, sim.sub);
        msh = setup_msh(nrb, sim.ngp);
        spc = sp_nurbs(nrb, msh);
        
        shm = op_u_v_tp(spc, spc, msh); % + op_gradu_gradv_tp(spc, spc, msh);
        sim.shm = [sim.shm; shm];
    end
    clear pc nrb msh spc shp;
    
    % Projection onto constraint space
    mat = sim.nsm' * sim.shm * sim.nsm;
    rhs = sim.nsm' * sim.drv;
    vec = mat \ rhs;
    
    % Solution vector
    sim.grd = sim.nsm * vec;
    sim.nrm = sqrt(vec' * mat * vec);
end
