% Set up multi-patch C^1 coupling
function [csm, lnk] = setup_c1c(sim, mdl)
    csm = zeros(sim.ndf);
    lnk = 1 : sim.ndf;
    for fc = 1 : numel(mdl.fcs)
        ifc = mdl.fcs{fc};
        pca = ifc(1, 1);
        pcb = ifc(2, 1);
        sda = ifc(1, 2);
        sdb = ifc(2, 2);
        
        msh = msh_eval_boundary_side(sim.msh{pca}, sda);
        bda = msh_boundary_side_from_interior(sim.msh{pca}, sda);
        bdb = msh_boundary_side_from_interior(sim.msh{pcb}, sdb);
        spa = sp_precompute(sim.spc{pca}.constructor(bda), bda, 'value', true, 'gradient', true);
        spb = sp_precompute(sim.spc{pcb}.constructor(bdb), bdb, 'value', true, 'gradient', true);
        
        if mdl.flp(fc)
            caa = op_gradu_n_gradv_n(spa, spa, msh, 1);
            cba = op_gradu_n_gradv_n_flip(spa, spb, msh, 1);
            cab = cba.';
            cbb = op_gradu_n_gradv_n(spb, spb, msh, 1);
        else
            caa = op_gradu_n_gradv_n(spa, spa, msh, 1);
            cba = op_gradu_n_gradv_n(spa, spb, msh, 1);
            cab = cba.';
            cbb = op_gradu_n_gradv_n(spb, spb, msh, 1);
        end
        
        csm(sim.dfs{pca}, sim.dfs{pca}) = csm(sim.dfs{pca}, sim.dfs{pca}) + caa;
        csm(sim.dfs{pca}, sim.dfs{pcb}) = csm(sim.dfs{pca}, sim.dfs{pcb}) - cab;
        csm(sim.dfs{pcb}, sim.dfs{pca}) = csm(sim.dfs{pcb}, sim.dfs{pca}) - cba;
        csm(sim.dfs{pcb}, sim.dfs{pcb}) = csm(sim.dfs{pcb}, sim.dfs{pcb}) + cbb;
        
        dfa = sim.spc{pca}.boundary(sda).dofs;
        dfb = sim.spc{pcb}.boundary(sdb).dofs;
        if mdl.flp(fc)
            dfb = reshape(flipud(reshape(dfb, [], sim.spc{pcb}.ncomp)), 1, []);
        end
        itf = [sim.dfs{pca}(dfa); sim.dfs{pcb}(dfb)].';
        for lk = 1 : size(itf, 1)
            mi = min(lnk(itf(lk, :)));
            ma = max(lnk(itf(lk, :)));
            if mi < ma
                lnk(lnk == ma) = mi;
            end
        end
    end
end
