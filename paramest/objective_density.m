% Integrated squared error with density as design variable
function [obj, drv] = objective_density(prb, sim)
    mev = msh_precompute(sim.msh{1});
    wjs = mev.quad_weights .* mev.jacdet;
    
    wds = squeeze(prb.wd(mev.quad_nodes(1, :, :), mev.quad_nodes(2, :, :)));
    mds = squeeze(prb.mh(mev.quad_nodes(1, :, :), mev.quad_nodes(2, :, :)));
    obj = 0.5 * sum(wjs .* (wds - mds).^2, 'all');
    
    if nargout >= 2
        sev = sp_precompute_param(sim.spc{1}, sim.msh{1});
        
        adj = sp_eval(sim.adj, sim.spc{1}, geo_load(sim.nrb{1}), {mev.qn{1}(:), mev.qn{2}(:)});
        adj = reshape(squeeze(permute(reshape(adj, [mev.nqn_dir(1), mev.nel_dir(1), mev.nqn_dir(2), mev.nel_dir(2)]), [1, 3, 2, 4])), [mev.nqn, mev.nel]);
        
        drv = zeros(sev.ndof, 1);
        for el = 1 : mev.nel
            for qn = 1 : mev.nqn
                wj = wjs(qn, el);
                for sh = 1 : sev.nsh(el)
                    dof = sev.connectivity(sh, el);
                    sfv = sev.shape_functions(qn, sh, el);
                    zfv = adj(qn, el);
                    drv(dof) = drv(dof) - wj * sfv * zfv;
                end
            end
        end
    end
    
    obj = obj * prb.ul^2;
end
