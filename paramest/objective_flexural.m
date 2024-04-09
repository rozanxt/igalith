% Integrated squared error with flexural parameter as design variable
function [obj, drv] = objective_flexural(prb, sim)
    mev = msh_precompute(sim.msh{1});
    wjs = mev.quad_weights .* mev.jacdet;
    
    wds = squeeze(prb.wd(mev.quad_nodes(1, :, :), mev.quad_nodes(2, :, :)));
    mds = squeeze(prb.mh(mev.quad_nodes(1, :, :), mev.quad_nodes(2, :, :)));
    obj = 0.5 * sum(wjs .* (wds - mds).^2, 'all');
    
    if nargout >= 2
        sev = sp_precompute_param(sim.spc{1}, sim.msh{1});
        
        val = sp_eval(sim.sol, sim.spc{1}, geo_load(sim.nrb{1}), {mev.qn{1}(:), mev.qn{2}(:)}, {'value', 'hessian', 'laplacian'});
        val{1} = reshape(squeeze(permute(reshape(val{1}, [mev.nqn_dir(1), mev.nel_dir(1), mev.nqn_dir(2), mev.nel_dir(2)]), [1, 3, 2, 4])), [mev.nqn, mev.nel]);
        val{2} = reshape(squeeze(permute(reshape(val{2}, [2, 2, mev.nqn_dir(1), mev.nel_dir(1), mev.nqn_dir(2), mev.nel_dir(2)]), [1, 2, 3, 5, 4, 6])), [2, 2, mev.nqn, mev.nel]);
        val{3} = reshape(squeeze(permute(reshape(val{3}, [mev.nqn_dir(1), mev.nel_dir(1), mev.nqn_dir(2), mev.nel_dir(2)]), [1, 3, 2, 4])), [mev.nqn, mev.nel]);
        
        adj = sp_eval(sim.adj, sim.spc{1}, geo_load(sim.nrb{1}), {mev.qn{1}(:), mev.qn{2}(:)}, {'value', 'hessian', 'laplacian'});
        adj{1} = reshape(squeeze(permute(reshape(adj{1}, [mev.nqn_dir(1), mev.nel_dir(1), mev.nqn_dir(2), mev.nel_dir(2)]), [1, 3, 2, 4])), [mev.nqn, mev.nel]);
        adj{2} = reshape(squeeze(permute(reshape(adj{2}, [2, 2, mev.nqn_dir(1), mev.nel_dir(1), mev.nqn_dir(2), mev.nel_dir(2)]), [1, 2, 3, 5, 4, 6])), [2, 2, mev.nqn, mev.nel]);
        adj{3} = reshape(squeeze(permute(reshape(adj{3}, [mev.nqn_dir(1), mev.nel_dir(1), mev.nqn_dir(2), mev.nel_dir(2)]), [1, 3, 2, 4])), [mev.nqn, mev.nel]);
        
        nus = squeeze(prb.nu(mev.quad_nodes(1, :, :), mev.quad_nodes(2, :, :)));
        
        drv = zeros(sev.ndof, 1);
        for el = 1 : mev.nel
            for qn = 1 : mev.nqn
                wj = wjs(qn, el);
                for sh = 1 : sev.nsh(el)
                    dof = sev.connectivity(sh, el);
                    sfv = sev.shape_functions(qn, sh, el);
                    wfl = val{3}(qn, el);
                    zfl = adj{3}(qn, el);
                    wfh = val{2}(:, :, qn, el);
                    zfh = adj{2}(:, :, qn, el);
                    nu = nus(qn, el);
                    drv(dof) = drv(dof) - wj * sfv * (nu * wfl * zfl + (1 - nu) * sum(wfh .* zfh, 'all'));
                end
            end
        end
    end
    
    obj = obj * prb.ul^2;
end
