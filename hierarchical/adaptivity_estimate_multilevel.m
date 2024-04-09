% Multi-level estimator for adaptive local refinement
function est = adaptivity_estimate_multilevel(sol, spc_hrc, msh_hrc, prb)
    npl = cumsum([0, msh_hrc.nel_per_level]);
    mrk_est = cell(msh_hrc.nlevels, 1);
    for lv = 1 : msh_hrc.nlevels
        if msh_hrc.msh_lev{lv}.nel > 0
            ie = (npl(lv) + 1) : npl(lv + 1);
            [~, ind, ~] = intersect(ie, 1 : msh_hrc.nel);
            mrk_est{lv} = msh_hrc.active{lv}(ind);
        end
    end
    
    [msh_est, new_cel] = hmsh_refine(msh_hrc, mrk_est);
    mrk_fcn = compute_functions_to_deactivate(msh_est, spc_hrc, mrk_est, 'elements');
    spc_est = hspace_refine(spc_hrc, msh_est, mrk_fcn, new_cel);
    
    sol_est = solve_plate_hier(spc_est, msh_est, prb);
    
    spc_ifm = hspace_in_finer_mesh(spc_hrc, msh_hrc, msh_est);
    
    [val1, ~] = hspace_eval_hmsh(sol, spc_ifm, msh_est);
    [val2, ~] = hspace_eval_hmsh(sol_est, spc_est, msh_est);
    
    npl = cumsum([0, msh_est.nel_per_level]);
    est = zeros(msh_est.nel, 1);
    for lv = 1 : msh_est.nlevels
        mev = msh_est.msh_lev{lv};
        if mev.nel > 0
            ie = (npl(lv) + 1) : npl(lv + 1);
            wd = mev.quad_weights .* mev.jacdet;
            est(ie) = sum(wd .* (val1(:, ie) - val2(:, ie)).^2);
        end
    end

    agl = [];
    npl = cumsum(msh_est.nel_per_level);
    for lv = 1 : msh_hrc.nlevels
        if msh_hrc.msh_lev{lv}.nel > 0
            % TODO: hmsh_get_children has weird behavior in the multipatch case
            [~, ~, coc] = hmsh_get_children(msh_est, lv, msh_hrc.active{lv});
            [~, ind] = ismember(coc(:), sort(coc(:)));
            agl = [agl; npl(lv) + ind];
        end
    end
    
    est = sum(reshape(est(agl), sum(msh_hrc.nsub), [])).';
end
