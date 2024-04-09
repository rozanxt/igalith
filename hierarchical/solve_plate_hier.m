% Flexural isostasy for a linearized Kirchhoff-Love plate on a hierarchical mesh
function sol = solve_plate_hier(spc_hrc, msh_hrc, prb)
    lpm = op_laplaceu_laplacev_hier(spc_hrc, spc_hrc, msh_hrc, prb.ld);
    hsm = op_gradgradu_gradgradv_hier(spc_hrc, spc_hrc, msh_hrc, prb.m2);
    bym = op_u_v_hier(spc_hrc, spc_hrc, msh_hrc, prb.bf);
    mat = lpm + hsm + bym;
    rhs = op_f_v_hier(spc_hrc, msh_hrc, prb.lf);
    
    [dbv, ddf] = sp_drchlt_l2_proj(spc_hrc, msh_hrc, @(x, y, ind) zeros(size(x)), []);
    rhs(ddf) = rhs(ddf) - mat(ddf, ddf) * dbv;
    dof = setdiff(1 : spc_hrc.ndof, ddf);
    
    sol = zeros(spc_hrc.ndof, 1);
    sol(dof) = mat(dof, dof) \ rhs(dof);
    sol(ddf) = dbv;
end
