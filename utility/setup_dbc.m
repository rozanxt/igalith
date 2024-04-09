% Set up degrees of freedom corresponding to the Dirichlet boundary
function ddf = setup_dbc(spc, dfs, dbd)
    ddf = [];
    switch size(dbd, 2)
        case 2
            for bc = 1:size(dbd, 1)
                pc = dbd(bc, 1);
                sd = dbd(bc, 2);
                bd = spc{pc}.boundary(sd);
                ddf = union(ddf, dfs{pc}(bd.dofs));
            end
        case 3
            for bc = 1:size(dbd, 1)
                pc = dbd(bc, 1);
                sd = dbd(bc, 2);
                cp = dbd(bc, 3);
                bd = spc{pc}.boundary(sd);
                ddf = union(ddf, dfs{pc}(bd.dofs(bd.comp_dofs{cp})));
            end
    end
end
