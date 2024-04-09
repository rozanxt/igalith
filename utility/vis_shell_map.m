% Visualize radial elevation for the spherical shell model
function vis_shell_map(prb, mdl, sim, vis)
    figure;
    for pc = 1 : numel(mdl.nrb)
        sol = sim.sol(sim.dfs{pc});
        spc = sim.spc{pc};
        geo = geo_load(sim.nrb{pc});
        dfm = geo_deform(sol, spc, geo);
        pts = eval_nrb(dfm.nurbs, vis.res);
        [X, Y, Z] = deal(squeeze(pts(1, :, :)), squeeze(pts(2, :, :)), squeeze(pts(3, :, :)));
        surf(X, Y, Z, prb.sr(X, Y, Z));
        hold on;
    end
    
    shading interp;
    axis equal;
    view(-162, 14);
    title('Deformation of the lithosphere');
    cb = colorbar;
    cb.Label.String = 'Elevation (km)';
end
