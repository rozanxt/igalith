% Visualize lithospheric depression in Central Java
function vis_central_java(prb, mdl, sim, vis)
    figure;
    for pc = 1 : numel(mdl.nrb)
        sol = sim.sol(sim.dfs{pc});
        spc = sim.spc{pc};
        geo = geo_load(sim.nrb{pc});
        [val, grd] = sp_eval(sol, spc, geo, vis.res);
        [X, Y] = deal(squeeze(grd(1, :, :)), squeeze(grd(2, :, :)));
        [lon, lat] = geo_coords(X, Y, prb.geod, prb.dmsz);
        contour(lon, lat, (val + prb.dp) * prb.ul, -31.6 : 0.2 : -30, 'ShowText', true, 'LabelSpacing', Inf, 'LineWidth', 2);
        hold on;
    end
    
    axis equal;
    title('Lithospheric depression (Vening-Meinesz)');
    xlabel('Longitude');
    ylabel('Latitude');
    ax = gca;
    ax.XAxis.TickLabelFormat = '%g°';
    ax.YAxis.TickLabelFormat = '%g°';
    cb = colorbar;
    cb.Label.String = 'Elevation (km)';
    clim([-31.6, -30]);
end
