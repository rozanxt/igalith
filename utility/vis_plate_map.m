% Visualize lithospheric depression for the plate model
function vis_plate_map(prb, mdl, sim, vis)
    figure;
    for pc = 1 : numel(mdl.nrb)
        sol = sim.sol(sim.dfs{pc});
        spc = sim.spc{pc};
        geo = geo_load(sim.nrb{pc});
        [val, grd] = sp_eval(sol, spc, geo, vis.res);
        [X, Y] = deal(squeeze(grd(1, :, :)), squeeze(grd(2, :, :)));
        [lon, lat] = geo_coords(X, Y, prb.geod, prb.dmsz);
        contour(lon, lat, (val + prb.dp) * prb.ul, 'ShowText', true, 'LabelSpacing', Inf, 'LineWidth', 2);
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
end
