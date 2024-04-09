% Visualize lithospheric depression for the plate model on a hierarchical mesh
function vis_plate_hier(prb, sim, vis)
    figure;
    [val, grd] = sp_eval(sim.sol, sim.spc_hrc, geo_load(sim.nrb), vis.res);
    [X, Y] = deal(squeeze(grd(1, :, :)), squeeze(grd(2, :, :)));
    [lon, lat] = geo_coords(X, Y, prb.geod, prb.dmsz);
    contour(lon, lat, (val + prb.dp) * prb.ul, 'ShowText', true, 'LabelSpacing', Inf, 'LineWidth', 2);
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
