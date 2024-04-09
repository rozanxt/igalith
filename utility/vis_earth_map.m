% Visualize global topographic load on a spherical Earth
function vis_earth_map(prb, mdl, vis)
    figure;
    for pc = 1 : numel(mdl.nrb)
        pts = eval_nrb(mdl.nrb{pc}, vis.res);
        [X, Y, Z] = deal(squeeze(pts(1, :, :)), squeeze(pts(2, :, :)), squeeze(pts(3, :, :)));
        surf(X, Y, Z, prb.tp(X, Y, Z) * prb.ul);
        hold on;
    end
    
    shading interp;
    axis equal;
    view(-162, 14);
    title('Topographic load');
    cb = colorbar;
    cb.Label.String = 'Rock-equivalent topography (km)';
    clim([min(min(prb.data)), max(max(prb.data))]);
end
