% Visualize geometry of the mesh
function vis_mesh_map(prb, mdl, sim, varargin)
    ip = inputParser;
    addParameter(ip, 'grid', false);
    parse(ip, varargin{:});
    
    for pc = 1 : numel(mdl.nrb)
        hold on;
        
        knt = mdl.nrb{pc}.knots;
        order = mdl.nrb{pc}.order;
        
        if ip.Results.grid
            p = nrbeval(mdl.nrb{pc}, {linspace(knt{1}(order(1)), knt{1}(end - order(1) + 1), sim.sub(1) + 1)...
                                      linspace(knt{2}(order(2)), knt{2}(end - order(2) + 1), sim.sub(2) + 1)});
            [X, Y, Z] = deal(squeeze(p(1, :, :)), squeeze(p(2, :, :)), squeeze(p(3, :, :)));
            [lon, lat] = geo_coords(X, Y, prb.geod, prb.dmsz);
            surf(lon, lat, Z, 'FaceColor', 'none', 'EdgeColor', '#808080');
        end
        
        p = nrbeval(mdl.nrb{pc}, {linspace(knt{1}(order(1)), knt{1}(end - order(1) + 1), 2)...
                                  linspace(knt{2}(order(2)), knt{2}(end - order(2) + 1), 2)});
        [X, Y, Z] = deal(squeeze(p(1, :, :)), squeeze(p(2, :, :)), squeeze(p(3, :, :)));
        [lon, lat] = geo_coords(X, Y, prb.geod, prb.dmsz);
        surf(lon, lat, Z, 'FaceColor', 'none', 'LineWidth', 1.5);
    end
end
