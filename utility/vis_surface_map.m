% Visualize map data using a surface plot
function vis_surface_map(X, Y, Z, varargin)
    ip = inputParser;
    addParameter(ip, 'fname', '');
    addParameter(ip, 'cname', '');
    addParameter(ip, 'view', '');
    parse(ip, varargin{:});
    
    figure;
    surf(X, Y, Z);
    shading interp;
    axis tight;
    aspect = daspect;
    pbaspect([aspect(1), aspect(2), 0.5 * min(aspect(1), aspect(2))]);
    
    switch ip.Results.view
        case 'top'
            view(0, 90)
        case 'side'
            view(0, 0)
    end
    
    title(ip.Results.fname);
    xlabel('Longitude');
    ylabel('Latitude');
    zlabel(ip.Results.cname);
    ax = gca;
    ax.XAxis.TickLabelFormat = '%g°';
    ax.YAxis.TickLabelFormat = '%g°';
    cb = colorbar;
    cb.Label.String = ip.Results.cname;
end
