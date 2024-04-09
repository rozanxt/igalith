% Visualize map data using a contour plot
function vis_contour_map(X, Y, Z, varargin)
    ip = inputParser;
    addParameter(ip, 'fname', '');
    addParameter(ip, 'cname', '');
    addParameter(ip, 'ShowText', true);
    addParameter(ip, 'LabelSpacing', Inf);
    addParameter(ip, 'LineWidth', 2);
    parse(ip, varargin{:});
    
    figure;
    contour(X, Y, Z, 'ShowText', ip.Results.ShowText, 'LabelSpacing', ip.Results.LabelSpacing, 'LineWidth', ip.Results.LineWidth);
    axis equal;
    title(ip.Results.fname);
    xlabel('Longitude');
    ylabel('Latitude');
    ax = gca;
    ax.XAxis.TickLabelFormat = '%g°';
    ax.YAxis.TickLabelFormat = '%g°';
    cb = colorbar;
    cb.Label.String = ip.Results.cname;
end
