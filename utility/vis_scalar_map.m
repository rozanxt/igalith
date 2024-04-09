% Visualize scalar data on a grid map (using a contour plot)
function vis_scalar_map(data, geod, lrbt, npts, varargin)
    [X, Y] = meshgrid(linspace(lrbt(1), lrbt(2), npts(1)), linspace(lrbt(3), lrbt(4), npts(2)));
    [lon, lat] = geo_coords(X, Y, geod, [lrbt(2) - lrbt(1), lrbt(4) - lrbt(3)]);
    vis_contour_map(lon, lat, data(X, Y), varargin{:});
end
