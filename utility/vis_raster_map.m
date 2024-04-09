% Visualize raster data on a grid map (using a contour plot)
function vis_raster_map(data, geod, varargin)
    mpsz = fliplr(size(data));
    grsz = [min(size(data)), min(size(data))];
    dmsz = mpsz ./ grsz;
    
    f = @(x, y) eval_uv(data, x, y, grsz);
    [X, Y] = meshgrid(linspace(0, dmsz(1), mpsz(1)), linspace(0, dmsz(2), mpsz(2)));
    [lon, lat] = geo_coords(X, Y, geod, dmsz);
    
    vis_contour_map(lon, lat, f(X, Y), varargin{:});
end
