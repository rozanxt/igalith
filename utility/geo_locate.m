% Restrict a raster tile map according to geographic coordinates
function data = geo_locate(data, geod)
    arcm = 60;
    lonp = 180;
    latp = 90;
    data = data(arcm * (geod.lon(1) + lonp) : arcm * (geod.lon(2) + lonp) - 1, ...
                arcm * (geod.lat(1) + latp) : arcm * (geod.lat(2) + latp) - 1);
end
