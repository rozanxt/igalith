% Convert simulation coordinates to geographic coordinates
function [lon, lat] = geo_coords(X, Y, geod, dmsz)
    [lon, lat] = deal(geod.lon(1) + (geod.sz(1) / dmsz(1)) * X, geod.lat(1) + (geod.sz(2) / dmsz(2)) * Y);
end
