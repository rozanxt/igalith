%% Topography and Mohorovicic depth data in Europe

% Geographic location
geod = geo_location('europe');

% Topographic elevation
file = fopen('../data/earth/Earth2014.BED2014.1min.geod.bin');
topo = fread(file, [21600, 10800], 'int16', 'ieee-be');
topo = 1e-3 * rot90(geo_locate(topo, geod));
fclose(file);

% Topographic load
file = fopen('../data/earth/Earth2014.RET2014.1min.geod.bin');
retl = fread(file, [21600, 10800], 'int16', 'ieee-be');
retl = 1e-3 * rot90(geo_locate(retl, geod));
fclose(file);

% Mohorovicic depth
[I, J, moho] = grdread2('../data/europe/Europe_moho_depth_2007.grd');
K = find(I == geod.lon(1)) : (find(I == geod.lon(2)) - 1);
L = find(J == geod.lat(1)) : (find(J == geod.lat(2)) - 1);
moho = -flipud(moho(L, K));

% Visualization of the data
vis_raster_map(topo, geod, 'fname', 'Topographic map (European Plate)', 'cname', 'Elevation (km)', 'ShowText', false);
vis_raster_map(retl, geod, 'fname', 'Topographic load (European Plate)', 'cname', 'Rock-equivalent topography (km)', 'ShowText', false);
vis_raster_map(moho, geod, 'fname', 'Moho depth map (European Plate)', 'cname', 'Elevation (km)', 'ShowText', true);
