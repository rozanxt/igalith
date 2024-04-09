% Visualize results of the Airy-Heiskanen model of local isostasy
function vis_airy_heiskanen(prb, iso, typ)
    [X, Y] = meshgrid(linspace(0, prb.dmsz(1), prb.mpsz(1)), linspace(0, prb.dmsz(2), prb.mpsz(2)));
    [lon, lat] = geo_coords(X, Y, prb.geod, prb.dmsz);

    switch typ
        case 'topo'
            vis_contour_map(lon, lat, iso.tp(X, Y), 'fname', 'Topographic load', 'cname', 'Rock-equivalent topography (km)', 'ShowText', false);
        case 'crth'
            vis_contour_map(lon, lat, iso.ct(X, Y), 'fname', 'Crustal thickness (Airy–Heiskanen)', 'cname', 'Thickness (km)', 'ShowText', false);
        case 'moho'
            vis_contour_map(lon, lat, iso.mh(X, Y), 'fname', 'Lithospheric depression (Airy–Heiskanen)', 'cname', 'Elevation (km)', 'ShowText', false);
    end
end
