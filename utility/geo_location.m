% Geographic coordinates of locations of interest
function geod = geo_location(name)
    switch name
        case 'europe'
            geod.lon = [-25, 25];
            geod.lat = [28, 78];
        case 'himalaya'
            geod.lon = [60, 120];
            geod.lat = [20, 50];
        case 'hawaii'
            geod.lon = [-165, -150];
            geod.lat = [13, 28];
        case 'indonesia'
            geod.lon = [90, 150];
            geod.lat = [-15, 15];
        case 'java_island'
            geod.lon = [105, 115];
            geod.lat = [-10, -5];
        case 'central_java'
            geod.lon = [109.5, 111.75];
            geod.lat = [-8.5, -6.25];
        case 'bandung_basin'
            geod.lon = [107.25, 108];
            geod.lat = [-7.25, -6.5];
        case 'bandung'
            geod.lon = [107.5, 107.75];
            geod.lat = [-7, -6.75];
    end
    geod.sz = [geod.lon(2) - geod.lon(1), geod.lat(2) - geod.lat(1)];
end
