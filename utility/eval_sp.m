% Evaluate a spline function with given degrees of freedom
function [val, grd] = eval_sp(dof, spc, geo, pts, varargin)
    ptc = cell(1, size(pts, 2));
    ind = cell(1, size(pts, 2));
    for j = 1 : size(pts, 2)
        [ptc{j}, ~, ind{j}] = uniquetol(pts(:, j));
    end
    [val, grd] = sp_eval(dof, spc, geo, ptc, varargin{:});
    isz = size(val);
    cs1 = repelem({':'}, ndims(val) - 2);
    cs2 = repelem({':'}, ndims(grd) - 2);
    val = val(cs1{:}, sub2ind(isz(end - 1 : end), ind{:}));
    grd = grd(cs2{:}, sub2ind(isz(end - 1 : end), ind{:}));
end
