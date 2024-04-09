% Output function for each step of the optimization procedure
function output_fcn(dsn, obj, dsc, nrm, stp)
    fprintf('function value: %f, gradient norm: %f, step size: %f\n', obj, nrm, stp);
end
