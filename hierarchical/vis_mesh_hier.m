% Plot hierarchical isogeometric mesh
function vis_mesh_hier(msh_hrc)
    figure;
    hmsh_plot_cells(msh_hrc);
    title('Hierarchical isogeometric mesh');
    pbaspect([1, 1, 0.5]);
    view(0, 90);
end
