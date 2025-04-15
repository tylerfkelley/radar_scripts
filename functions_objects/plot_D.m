function plot_D(D,x,t,title_str)
    % create colorbar limits
    clim_abs_max = 250;
    threshold    = 0.99;
    D_sort = sort(D(:),'ascend');
    clims  = [D_sort(round(numel(D_sort)*(1-threshold))), D_sort(round(numel(D_sort)*threshold))];
    clims  = [max(-clim_abs_max, clims(1)),min(clim_abs_max, clims(2))];

    figure(); colormap(gray(2^12));
    imagesc(x, t.*1e9, D);
    colorbar;
    caxis(clims);
    title(title_str);
    xlabel("x (m)");
    ylabel("t (ns)");
end
