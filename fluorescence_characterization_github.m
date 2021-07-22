close all

%% 
segment_folder_info = dir("/Users/ctzouanas/Documents/MIT/Shalek/Spacecat/190720 - Full Timecourse Segmentation");
cell_count = zeros(1,10);
% At 40x magnification, 2048 pixels across (i.e. the dimension of 1 FOV)
% equals 333 um.
um_per_pixel = 333/2048;
% Radius of the region being photoactivated. One image had a slightly wider
% radius to enable sufficient cells to be in the region-of-interest
radius_value = um_per_pixel * [192 192 250 192 192 192 192 192 192 192];
% radius_value = zeros(1,10) + 192;
segment_cell = cell(1,10);
radius_inside = [];
dF_F0_inside = [];
radius_outside = [];
dF_F0_outside = [];

for ii = 4:size(segment_folder_info,1)
    segment_string_ii = strcat(segment_folder_info(ii).folder, '/', segment_folder_info(ii).name);
    analysis_spreadsheet_ii = readmatrix(segment_string_ii, 'NumHeaderLines', 1);
    % Use the below line to pick out the mean intensity for a cell
%     mean_fluor_mat = analysis_spreadsheet_ii(:,21:26);
    % Use the below line to pick out the integrated intensity for a cell
    mean_fluor_mat = analysis_spreadsheet_ii(:,3:8);
    x_locat_ii = analysis_spreadsheet_ii(:,51);
    y_locat_ii = analysis_spreadsheet_ii(:,52);
    radius_ii = um_per_pixel*sqrt((x_locat_ii-1024).^2 + (y_locat_ii-1024).^2);
    cell_count(ii-3) = numel(radius_ii(radius_ii < radius_value(ii-3)));
    repeated_init = repmat(mean_fluor_mat(:,1),1,6);
    % % Use the below line to calculate dF/F0
        dF_F0_mat_ii = (mean_fluor_mat - repeated_init)./(repeated_init);
    % % Use the below line to calculate change in mean or integrated,
    % % depending on which columns were picked out above
%         dF_F0_mat_ii = (mean_fluor_mat - repeated_init);
    % % Use the below line to just plot the mean or integrated fluorescence,
    % % depending on which columns were picked out above
%     dF_F0_mat_ii = (mean_fluor_mat);
    inside_radius = radius_ii < radius_value(ii-3);
    outside_radius = radius_ii > radius_value(ii-3);
    dF_F0_mat_ii_inside = dF_F0_mat_ii(inside_radius,:);
    dF_F0_mat_ii_outside = dF_F0_mat_ii(outside_radius,:);

    
    radius_inside = [radius_inside; radius_ii(radius_ii < radius_value(ii-3))];
    radius_outside = [radius_outside; radius_ii(radius_ii > radius_value(ii-3))];
    % % Use the below two lines if trying to plot actual mean or integrated
    % % fluorescence (i.e. not the change)
    dF_F0_inside = [dF_F0_inside; dF_F0_mat_ii_inside(:,1:end)];
    dF_F0_outside = [dF_F0_outside; dF_F0_mat_ii_outside(:,1:end)];    
    %     % Use the below two lines if trying to plot change in fluorescence
    %     dF_F0_inside = [dF_F0_inside; dF_F0_mat_ii_inside(:,2:end)];
    %     dF_F0_outside = [dF_F0_outside; dF_F0_mat_ii_outside(:,2:end)];
    
end

% Uncomment the below three lines if doing a difference
dF_F0_inside = dF_F0_inside(:,2:end);
dF_F0_outside = dF_F0_outside(:,2:end);
timepoint_order = [1 2 5 3 4];
% random_indices = randi([1 numel(dF_F0_outside)], [1, numel(dF_F0_inside)]);
% random_indices = 1:numel(dF_F0_outside);

y_shift_vec = linspace(0, 1.5*(numel(timepoint_order)-1), (numel(timepoint_order)));
figure
hold on
set(gca, 'FontSize', 14)
color_vec_in = [0 0.8 0.26]*0.8;
color_vec_out = [0.35 0.35 0.35];
for ii = 1:numel(timepoint_order)
    index_ii = timepoint_order(ii);
    [ks_outside, xi_outside, bw_outside] = ksdensity(dF_F0_outside(:,index_ii));
    [ks_inside, xi_inside, bw_inside] = ksdensity(dF_F0_inside(:,index_ii), 'Bandwidth', 0.5);
    if min(xi_outside) > -0.1
        xi_outside = [-0.1 xi_outside];
        ks_outside = [0 ks_outside];
    end
    if min(xi_inside) > -0.1
        xi_inside = [-0.1 xi_inside];
        ks_inside = [0 ks_inside];
    end
    if max(xi_outside) < 7
        xi_outside = [xi_outside 7];
        ks_outside = [ks_outside 0];
    end
    if max(xi_inside) < 7
        xi_inside = [xi_inside 7];
        ks_inside = [ks_inside 0];
    end
    
    plot(xi_inside, ks_inside/max(ks_inside) + y_shift_vec(ii), 'Linewidth', 3, 'Color', color_vec_in)
    plot(xi_outside, ks_outside/max(ks_outside) + y_shift_vec(ii), 'Linewidth', 3, 'Color', color_vec_out)
    
    disp(timepoint_order(ii))
    [h,p] = ttest2(dF_F0_outside(:,index_ii), dF_F0_inside(:,index_ii))
end
xlim([-0.1 7])
ylabel('Density (A.U.)')
yticks(0.5 + y_shift_vec)
yticklabels({'0''' '10''' '60''' '120''' '180'''})
xlabel('\DeltaF/F_0')
set(gca, 'TickLength',[0 0])
set(gca, 'FontName', 'Myriad Pro')
set(gca, 'FontSize', 16)


% % legend(axes_dFF0_inside, '0^-', '0^+', '10', '60', '120', '180')
% % legend(axes_dFF0_outside, '0^-', '0^+', '10', '60', '120', '180')
% legend(axes_dFF0_inside, '0^+ min', '10 min', '60 min', '120 min', '180 min')

% figure
% hold on
% scatter(radius_outside, dF_F0_outside(:,5), 30, [0 0 0], 'filled');
% alpha 0.3
% scatter(radius_inside, dF_F0_inside(:,5), 30, [0 0.75 0.25], 'filled');
% xlabel('Distance from Center of Excitation Region (pixels)')
% ylabel('Change in Fluorescence (\DeltaF/F_0)')
% set(gca, 'FontSize', 16)
% legend('Unphotoactivated', 'Photoactivated')
% legend boxoff

% nbins_outside = 25;
% combined_outside_dFF0_mat_180 = [radius_outside dF_F0_outside(:,5)];
% [n_outside,c_outside] = hist3(combined_outside_dFF0_mat_180,'Edges',{linspace(0, 1300*um_per_pixel, nbins_outside) linspace(0, 9, nbins_outside)});
% [xold_outside, yold_outside] = meshgrid(c_outside{1}, c_outside{2});
% meshgrid_x_outside = linspace(min(c_outside{1}), max(c_outside{1}), 600);
% meshgrid_y_outside = linspace(min(c_outside{2}),max(c_outside{2}),600);
% [xnew_outside, ynew_outside] = meshgrid(meshgrid_x_outside, meshgrid_y_outside);
% n_interp_outside = interp2(xold_outside,yold_outside,(n_outside/max(max(n_outside))),xnew_outside,ynew_outside, 'makima');
% 
% nbins_inside = 16;
% combined_inside_dFF0_mat_180 = [radius_inside dF_F0_inside(:,5)];
% [n_inside,c_inside] = hist3(combined_inside_dFF0_mat_180,'Edges',{linspace(0, 1300*um_per_pixel, nbins_inside) linspace(0, 9, nbins_inside)});
% [xold_inside, yold_inside] = meshgrid(c_inside{1}, c_inside{2});
% meshgrid_x_inside = linspace(min(c_inside{1}), max(c_inside{1}), 600);
% meshgrid_y_inside = linspace(min(c_inside{2}),max(c_inside{2}),600);
% [xnew_inside, ynew_inside] = meshgrid(meshgrid_x_inside, meshgrid_y_inside);
% n_interp_inside = interp2(xold_inside,yold_inside,(n_inside/max(max(n_inside))),xnew_inside,ynew_inside, 'makima');
% 
% number_contour_levels = 7;
% contour_levels_vec = (1/number_contour_levels):(1/number_contour_levels):1;

% figure
% hold on
% % scatter(radius_inside, dF_F0_inside(:,5), 30, [0 0.3 0.1], 'filled');
% % scatter(radius_outside, dF_F0_outside(:,5), 30, [0 0 0], 'filled');
% % alpha 0.7
% [contour_matrix_in, contour_object_in] = contour(xnew_inside,ynew_inside,(n_interp_inside)', contour_levels_vec, 'Linewidth', 2, 'Showtext', 'Off');
% [contour_matrix_out, contour_object_out] = contour(xnew_outside,ynew_outside,(n_interp_outside)', contour_levels_vec, 'Linewidth', 2, 'Showtext', 'Off');
% contour_object_in.LineColor = [0 0.5 0.16];
% contour_object_out.LineColor = [0 0 0];
% % clabel(contour_matrix_in, contour_object_in, 'FontSize', 12)
% xticks([0 50 100 150 200])
% xticklabels({'0' '50' '100' '150' '200'})
% yticks([0 1 2 3 4 5])
% yticklabels({'0' '1' '2' '3' '4' '5'})
% legend('Photoactivated', 'Unhotoactivated')
% set(gca, 'Fontsize', 16)
% ylabel('\DeltaF/F_0')
% xlabel('Distance from Center of Photoactivation Region (\mum)')
% ylim([0 5])
% set(gca, 'FontName', 'Myriad Pro')

% % Add a little padding so that both the inside and outside isolated regions
% % have the same domain
% xi_outside_2 = [xi_outside xi_inside(end)];
% ks_outside_2 = [ks_outside ks_outside(2)];
% xi_inside_2 = [xi_outside(1) xi_inside];
% ks_inside_2 = [ks_outside(1) ks_inside];
% 
% ks_outside_integ = 1 - cumsum(ks_outside_2)*(xi_outside(2)-xi_outside(1));
% ks_inside_integ = 1 - cumsum(ks_inside_2)*(xi_inside_2(3)-xi_inside_2(2));
% 
% distribution_domain = linspace(xi_outside_2(1), xi_outside_2(end));
% pchip_outside = pchip(xi_outside_2, ks_outside_integ, distribution_domain);
% pchip_inside = pchip(xi_inside_2, ks_inside_integ, distribution_domain);
% isolation_probability = (pchip_inside) ./ (pchip_outside + pchip_inside);

% figure
% hold on
% plot(xi_outside_2, ks_outside_2, 'Linewidth', 3, 'Color', 'k')
% plot(xi_inside_2, ks_inside_2, 'Linewidth', 3, 'Color', [0 0.75 0.25])
% set(gca, 'Fontsize', 16)
% legend('Unphotoactivated', 'Photoactivated')
% xlabel('Change in Fluorescence (\DeltaF/F_0)')
% legend boxoff
% ylim([0 1.3])

% figure
% hold on
% plot(xi_outside_2, ks_outside_integ)
% plot(xi_inside_2, ks_inside_integ)
% legend('Outside Cumulative Distribution', 'Inside Cumulative Distribution')
% xlabel('Change in Fluorescence (dF/F0)')
% ylabel('Proportion of Cells with Fluorescence Change > X')

% figure
% plot(distribution_domain, isolation_probability, 'Linewidth', 3, 'Color', [0 0.75 0.25])
% xlabel('Change in Fluorescence (\DeltaF/F_0)')
% ylabel('Isolation Specificity')
% xlim([-Inf 5])
% set(gca, 'FontSize', 16)

