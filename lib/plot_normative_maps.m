
switch type_power_spectrum
    % for Complete PSD

    case 'complete'

        %climits
        climits=[nanmean(normative_maps,1)'-(nanstd(normative_maps,1))' nanmean(normative_maps,1)'+(nanstd(normative_maps,1))'];%[0.155 0.19];%Delta

        %feature names
        feat_names={'Delta Bandpower','Theta Bandpower','Alpha Bandpower','Beta Bandpower'};

        for i=1:size(feat_names,2)
            bp_feat=feat_names{i};

            h=vis_norm_map_on_brain_for_fooof_paper(normative_maps(:,i),feat_names(i),atlas_tbl,'AtlasIndex',parc,'MarkerSize',200,'Colormap',sort(colormap_idx{i},'descend'),...
                'View',{'left','top','right'},'climits',climits(i,:),'PlotBrain',true,'FontSize',20);

            %Save
            %set(gcf,'renderer','painters');

            saveas(h, [analysis_location,figdir_1, bp_feat,'_parc_',num2str(parc),'.pdf']); % specify the full path
            close all
        end

    case 'periodic'
        % PLOT flat PSD
        %climits
        climits=[nanmean(normative_maps,1)'-(nanstd(normative_maps,1))' nanmean(normative_maps,1)'+(nanstd(normative_maps,1))'];%[0.155 0.19];%Delta

        %feature names
        feat_names={'Periodic Delta Bandpower','Periodic Theta Bandpower','Periodic Alpha Bandpower','Periodic Beta Bandpower'};

        for i=1:size(feat_names,2)
            bp_feat=feat_names{i};

            h=vis_norm_map_on_brain_for_fooof_paper(normative_maps(:,i),feat_names(i),atlas_tbl,'AtlasIndex',parc,'MarkerSize',200,'Colormap',sort(colormap_idx{i},'descend'),...
                'View',{'left','top','right'},'climits',climits(i,:),'PlotBrain',true,'FontSize',20);

            %Save
            %set(gcf,'renderer','painters');
            saveas(h, [analysis_location,figdir_1,bp_feat,'_parc_',num2str(parc),'.pdf']); % specify the full path
            close all
        end

    case 'aperiodic'
        climits=[nanmean(normative_maps,1)-(nanstd(normative_maps,1)) nanmean(normative_maps,1)+(nanstd(normative_maps,1))];

        feat_names={'Aperiodic exponent'};
        h=vis_norm_map_on_brain_for_fooof_paper(normative_maps,feat_names,atlas_tbl,'AtlasIndex',parc,'MarkerSize',200,'Colormap',sort(colormap_idx{5},'descend'),...
            'View',{'left','top','right'},'climits',climits,'PlotBrain',true,'FontSize',20);

        %Save
        %set(gcf,'renderer','painters');
        saveas(h, [analysis_location,figdir_1, feat_names{size(feat_names)},'_parc_',num2str(parc),'.pdf']); % specify the full path
        close all
end