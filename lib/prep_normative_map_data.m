
switch type_power_spectrum
    case 'complete'

        %Complete PSD

        normative_maps=nan(max(rois),n_bands);
        for i=1:max(rois)
            locs=find(rois(:)==i);
            normative_maps(i,:)=mean(rel_bp_complete(locs,:),1);

        end

    case 'periodic'
        % flat PSD

        normative_maps=nan(max(rois),n_bands);
        for i=1:max(rois)
            locs=find(rois(:)==i);
            normative_maps(i,:)=mean(rel_bp_periodic(locs,:),1);

        end
    case 'aperiodic'
        % Aperiodic
        normative_maps=nan(max(rois),1);
        for i=1:max(rois)
            locs=find(rois(:)==i);
            normative_maps(i)=mean(normative_table.aperiodic_cmps_2(locs));

        end

end
