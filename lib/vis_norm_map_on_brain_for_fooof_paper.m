function f = vis_norm_map_on_brain_for_fooof_paper(norm_meas,feat_names,atlas,options)
% VIS_NORM_MAP_ON_BRAIN Plot normative map (all features) on brain surface
% in one or more views.
%
% Creates a new  figure with views x features subplots.
%
% Calls vis_plot_abnormalities_on_brain with options that work well for
% this application.
%
%   f = VIS_NORM_MAP_ON_BRAIN(norm_meas,feat_names,atlas) returns the
%   figure handle, f, to the plots of the normative map, norm_meas (ROIs x
%   features matrix). Each column of subplots corresponds to a different
%   feautre whose names are given in the cell array feat_names. ROI info is
%   passed via atlas, a table that must contain variables "names" (cell
%   arrays of ROI names) and "xyz" (ROI xyz coordinates). By default, the
%   first row of atlas is used for ROI info, but a different row can be
%   specified using the optional argument "AtlasIndex".
%
%   See argument validation comments for Name/Value pairs that provide
%   limited customisation (e.g., views to plot).
%
% See also VIS_PLOT_ABNORMALITES_ON_BRAIN.
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023

arguments
    norm_meas {mustBeNumeric}
    feat_names cell
    atlas table
    options.View cell = {'top','anterior','left','right'};  % must be subset of the default views
    options.MarkerSize (1,1) = 200                          % size to plot ROI markers
    options.Colormap = 'parula';    % colormap color
    options.ColorbarLocation = 'southoutside';              % colorbar location; see colorbar options
    options.PlotBrain = false        % whether to plot brain surface (can be useful to just export markers in epsc format - surface needs to be exported as non-vector graphic)
    options.FontSize = 32           % fontsize (colorbar labels only)
    options.TitleFontSize = 16;     % size of title (feature names)
    options.AtlasIndex = 4;
    options.MarkerEdgeOff = false ;  % whether to turn off plotting marker edges in a different colour; LineWidthSpared will still influence total marker size
    options.LineWidth = 0.5;        % marker linewidth
    options.MarkerEdgeColor = [0 0 0];
    options.climits = [0 4] % marker color if MarkerEdgeOff = false
    options.resected = [] %add resected values if needed
end


% number of views to plot brain
n_view = length(options.View);

% start figure and set size
f=figure();
set(f,'units','centimeters','position',[2 2 10+7 10+n_view*6])

% counter for subplots
my_count=1;

% plot for each view and feature
for i=1:n_view


    % which rois to plot (0 = plot, NaN = missing from map)
    x = single(isnan(norm_meas(:,:)));
    x(x==1) = NaN;
    
    % Add resected ROIs
    x(options.resected) = 1;
    % only plot colorbar for last view
    if i == n_view
        cbar = true;
    else
        cbar = false;
    end

    ax = subplot(n_view,1,my_count);
    pos = ax.Position;

    % call vis_plot_abnormalities_on_brain
    vis_plot_abnormalities_on_brain_fooof_paper(norm_meas(:,:),x,...
        atlas.name{options.AtlasIndex},atlas.xyz{options.AtlasIndex},...
        'SparedMarkerEdgeOff',false,'SparedColor',[1,1,1]',...
        'PlotNewFig',false,'MarkerSize',options.MarkerSize,'View',options.View{i},...
        'PlotColorbar',cbar,'Colormap',options.Colormap,'ColorbarLocation',...
        options.ColorbarLocation,'FontSize',options.FontSize,'PlotBrain',options.PlotBrain,...
        'SparedMarkerEdgeOff',options.MarkerEdgeOff,'LineWidthSpared',2,'CLim',options.climits);

    % return subplot to full size if has colorbar
    if i == n_view
        ax.Position = pos;
    end

    % title
    if i == 1
        title(feat_names,'FontSize',options.TitleFontSize)
    end

    % subplot counter
    my_count = my_count+1;

end
end