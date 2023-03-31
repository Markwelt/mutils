classdef raincloud_plot < handle
    % raincloud_plot creates raincloud plots for some data
    %   A raincloud plot is an easy to read substitute for 
    %   a box plot that replaces the box shape with a kernel
    %   density estimate of the data (cloud), and optionally 
    %   adjoins the data points (rain).
    %
    %   Additional constructor parameters include the width 
    %   of the plot, the bandwidth of the kernel density 
    %   estimation, and the X-axis position of the raincloud plot.
    %
    %   Use raincloud_plot for a <a href="matlab:help('boxplot')">boxplot</a>-like wrapper for
    %   interactive plotting.
    %
    % raincloud_plot Properties:
    %    raincloudColor - Fill color of the raincloud area and data points.
    %                  Defaults to the next default color cycle.
    %    raincloudAlpha - Transparency of the ciolin area and data points.
    %                  Defaults to 0.3.
    %    EdgeColor   - Color of the raincloud area outline.
    %                  Defaults to [0.5 0.5 0.5]
    %    BoxColor    - Color of the box, whiskers, and the outlines of
    %                  the median point and the notch indicators.
    %                  Defaults to [0.5 0.5 0.5]
    %    MedianColor - Fill color of the median and notch indicators.
    %                  Defaults to [1 1 1]
    %    ShowData    - Whether to show data points.
    %                  Defaults to true
    %    ShowNotches - Whether to show notch indicators.
    %                  Defaults to false
    %    ShowMean    - Whether to show mean indicator.
    %                  Defaults to false
    %
    % raincloud_plot Children:
    %    ScatterPlot - <a href="matlab:help('scatter')">scatter</a> plot of the data points
    %    raincloudPlot  - <a href="matlab:help('fill')">fill</a> plot of the kernel density estimate
    %    BoxPlot     - <a href="matlab:help('fill')">fill</a> plot of the box between the quartiles
    %    WhiskerPlot - line <a href="matlab:help('plot')">plot</a> between the whisker ends
    %    MedianPlot  - <a href="matlab:help('scatter')">scatter</a> plot of the median (one point)
    %    NotchPlots  - <a href="matlab:help('scatter')">scatter</a> plots for the notch indicators
    %    MeanPlot    - line <a href="matlab:help('plot')">plot</a> at mean value
    %
    % Example usage:
    %   rcp = raincloud_plot(randi(1000,[1, 500]),1,'EdgeColor',[0.5 0 0.2]);
    %   rcpSizes = ones([1, 500])*20;
    %   rcp.ScatterPlot.SizeDataSource = 'rcpSizes';
    %   refreshdata(rcp.ScatterPlot) %  bigger scatter points
    %   rcp.WhiskerPlot.Visible = 'off';
    %

    % Adapted by MaFu from violin_plot.m (Matlab Exchange) Violin.m (from github: Violinplot-Matlab)
    % Copyright (c) 2016, Bastian Bechtold
    % This code is released under the terms of the BSD 3-clause license

    properties
        ScatterPlot % scatter plot of the data points
        raincloudPlot  % fill plot of the kernel density estimate
        BoxPlot     % fill plot of the box between the quartiles
        WhiskerPlot % line plot between the whisker ends
        MedianPlot  % scatter plot of the median (one point)
        NotchPlots  % scatter plots for the notch indicators
        MeanPlot    % line plot of the mean (horizontal line)
    end

    properties (Dependent=true)
        raincloudColor % fill color of the raincloud area and data points
        raincloudAlpha % transparency of the raincloud area and data points
        EdgeColor   % color of the raincloud area outline
        BoxColor    % color of box, whiskers, and median/notch edges
        MedianColor % fill color of median and notches
        ShowData    % whether to show data points
        ShowNotches % whether to show notch indicators
        ShowMean    % whether to show mean indicator
    end

    methods
        function obj = raincloud_plot(data, pos, varargin)
            %raincloud_plot plots a raincloud plot of some data at pos
            %   raincloud(DATA, POS) plots a raincloud at x-position POS for
            %   a vector of DATA points.
            %
            %   raincloud(..., 'PARAM1', val1, 'PARAM2', val2, ...)
            %   specifies optional name/value pairs:
            %     'Width'        Width of the raincloud in axis space.
            %                    Defaults to 0.3
            %     'Bandwidth'    Bandwidth of the kernel density
            %                    estimate. Should be between 10% and
            %                    40% of the data range.
            %     'raincloudColor'  Fill color of the raincloud area and
            %                    data points. Defaults to the next
            %                    default color cycle.
            %     'raincloudAlpha'  Transparency of the raincloud area and
            %                    data points. Defaults to 0.3.
            %     'EdgeColor'    Color of the raincloud area outline.
            %                    Defaults to [0.5 0.5 0.5]
            %     'BoxColor'     Color of the box, whiskers, and the
            %                    outlines of the median point and the
            %                    notch indicators. Defaults to
            %                    [0.5 0.5 0.5]
            %     'MedianColor'  Fill color of the median and notch
            %                    indicators. Defaults to [1 1 1]
            %     'ShowData'     Whether to show data points.
            %                    Defaults to true
            %     'ShowNotches'  Whether to show notch indicators.
            %                    Defaults to false
            %     'ShowMean'     Whether to show mean indicator.
            %                    Defaults to false

            args = obj.checkInputs(data, pos, varargin{:});
            data = data(not(isnan(data)));
            if numel(data) == 1
                obj.MedianPlot = scatter(pos, data, 'filled');
                obj.MedianColor = args.MedianColor;
                obj.MedianPlot.MarkerEdgeColor = args.EdgeColor;
                return
            end

            hold('on');

            % calculate kernel density estimation for the raincloud
            [density, value] = ksdensity(data, 'bandwidth', args.Bandwidth);
            density = density(value >= min(data) & value <= max(data));
            value = value(value >= min(data) & value <= max(data));
            value(1) = min(data);
            value(end) = max(data);
            ndens = numel(value);

            % all data is identical
            if min(data) == max(data)
                density = 1;
            end

            if isempty(args.Width)
                width = 0.3/max(density);
            else
                width = args.Width/max(density);
            end
            widthdens = width*max(density);

            % plot the data points within the raincloud area
            jitter =  0.5*widthdens.*(rand(size(data))-0.5);  %vi0lin had instead 2*(rand()) and jitter.*jitterstrength
            obj.ScatterPlot = ...
                scatter(pos - 0.4*widthdens + jitter, data, 20, 'filled'); %MaFu addition

            % plot the data mean
            meanValue = mean(value);
            if length(density) > 1
                meanDensity = interp1(value, density, meanValue);
            else % all data is identical:
                meanDensity = density;
            end
            obj.MeanPlot = plot([pos pos+meanDensity*width], ...
                                [meanValue meanValue]);
            obj.MeanPlot.LineWidth = 0.75;

            % plot the raincloud
            obj.raincloudPlot =  ... % plot color will be overwritten later
                fill([pos+density*width ones([1 ndens])*pos], ... violin had instead pos-density(end:-1:1)*width
                     [value value(end:-1:1)], [1 1 1], 'Marker', 'none', 'LineStyle', '-'); %MaFu addition

            % plot the mini-boxplot within the raincloud
            quartiles = quantile(data, [0.25, 0.5, 0.75]);
            obj.BoxPlot = ... % plot color will be overwritten later
                fill([pos-0.005 pos+0.005 pos+0.005 pos-0.005], ...
                     [quartiles(1) quartiles(1) quartiles(3) quartiles(3)], ...
                     [1 1 1], 'Marker', 'none', 'LineStyle', '-'); %MaFu addition
            IQR = quartiles(3) - quartiles(1);
            lowhisker = quartiles(1) - 1.5*IQR;
            lowhisker = max(lowhisker, min(data(data > lowhisker)));
            hiwhisker = quartiles(3) + 1.5*IQR;
            hiwhisker = min(hiwhisker, max(data(data < hiwhisker)));
            if ~isempty(lowhisker) && ~isempty(hiwhisker)
                obj.WhiskerPlot = plot([pos pos], [lowhisker hiwhisker], 'Marker', 'none', 'LineStyle', '-'); %MaFu addition
            end
            obj.MedianPlot = scatter(pos, quartiles(2), [], [1 1 1], 'filled');

            obj.NotchPlots = ...
                 scatter(pos, quartiles(2)-1.57*IQR/sqrt(length(data)), ...
                         [], [1 1 1], 'filled', '^');
            obj.NotchPlots(2) = ...
                 scatter(pos, quartiles(2)+1.57*IQR/sqrt(length(data)), ...
                         [], [1 1 1], 'filled', 'v');

            obj.EdgeColor = args.EdgeColor;
            obj.BoxColor = args.BoxColor;
            obj.MedianColor = args.MedianColor;
            if not(isempty(args.raincloudColor))
                obj.raincloudColor = args.raincloudColor;
            else
                obj.raincloudColor = obj.ScatterPlot.CData;
            end
            obj.raincloudAlpha = args.raincloudAlpha;
            obj.ShowData = args.ShowData;
            obj.ShowNotches = args.ShowNotches;
            obj.ShowMean = args.ShowMean;
        end

        function set.EdgeColor(obj, color)
            obj.raincloudPlot.EdgeColor = color;
        end

        function color = get.EdgeColor(obj)
            color = obj.raincloudPlot.EdgeColor;
        end

        function set.MedianColor(obj, color)
            obj.MedianPlot.MarkerFaceColor = color;
            if not(isempty(obj.NotchPlots))
                obj.NotchPlots(1).MarkerFaceColor = color;
                obj.NotchPlots(2).MarkerFaceColor = color;
            end
        end

        function color = get.MedianColor(obj)
            color = obj.MedianPlot.MarkerFaceColor;
        end

        function set.BoxColor(obj, color)
            obj.BoxPlot.FaceColor = color;
            obj.BoxPlot.EdgeColor = color;
            obj.WhiskerPlot.Color = color;
            obj.MedianPlot.MarkerEdgeColor = color;
            obj.NotchPlots(1).MarkerEdgeColor = color;
            obj.NotchPlots(2).MarkerEdgeColor = color;
        end

        function color = get.BoxColor(obj)
            color = obj.BoxPlot.FaceColor;
        end

        function set.raincloudColor(obj, color)
            obj.raincloudPlot.FaceColor = color;
            obj.ScatterPlot.MarkerFaceColor = color*0.7; %MaFu addition
            obj.MeanPlot.Color = color;
        end

        function color = get.raincloudColor(obj)
            color = obj.raincloudPlot.FaceColor;
        end

        function set.raincloudAlpha(obj, alpha)
            obj.ScatterPlot.MarkerFaceAlpha = alpha;
            obj.raincloudPlot.FaceAlpha = alpha;
        end

        function alpha = get.raincloudAlpha(obj)
            alpha = obj.raincloudPlot.FaceAlpha;
        end

        function set.ShowData(obj, yesno)
            if yesno
                obj.ScatterPlot.Visible = 'on';
            else
                obj.ScatterPlot.Visible = 'off';
            end
        end

        function yesno = get.ShowData(obj)
            yesno = logical(strcmp(obj.NotchPlots(1).Visible, 'on'));
        end

        function set.ShowNotches(obj, yesno)
            if yesno
                obj.NotchPlots(1).Visible = 'on';
                obj.NotchPlots(2).Visible = 'on';
            else
                obj.NotchPlots(1).Visible = 'off';
                obj.NotchPlots(2).Visible = 'off';
            end
        end

        function yesno = get.ShowNotches(obj)
            yesno = logical(strcmp(obj.ScatterPlot.Visible, 'on'));
        end

        function set.ShowMean(obj, yesno)
            if yesno
                obj.MeanPlot.Visible = 'on';
            else
                obj.MeanPlot.Visible = 'off';
            end
        end

        function yesno = get.ShowMean(obj)
            yesno = logical(strcmp(obj.MeanPlot.Visible, 'on'));
        end
    end

    methods (Access=private)
        function results = checkInputs(~, data, pos, varargin)
            isscalarnumber = @(x) (isnumeric(x) & isscalar(x));
            p = inputParser();
            p.addRequired('Data', @isnumeric);
            p.addRequired('Pos', isscalarnumber);
            p.addParameter('Width', [], isscalarnumber);
            p.addParameter('Bandwidth', [], isscalarnumber);
            iscolor = @(x) (isnumeric(x) & length(x) == 3);
            p.addParameter('raincloudColor', [], iscolor);
            p.addParameter('BoxColor', [0.5 0.5 0.5], iscolor);
            p.addParameter('EdgeColor', [0.5 0.5 0.5], iscolor);
            p.addParameter('MedianColor', [1 1 1], iscolor);
            p.addParameter('raincloudAlpha', 0.3, isscalarnumber);
            isscalarlogical = @(x) (islogical(x) & isscalar(x));
            p.addParameter('ShowData', true, isscalarlogical);
            p.addParameter('ShowNotches', false, isscalarlogical);
            p.addParameter('ShowMean', false, isscalarlogical);

            p.parse(data, pos, varargin{:});
            results = p.Results;
        end
    end
end
