%%% 
% Example file (with functions included) to calculate climatology, 10% &
% 90% thresholds, and identify periods of extreme anomalies. In this
% example, temperature and salinity are the properties of interest.
% However, you can use it on any type of numerical data. 

% Author: Conor Padmanabhan, 2024
% Based on detect.m by Zhao and Marin (2019)

% For questions, email me at conorpad@gmail.com

% Stay curious!
%%%
load('/Users/conor/OneDrive - Bowdoin College/WHOI2024/Datasets/GoM Mooring Data/M01/filtered_data_and_anomaly_data.mat');

% Example datatable 
%    - Can have multiple measurements per day
%    - Lat and Lon are not required, but column names must match exactly

% FilteredDataTable:
%     lat     Lon  Depth        Date          Temperature           Salinity
%    43.49	-67.88	100	    '09-Jul-2003'	5.04530000686646	33.1260871887207
%    43.49	-67.88	150	    '09-Jul-2003'	5.73729991912842	33.5675926208496
%    43.49	-67.88	200	    '09-Jul-2003'	7.23169994354248	34.1415977478027
%    43.49	-67.88	20	    '09-Jul-2003'	8.29139995574951	32.4734497070313
%    43.49	-67.88	50	    '09-Jul-2003'	4.73479986190796	32.6783142089844
%    43.48	-67.88	250	    '09-Jul-2003'	7.28999996185303	34.1770324707031
%    43.49	-67.88	1	    '09-Jul-2003'	16.0174007415772	32.2498474121094
%    43.49	-67.88	100	    '09-Jul-2003'	5.07420015335083	33.1242866516113
%    43.49	-67.88	150	    '09-Jul-2003'	5.69659996032715	33.5404090881348
%    43.49	-67.88	200	    '09-Jul-2003'	7.05859994888306	34.0564575195313
%    43.49	-67.88	20	    '09-Jul-2003'	8.31690025329590	32.4718742370606
%    43.49	-67.88	50	    '09-Jul-2003'	4.77110004425049	32.6658210754395
%    43.49	-67.88	250	    '09-Jul-2003'	7.28000020980835	34.1705589294434
%    43.49	-67.88	1	    '09-Jul-2003'	15.8566999435425	32.2484359741211
%    43.49	-67.88	1	    '10-Jul-2003'	15.6904001235962	32.2936134338379
%    43.49	-67.88	100	    '10-Jul-2003'	4.95839977264404	33.1168975830078
%    43.49	-67.88	150	    '10-Jul-2003'	6.19180011749268	33.6583518981934
%    43.49	-67.88	200	    '10-Jul-2003'	7.21840000152588	34.1324653625488
%    43.49	-67.88	20	    '10-Jul-2003'	8.67109966278076	32.4791221618652
%    43.49	-67.88	50	    '10-Jul-2003'	4.80089998245239	32.8365020751953
%    43.48	-67.88	250 	'10-Jul-2003'	7.30999994277954	34.2019309997559


% Example output:
%    Onset          End        Dur.      Max Int.           Mean Int.            Int. Var.          Cum. Int.     Depth
% '01-Jan-2007'	'18-Jan-2007'	18	1.72306009686281	1.58586823818906	0.0667496518191561	28.5456282874031	1

% Define the valid depth index (if you want to exclude any depths for any reason)
validDepthIdx = filteredDataTable.depth ~= 10;
filteredDataTable = filteredDataTable(validDepthIdx, :);

% Filter out bad salinity data based on the thresholds
salinityThreshold200m = 32.8;
salinityThreshold250m = 32.8;
filteredDataTable.salinity(filteredDataTable.depth == 200 & filteredDataTable.salinity < salinityThreshold200m) = NaN;
filteredDataTable.salinity(filteredDataTable.depth == 250 & filteredDataTable.salinity < salinityThreshold250m) = NaN;

% Given the format of data, perform neccecary conversions to convert to datetime
% Convert time from modified Julian days to datetime
timeOrigin = datetime(1858, 11, 17, 0, 0, 0);
filteredDataTable.time = timeOrigin + days(filteredDataTable.time);
filteredDataTable.time = dateshift(filteredDataTable.time, 'start', 'day');

% Ensure cli_start, cli_end, mhw_start, and mhw_end only consider the date part
cli_start = dateshift(filteredDataTable.time(1), 'start', 'day');
cli_end = dateshift(filteredDataTable.time(end), 'start', 'day');
mhw_start = cli_start;
mhw_end = cli_end;

% Call the function to detect depth-specific extremes for temperature
[extremes_temp_table, mclim_temp, positive_threshold_temp, negative_threshold_temp] = detect_extremes(filteredDataTable.temperature, filteredDataTable.time, filteredDataTable.depth, cli_start, cli_end, mhw_start, mhw_end);

% Call the function to detect depth-specific extremes for salinity
[extremes_sal_table, mclim_sal, positive_threshold_sal, negative_threshold_sal] = detect_extremes(filteredDataTable.salinity, filteredDataTable.time, filteredDataTable.depth, cli_start, cli_end, mhw_start, mhw_end);

% Display the results
disp('Temperature Extremes:');
disp(extremes_temp_table);

disp('Salinity Extremes:');
disp(extremes_sal_table);




%% Functions
% This function calculates the daily mean of the input data, computes climatology and thresholds,
% and detects events (both positive and negative anomalies) based on the thresholds.
% Inputs:
%   - data: time series data (temperature or salinity)
%   - time: corresponding time series
%   - depths: depth values corresponding to the data
%   - cli_start: start date for climatology calculation
%   - cli_end: end date for climatology calculation
%   - mhw_start: start date for event detection
%   - mhw_end: end date for event detection
% Outputs:
%   - events_table: table of detected events with metrics
%   - mclim: climatology values
%   - positive_threshold: threshold for positive anomalies
%   - negative_threshold: threshold for negative anomalies
function [events_table, mclim, positive_threshold, negative_threshold] = detect_extremes(data, time, depths, cli_start, cli_end, mhw_start, mhw_end)
    % Constants
    vWindowHalfWidth = 5;
    vsmoothPercentileWidth = 31;
    vminDuration = 5;
    vmaxGap = 2;
    vThreshold = 0.9;

    % Filter valid data
    validIdx = ~isnan(data);
    data = data(validIdx);
    time = time(validIdx);
    depths = depths(validIdx);

    % Unique depths
    uniqueDepths = unique(depths);

    % Initialize outputs
    events_table = [];
    mclim = [];
    positive_threshold = [];
    negative_threshold = [];

    % Process each depth separately
    for d = 1:length(uniqueDepths)
        depthIdx = depths == uniqueDepths(d);
        climData = data(depthIdx);
        climTime = time(depthIdx);

        % Calculate daily mean
        [dailyMeanData, daysOfTimeSeries] = calculate_daily_mean(climData, climTime);

        % Calculate climatology and thresholds for this depth
        [mclim_d, positive_threshold_d, negative_threshold_d] = calculate_climatology_threshold(dailyMeanData, daysOfTimeSeries, cli_start, cli_end, vWindowHalfWidth, vThreshold, vsmoothPercentileWidth);

        % Detect events for this depth using daily mean data
        eventIndex = (daysOfTimeSeries >= mhw_start & daysOfTimeSeries <= mhw_end);
        eventData = dailyMeanData(eventIndex);
        eventTime = daysOfTimeSeries(eventIndex);

        % Detect positive anomalies
        positive_events = detect_events(eventData, eventTime, mclim_d, positive_threshold_d, vminDuration, vmaxGap, 'positive');
        % Detect negative anomalies
        negative_events = detect_events(eventData, eventTime, mclim_d, negative_threshold_d, vminDuration, vmaxGap, 'negative');

        % Combine positive and negative events
        events = [positive_events; negative_events];

        % Collect event details for this depth
        for e = 1:size(events, 1)
            event_details = events(e, :);
            if length(event_details) == 8  % Ensure event_details has the expected length
                % Append the new event details to the events_table
                events_table = [events_table; event_details];
                
                % Update the depth for the newly added row
                events_table(end, 8) = uniqueDepths(d);
            else
                fprintf('Event details dimensions are inconsistent for depth %f: %d\n', uniqueDepths(d), length(event_details));
            end
        end

        % Store climatology and thresholds
        mclim(d, :) = mclim_d;
        positive_threshold(d, :) = positive_threshold_d;
        negative_threshold(d, :) = negative_threshold_d;
    end

    % Convert to table format with appropriate column names
    if ~isempty(events_table)
        events_table = array2table(events_table, 'VariableNames', {'onset', 'end', 'duration', 'max_intensity', 'mean_intensity', 'intensity_variance', 'cumulative_intensity', 'depth'});
        % Convert onset and end to datetime format
        events_table.onset = datetime(events_table.onset, 'ConvertFrom', 'datenum');
        events_table.end = datetime(events_table.end, 'ConvertFrom', 'datenum');
    else
        events_table = array2table(zeros(0, 8), 'VariableNames', {'onset', 'end', 'duration', 'max_intensity', 'mean_intensity', 'intensity_variance', 'cumulative_intensity', 'depth'});
    end
end

% This function calculates the daily mean of the input data, ignoring the times of day.
% Inputs:
%   - data: time series data (temperature or salinity)
%   - time: corresponding time series
% Outputs:
%   - dailyMeanData: daily mean values
%   - daysOfTimeSeries: unique days in the time series
function [dailyMeanData, daysOfTimeSeries] = calculate_daily_mean(data, time)
    % Calculate daily means, ignoring the times of day
    [uniqueDays, ~, idx] = unique(time);
    dailyMeanData = accumarray(idx, data, [], @mean);
    daysOfTimeSeries = uniqueDays;
end

% This function calculates the climatology and thresholds for the input data.
% Inputs:
%   - data: daily mean values
%   - time: unique days in the time series
%   - cli_start: start date for climatology calculation
%   - cli_end: end date for climatology calculation
%   - vWindowHalfWidth: half-width of the moving window for smoothing
%   - vThreshold: threshold percentile for detecting extremes
%   - vsmoothPercentileWidth: width of the moving window for smoothing
% Outputs:
%   - clim: climatology values
%   - positive_thr: threshold for positive anomalies
%   - negative_thr: threshold for negative anomalies
function [clim, positive_thr, negative_thr] = calculate_climatology_threshold(data, time, cli_start, cli_end, vWindowHalfWidth, vThreshold, vsmoothPercentileWidth)
    % Filter data within the climatology period
    cliIdx = (time >= cli_start) & (time <= cli_end);
    climData = data(cliIdx);
    climTime = time(cliIdx);

    % Initialize climatology and threshold arrays
    clim = NaN(1, 366);
    positive_thr = NaN(1, 366);
    negative_thr = NaN(1, 366);

    % Get unique days of the year
    uniqueDays = day(climTime, 'dayofyear');

    % Iterate through each day of the year
    for i = 1:366
        % Find indices of the current day in a leap year
        current_day_idx = (uniqueDays == i);

        % Collect data within the window for the current day
        window_data = [];
        for w = -vWindowHalfWidth:vWindowHalfWidth
            window_day = mod(i + w - 1, 366) + 1;  % Wrap around to handle days outside the range
            window_data = [window_data; climData(uniqueDays == window_day)];
        end

        % Calculate the climatology and threshold for the current day
        if ~isempty(window_data)
            positive_thr(i) = quantile(window_data, vThreshold);
            negative_thr(i) = quantile(window_data, 1 - vThreshold);
            clim(i) = mean(window_data, 'omitnan');
        end
    end

    % Handle Feb 29 by averaging Feb 28 and Mar 1
    if length(clim) >= 60
        positive_thr(60) = mean(positive_thr([59 61]), 'omitnan');
        negative_thr(60) = mean(negative_thr([59 61]), 'omitnan');
        clim(60) = mean(clim([59 61]), 'omitnan');
    end

    % Smooth the climatology and threshold
    clim = smoothdata(clim, 'movmean', vsmoothPercentileWidth);
    positive_thr = smoothdata(positive_thr, 'movmean', vsmoothPercentileWidth);
    negative_thr = smoothdata(negative_thr, 'movmean', vsmoothPercentileWidth);
end

% This function detects extreme events (both positive and negative anomalies) in the input data.
% Inputs:
%   - data: time series data (temperature or salinity)
%   - time: corresponding time series
%   - clim: climatology values
%   - thr: threshold values (positive or negative)
%   - vminDuration: minimum duration of an event
%   - vmaxGap: maximum gap between consecutive event days
%   - type: type of threshold ('positive' or 'negative')
% Outputs:
%   - events: detected events with metrics
function events = detect_events(data, time, clim, thr, vminDuration, vmaxGap, type)
    % Detect MHW or salinity extremes events
    timeSeriesDays = day(datetime(datevec(time)), 'dayofyear');
    unique_days = unique(timeSeriesDays);
    events = [];

    % Initialize variables to track potential events
    potential_start_day_indicator = 0;
    potential_events = [];
    extremeEventIndicator = zeros(size(data));

    % Loop over each unique day
    for day_of_year = 1:length(unique_days)
        if strcmp(type, 'positive')
            day_condition = data(timeSeriesDays == day_of_year) >= thr(day_of_year);
        elseif strcmp(type, 'negative')
            day_condition = data(timeSeriesDays == day_of_year) <= thr(day_of_year);
        end
        extremeEventIndicator(timeSeriesDays == day_of_year) = day_condition;
    end

    % Event detection logic
    for n = 1:length(extremeEventIndicator)
        if potential_start_day_indicator == 0 && extremeEventIndicator(n) == 1
            start_here = datenum(time(n));
            potential_start_day_indicator = 1;
        elseif potential_start_day_indicator == 1 && extremeEventIndicator(n) == 0
            end_here = datenum(time(n)) - 1;
            potential_start_day_indicator = 0;
            potential_events = [potential_events; [start_here end_here]];
        elseif n == length(extremeEventIndicator) && extremeEventIndicator(n) == 1
            end_here = datenum(time(n));
            potential_events = [potential_events; [start_here end_here]];
        end
    end

    % Filter events by duration and gap
    if ~isempty(potential_events)
        potential_events = potential_events((potential_events(:, 2) - potential_events(:, 1) + 1) >= vminDuration, :);

        if ~isempty(potential_events)
            gaps = (potential_events(2:end, 1) - potential_events(1:end-1, 2) - 1);

            % combine events with gaps of <= 2 days
            while min(gaps) <= vmaxGap
                potential_events(find(gaps <= vmaxGap), 2) = potential_events(find(gaps <= vmaxGap) + 1, 2);
                loc_should_del = (find(gaps <= vmaxGap) + 1);
                loc_should_del = loc_should_del(~ismember(loc_should_del, find(gaps <= vmaxGap)));
                potential_events(loc_should_del, :) = [];
                gaps = (potential_events(2:end, 1) - potential_events(1:end-1, 2) - 1);
            end

            % Add a unique event identifier
            event_id = 1;
            
            % Collect event metrics
            for event = 1:size(potential_events, 1)
                current_event_start_end = potential_events(event, :);
            
                % Convert from datenum to day of the year so we can retrieve from clim
                event_start_day = day(datetime(current_event_start_end(1), 'ConvertFrom', 'datenum'), 'dayofyear');
                event_end_day = day(datetime(current_event_start_end(2), 'ConvertFrom', 'datenum'), 'dayofyear');
            
                % Convert current_event_start_end(1) and current_event_start_end(2) to datetime format
                % This is to get the raw data from the original time
                % series, as it needs to be compared to the time array (full data set)
                event_start_datetime = datetime(current_event_start_end(1), 'ConvertFrom', 'datenum');
                event_end_datetime = datetime(current_event_start_end(2), 'ConvertFrom', 'datenum');
            
                % Check if year is a leap year
                is_leap_year = @(year) mod(year, 4) == 0 & (mod(year, 100) ~= 0 | mod(year, 400) == 0);
                leap_year = is_leap_year(year(event_start_datetime));
            
                % Get climatology in relevant time range
                mcl = clim(event_start_day:event_end_day);
            
                % Adjust mcl indices if it's a non-leap year and event spans beyond February 28
                if ~leap_year && event_start_day <= 59 && event_end_day > 59
                    event_start_day_indices = event_start_day:event_end_day;
                    mcl_adjusted_indices = event_start_day_indices;
                    mcl_adjusted_indices(event_start_day_indices > 59) = mcl_adjusted_indices(event_start_day_indices > 59) + 1;
                    mcl = clim(mcl_adjusted_indices);
                end
            
                % Find the indices in 'time' corresponding to the start and end dates
                event_start_index = find(time >= event_start_datetime, 1, 'first');
                event_end_index = find(time <= event_end_datetime, 1, 'last');
            
                % Get anomalies across event
                disp(['Processing Event ID: ', num2str(event_id)]);
            
                % Check if lengths match before calculating anomalies
                if (event_end_index - event_start_index + 1) == length(mcl)
                    event_anomalies = data(event_start_index:event_end_index) - mcl';
            
                    % Ensure lengths match
                    if length(event_anomalies) == length(mcl')
                        max_anomaly = max(event_anomalies);
                        mean_anomaly = mean(event_anomalies);
                        std_anomalies = std(event_anomalies);
                        cum_anomaly = sum(event_anomalies);
            
                        % Preallocate the event_details array to ensure consistent dimensions
                        event_details = NaN(1, 8);
                        event_details(1:8) = [datenum(event_start_datetime), datenum(event_end_datetime), length(event_anomalies), max_anomaly, mean_anomaly, std_anomalies, cum_anomaly, NaN];

                        events = [events; event_details];
                    else
                        disp('Length mismatch between event anomalies and mcl. Skipping this event. (likely data gap)');
                    end
                else
                    disp('Mismatch in mcl and event duration lengths. Skipping this event. (likely data gap)');
                end
            
                % Increment the event identifier
                event_id = event_id + 1;
            end
        end
    end
end