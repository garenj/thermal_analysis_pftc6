
% JCG Test
% [Tkelvin_aligned_calibrated, finalstats] = calibrate_thermal(aa, stats, 263, 343, 0, "none", 0, '/Users/joe/Desktop/flir_thermal_control-jcg-update/results_test_2/out_200916-160200_0-visible.png', 'jcg_test_output.mat');

% [Tkelvin_aligned_calibrated, finalstats] = calibrate_thermal(aa, stats, 263, 343, 0, "none", 0, '/Users/joe/Desktop/flir_thermal_control-jcg-update/results_test_4/out_200924-133452_0-visible.png', 'jcg_test_output.mat');

% [Tkelvin_aligned_calibrated, finalstats] = calibrate_thermal(aa, stats, 263, 343, 1, "none", 0, '/Users/joe/Desktop/thermal_test_run/out_201008-150019_0-visible.png', '/Users/joe/Desktop/thermal_test_run/jcg_test_output.mat');

% assumes that stats (from camera) are in the format of 
% [hours, minutes, seconds, temp_black, temp_refl, temp_sky, ...]

function [Tkelvin_aligned_calibrated, finalstats] = calibrate_thermal(image_array, stats, bound_temp_lo, bound_temp_hi, dogroundcalibration, xlsinputname, time_offset, file_visible_lores, outputname)
%    if dogroundcalibration==1
%        xls_raw = xlsread(xlsinputname);
%        xls_time = (xls_raw(:,1)*60 + xls_raw(:,2)) * 60; % convert to seconds
%        xls_time_elapsed = xls_time - xls_time(1);
%        xls_temp_black = xls_raw(:,4) + 273.15;
%        xls_temp_refl = xls_raw(:,5) + 273.15;
%
%    end
    
    num_files = size(image_array,3);
    
    % calculate time deltas
    time_raw = nan([length(stats), 6]);
    time_elapsed = zeros([length(stats), 1]);
    for i=1:length(stats)        
        %try
            %thisdate = stats{i}.gps_utc;
            %time_raw(i,:) = datevec(thisdate, 'yyyy-mm-ddTHH:MM:SS.000Z');
        %end
        try
            time_raw(i,:) = datevec(stats{i}.Date,'yymmdd-HHMMSS');
        end
    end
    % convert to elapsed time
    for i=1:length(stats)
        try
            time_elapsed(i) = etime(time_raw(i,:),time_raw(1,:));
        end
    end
    
    time_elapsed
    
    %return;
    
    image_mean = zeros([size(image_array,1) size(image_array,2)]);
    for i=1:num_files
        image_mean = image_mean + double(image_array(:,:,i));
        fprintf('%d\n',i);
    end
    image_mean = image_mean / num_files;
    
    % I think this is only used in the image alignment thing - can we just
    % use the first image? - no image_mean is used for selecting black
    % reference roi
    %image_mean = double(image_array(:,:,1));
    
    % fill in other stats
    temp_atm = 293*ones([length(stats), 1]);
    temp_reflected = 293*ones([length(stats), 1]);
    temp_external_optics = 293.*ones([length(stats), 1]);
    relative_humidity = 0.5*ones([length(stats), 1]);
    emissivity = 0.97*ones([length(stats), 1]);
    distance_focal = 2*ones([length(stats), 1]);
    light_ir = 0*ones([length(stats), 1]);
    light_vis = 0*ones([length(stats), 1]);
    light_uv = 0*ones([length(stats), 1]);
    pressure = NaN([length(stats), 1]);
    temp_sens = 293*ones([length(stats), 1]);
    latitude = NaN([length(stats), 1]);
    longitude = NaN([length(stats), 1]);
    
    %%%%%%%%%%%
    %temp_ambient = 293*ones([length(stats), 1]);
    temp_soil_c = 293*ones([length(stats), 1]);
    ppfd_mV = 0*ones([length(stats), 1]);
    ppfd_umol_m2_s = 0*ones([length(stats), 1]);
    %temp_soil_2 = 293*ones([length(stats), 1]);
    %temp_soil_3 = 293*ones([length(stats), 1]);
    %temp_soil_mean = 293*ones([length(stats), 1]);

    %%%%%%%%%%%
    X_in = 0*ones([length(stats), 1]);
    alpha_1_in = 0*ones([length(stats), 1]);
    alpha_2_in = 0*ones([length(stats), 1]);
    beta_1_in = 0*ones([length(stats), 1]);
    beta_2_in = 0*ones([length(stats), 1]);
    ext_optics_trans_in = 0*ones([length(stats), 1]);
    ext_optics_temp_in = 0*ones([length(stats), 1]);
    J_0_in = 0*ones([length(stats), 1]);
    J_1_in = 0*ones([length(stats), 1]);
    R_in = 0*ones([length(stats), 1]);
    F_in = 0*ones([length(stats), 1]);
    B_in = 0*ones([length(stats), 1]);
    
    for i=1:length(stats)
        try
            %temp_atm(i) = stats{i}.wx_temp_air_c + 273.15; %in C, convert to K
            temp_atm(i) = stats{i}.tc_amb_c + 273.15;
            print(temp_atm(i));
        end   
        
        try
            %temp_reflected(i) = temp_refl(i); % What is happening here? This variable isn't defined. Do they mean the xls calib value?
            temp_reflected(i) = stats(i).ReflectedTemperature
            print(temp_reflected(i));
        end  

        try
            temp_external_optics(i) = stats{i}.TSens;
            print(temp_external_optics(i));
        end 
        
        % Humidity sensor broken
        %try
            %relative_humidity(i) = stats{i}.wx_rel_hum;
            %relative_humidity(i) = stats{i}.RelativeHumidity;
            %print(relative_humidity(i));
        %end
        
        try
            distance_focal(i) = stats{i}.FocusDistance;
            if (distance_focal(i) > 2)
                distance_focal(i) = 2;
            end
        end
        
        try
            light_ir(i) = stats{i}.wx_ir_lux;
        end
        
        try
            light_vis(i) = stats{i}.wx_vis_lux;
        end
        
        try
            light_uv(i) = stats{i}.wx_uv_lux;
        end
        
        try
            pressure(i) = stats{i}.wx_pressure_hpa;
        end
        
        try
            temp_sens(i) = stats{i}.TSens;
        end
        
        try
            latitude(i) = stats{i}.gps_latitude;
        end
        
        try
            longitude(i) = stats{i}.gps_longitude;
        end
        
        try % my test block jcg
            X_in(i) = stats{i}.X;
            alpha_1_in(i) = stats{i}.alpha1;
            alpha_2_in(i) = stats{i}.alpha2;
            beta_1_in(i) = stats{i}.beta1;
            beta_2_in(i) = stats{i}.beta2;
            ext_optics_trans_in(i) = stats{i}.ExtOpticsTransmission;
            ext_optics_temp_in(i) = stats{i}.ExtOpticsTemperature;
            J_0_in(i) = stats{i}.J0;
            J_1_in(i) = stats{i}.J1;
            R_in(i) = stats{i}.R;
            F_in(i) = stats{i}.F;
            B_in(i) = stats{i}.B;
        end
        try
            %temp_ambient(i) = stats{i}.tc_amb_c + 273.15;
            temp_soil_c(i) = stats{i}.tc_soil1_c + 273.15;
            ppfd_mV(i) = stats{i}.ppfd_mV_raw;
            ppfd_umol_m2_s(i) = stats{i}.ppfd_umol_m2_s;
            %temp_soil_2(i) = stats{i}.tc_soil2_c + 273.15;
            %temp_soil_3(i) = stats{i}.tc_soil3_c + 273.15;
            %temp_soil_mean(i) = (stats{i}.tc_soil1_c+stats{i}.tc_soil2_c+stats{i}.tc_soil3_c)/3 + 273.15;
        end
            
    end
    
    % convert radiometric counts to temperatures using trefl values
    %temp_array = repmat(single(0), size(image_array));
    temp_array = double(repmat(uint16(0), size(image_array)));
    for i=1:length(stats)
        fprintf('*');
        ithis = image_array(:,:,i);
        %ithis(ithis==0) = NaN; % remove any extraneous values values
%        temp_this = calibrated_temperature_simple(...
%                double(image_array(:,:,i)), ...
%                temp_atm(i), ...
%                temp_reflected(i), ...
%                temp_external_optics(i), ...
%                relative_humidity(i), ...
%                emissivity(i), ...
%                distance_focal(i) ...
%            );

            temp_this = calibrated_temperature(...
                 double(image_array(:,:,i)), ... %double(counts), ... % input data
                 relative_humidity(i), ... %relative_humidity, ... % relative humidity (fraction)
                 temp_atm(i), ... %temp_atm, ... % atmospheric temp (K)
                 distance_focal(i), ...%distance_focal, ... % object distance (m)
                 X_in(i), ... %1.9, ... % X
                 alpha_1_in(i), ... %0.006569, ... % alpha_1
                 beta_1_in(i), ... %-0.002276, ... % beta_1
                 alpha_2_in(i), ... %0.01262, ... % alpha_2
                 beta_2_in(i), ... %-0.00667, ... % beta_2
                 emissivity(i), ... %emissivity, ... % emissivity
                 ext_optics_trans_in(i), ... %1, ... % external optics transmission
                 temp_reflected(i), ... %temp_reflected, ... % reflected (ambient) temperature
                 ext_optics_temp_in(i), ... %temp_external_optics, ... % external optics temperature
                 J_0_in(i), ... %4214, ... % J0
                 J_1_in(i), ... %69.62449646, ... % J1
                 R_in(i), ... %16671, ... % R
                 F_in(i), ... %1, ... % F
                 B_in(i) ... %1430.099976... % B
            );

%        disp(temp_this(250:260,250:260)));
%        temp_this = double(image_array(:,:,i)) / 100.0;
        % remove further extraneous values
        %temp_this(temp_this < bound_temp_lo) = NaN;
        %temp_this(temp_this > bound_temp_hi) = NaN;
        temp_this(temp_this < bound_temp_lo) = 0;
        temp_this(temp_this > bound_temp_hi) = 0;
        %temp_array(:,:,i) = temp_this;
        temp_array(:,:,i) = temp_this*100;
%        disp(temp_this);
    end
    fprintf('\n');
    
    if (dogroundcalibration==0)
        temp_black = NaN([length(time_elapsed) 1]);
    end
    
    if (dogroundcalibration==1)
        % choose region of interest
        
        f1 = figure('Name','Select ROI for black reference');
        %image_mean
        imshow(rescale_image_quantile(image_mean,0.4,0.85));
        bw = roipoly;
        pixels_keep = bw>0;
        %print(pixels_keep);
        close(f1);

        % calculate stats in this region
        temperature_stats = zeros([num_files 5]);
        for i=1:num_files
            %temp_this = temp_array(:,:,i);
            temp_this = double(temp_array(:,:,i)) / 100;
            pixels_this = temp_this(pixels_keep);
            % keep non-NA pixels
            %pixels_this = pixels_this(~isnan(pixels_this));
            pixels_this = pixels_this(pixels_this>0);

            temperature_stats(i,1) = mean(pixels_this);
            temperature_stats(i,2) = std(pixels_this);
            temperature_stats(i,3) = quantile(pixels_this,0.05);
            temperature_stats(i,4) = quantile(pixels_this,0.5);
            temperature_stats(i,5) = quantile(pixels_this,0.95); 

            fprintf('.');
        end
        fprintf('\n')
        

        
        regressionfitagain = true;
        while (regressionfitagain==true)
%            inputans = inputdlg({'Start time (s)','Stop time (s)','Delta (s)'},'Select index range', 1, {'0','45000',sprintf('%d', time_offset)}); 
%            start_time = str2num(inputans{1});
%            stop_time = str2num(inputans{2});
%            time_offset = str2num(inputans{3});
%            index_start = find(time_elapsed >= (start_time),1,'first');
%            index_stop = find(time_elapsed < (stop_time),1,'last');

            % get values of the xls temperatures at these times
%            if dogroundcalibration==1
%                temp_black = interp1(xls_time_elapsed - time_offset, xls_temp_black, time_elapsed,'spline',NaN);
%            end

            % get those temp_black values
            temp_black = 293*ones([length(stats), 1]);
            for i=1:length(stats)
                temp_black(i) = stats{i}.tc_black_c+273.15;
            end


                        
            ts = temperature_stats(:,4); %(index_start:index_stop,4);
            [b,~,~,~,~] = regress(temp_black, [ones(size(ts)) ts ts.^2 ]);
            
            % show plots
            f2 = figure('Name','Black ref (black) uncalibrated (green) xlstime (blue) recalibrated (red)');
            plot(time_elapsed, temp_black,'-k'); hold on;
            plot(time_elapsed, temperature_stats(:,4),'-g'); hold on;

%            plot(xls_time_elapsed - time_offset, xls_temp_black,'-b'); hold on; 
%            plot(time_elapsed(index_start:index_stop), b(1) + b(2) * ts + b(3) * ts.^2, '-r');
%            plot(time_elapsed - time_offset, temp_black,'-b'); hold on; 
            plot(time_elapsed, b(1) + b(2) * ts + b(3) * ts.^2, '-r');
            
            drawnow; 
            
            regressionfitagain = ~strcmp(questdlg('Done?','Fit','Yes','No','Yes'),'Yes');
            
            close(f2);
        end

        Tkelvin_aligned_calibrated = uint16(temp_array);
        for i=1:num_files
            temp_this = double(temp_array(:,:,i))/100;
            Tkelvin_aligned_calibrated(:,:,i) = uint16(100*( b(1) + b(2) * temp_this + b(3) * temp_this.^2 ) );
            fprintf('|');
        end
        fprintf('\n');
        %Tkelvin_aligned_calibrated = b(1) + (b(2)*100)*temp_array + (b(3)*100^2)*temp_array.^2;
        %Tkelvin_aligned_calibrated = temp_array;

        temp_new = zeros([num_files 1]);
        for i=1:num_files
            temp_this = double(Tkelvin_aligned_calibrated(:,:,i))/100;
            pixels_this = temp_this(pixels_keep);
            pixels_this = pixels_this(pixels_this>0);
            temp_new(i) = quantile(pixels_this,0.5);
            fprintf('/');
        end
        fprintf('\n');


        plot(time_elapsed, temp_new,'-r')


    else
        Tkelvin_aligned_calibrated = temp_array;
        %disp(temp_array);
    end
    
    
%     mediantemp = zeros([size(Tkelvin_aligned_calibrated,3) 1]);
%     mediantime = zeros([size(Tkelvin_aligned_calibrated,3) 1]);
%     for i=1:size(Tkelvin_aligned_calibrated,3)
%         pixels_all = double(Tkelvin_aligned_calibrated(:,:,i))/100;
%         %pixels_all = pixels_all(~isnan(pixels_all));
%         pixels_all = pixels_all(pixels_all>0);
%         
%         mediantemp(i) = median(pixels_all);
%         mediantime(i) = time_elapsed(i);
%         fprintf('-');
%     end
%     fprintf('\n');
    
    %f4 = figure('Name','Image median (magenta)');
    %plot(time_elapsed, mediantemp,'-m');
    
    finalstats = table(time_elapsed, temp_black, temp_atm, temp_soil_c, ppfd_mV, ppfd_umol_m2_s, temp_reflected, temp_external_optics, relative_humidity, emissivity, distance_focal, light_vis, light_ir, light_uv, pressure, temp_sens, latitude, longitude, time_raw); 

    % do visible alignment
    image_visible_lores = imread(file_visible_lores);
    points_thermal_lores = [];
    points_visible_lores = [];
    whichim = floor(size(Tkelvin_aligned_calibrated,3)/2); % another attempt
    image_thermal_representative = double(Tkelvin_aligned_calibrated(:,:,whichim))/100;
    image_thermal_representative = imresize(image_thermal_representative, 2);
    image_thermal_representative = rescale_image_quantile(image_thermal_representative, 0.05, 0.95);
    image_thermal_representative = imsharpen(image_thermal_representative,'Radius',2,'Amount',1.5);
    image_thermal_representative = ind2rgb(floor(255*image_thermal_representative),hot(255));
    
    [image_fused_lores, image_visible_lores_registered, points_thermal_lores, points_visible_lores] = image_align(image_thermal_representative, image_visible_lores, points_thermal_lores, points_visible_lores); 
    
    dosave = questdlg('Save matrix of output','Do save?','yes','no','yes');
    if (strcmp(dosave,'yes'))
        if dogroundcalibration==0
            b = NaN; % this allows for partial loading
        end
        
        save(outputname, 'Tkelvin_aligned_calibrated', 'finalstats','image_visible_lores_registered','b','time_offset', '-v7.3'); % this allows for partial loading
    end
end