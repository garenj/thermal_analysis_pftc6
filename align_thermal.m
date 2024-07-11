% [cm mm aa stats]= align_thermal('/Users/benjaminblonder/Documents/rmbl/rmbl 2016/rmbl thermal ecology/thermal data/pfeiler jun 27/thermal/combined/', 2, 1, 0.5, 200);

% JCG Test
% [cm mm aa stats]= align_thermal('/Users/joe/Desktop/flir_thermal_control-jcg-update/test_radiometric_2/', 2, 1, 0.01, 1);
% [cm mm aa stats]= align_thermal('/Users/joe/Desktop/thermal_test_run/', 12, 1, 0.01, 1);

% interval_keyframes
% interval_frames
% threshold
% keyframe_id
function [correctedMean movMean array_aligned stats] = align_thermal(folder_in_thermal_timeseries, interval_keyframes, interval_frames, threshold, keyframe_id)
    addpath('npy-matlab-master');
     
    % set alignment parameters
    if (nargin < 2)
        interval_keyframes = 10;
        interval_frames = 1;
    end

    
    % load in thermal images
    %files_thermal_timeseries_npy = dir([folder_in_thermal_timeseries '*.npy']);
    files_thermal_timeseries_npy = dir(strcat(folder_in_thermal_timeseries, "*.npy"))
    
    
    keepfiles = ones([length(files_thermal_timeseries_npy) 1]);
    for i=1:length(files_thermal_timeseries_npy)
        fn = fullfile(folder_in_thermal_timeseries, files_thermal_timeseries_npy(i).name);
        try
            im_current = double(readNPY(fn));

            zeros = nnz(im_current==0);
            if (zeros > 50)
                fprintf('drop %d %d %s\n', i, zeros, fn);
                keepfiles(i) = 0;
            else
                fprintf(' *%d* %s\n',i, fn);
            end
            fprintf('\n');
        catch
            fprintf('drop %d %d %s\n', i, zeros, fn);
            keepfiles(i) = 0;           
        end
        
    end
    files_thermal_timeseries_npy = files_thermal_timeseries_npy(keepfiles>0);
    numfiles_thermal_timeseries = length(files_thermal_timeseries_npy);
    
    % get first frame
    movMean = double(readNPY(fullfile(folder_in_thermal_timeseries, files_thermal_timeseries_npy(keyframe_id).name)));
    movMean = imgaussfilt(rescale_image_quantile(movMean, 0.01, 0.99),2);
    imgB = movMean;
    imgBp = imgB;
    correctedMean = imgBp;
    Hcumulative = eye(3);
    count_average = 0;
    
    transform_current = Hcumulative;
    
    indexvals = 1:interval_frames:numfiles_thermal_timeseries;        
    % allocate memory
    array_aligned = repmat(uint16(zeros(1)),[480 640 length(indexvals)]);

    stats = cell([length(indexvals) 1]);
    
    figure;
    for i=1:length(indexvals)
        fn = fullfile(folder_in_thermal_timeseries, files_thermal_timeseries_npy(indexvals(i)).name);

        fn_stats = strrep(fn,'infrared-data.npy', 'stats.csv');

        try
            stats{i} = readtable(fn_stats);
            stats{i}.gps_time = {stats{i}.gps_time};        
        end
        
        try
            if (mod(indexvals(i),interval_keyframes)==1)
                fprintf('*');
                % Move old frames
                imgA = imgB; % z^-1
                imgAp = imgBp; % z^-1
                % Read in new frame
                imgB = double(readNPY(fn));
                imgB_untransformed = imgB;
                imgB = imgaussfilt(rescale_image_quantile(imgB, 0.01, 0.99),2);
                movMean = movMean + imgB;
                %disp(imgA);
                %disp(imgB);

                % do stabilization transform and warp
                %H = cvexEstStabilizationTform(imgA,imgB,threshold);
                %HsRt = cvexTformToSRT(H);
                %Hcumulative = HsRt * Hcumulative;
                %imgBp = imwarp(imgB,affine2d(Hcumulative),'OutputView',imref2d(size(imgB)));
                % add new value to mean
                imgBp = imgB; % THIS GETS RID OF THE STABILIZATION
                correctedMean = correctedMean + imgBp;
                
                transform_current = Hcumulative;
                
                count_average = count_average + 1;
                fprintf('%d', count_average);
            end
        end
        
        im_current = double(readNPY(fn));

        % transform the raw image too
        %im_current_warped = imwarp(im_current,affine2d(transform_current),'OutputView',imref2d(size(im_current)));         
        im_current_warped = im_current; % THIS GETS RID OF THE STABILIZATION
        imshow(rescale_image_quantile(im_current_warped, 0.05, 0.95));
        % reconvert back to uint16
        array_aligned(:,:,i) = uint16(im_current_warped);
        
        fprintf('%d ', indexvals(i))

    end
    
    correctedMean = correctedMean/(count_average);
    movMean = movMean/(count_average);

    fprintf('\n')
    fprintf('%d', count_average)
end



