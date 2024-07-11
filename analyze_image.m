%x = analyze_image('/Users/benjaminblonder/Documents/rmbl/rmbl 2016/rmbl thermal ecology/thermal data/baldy 28 jul/baldy_partone_2016_07_28.mat', '/Users/benjaminblonder/Documents/rmbl/rmbl 2016/rmbl thermal ecology/thermal data/baldy 28 jul/visible/ERIUMB1-1 (0142).jpg','~/Downloads/out/');

%x = analyze_image("jcg_test_output.mat", '/Users/joe/Desktop/flir_thermal_control-jcg-update/test_radiometric_2/out_200928-165122_0-visible.png',"test_out_folder")

%x = analyze_image("/Users/joe/Desktop/thermal_test_run/jcg_test_output.mat", '/Users/joe/Desktop/thermal_test_run/out_201008-150019_0-visible.png',"test_out_folder/")
%x = analyze_image('/Users/joe/Desktop/licor_thermal_data/results_200731-103709/licor_6800_1.mat', "/Users/joe/Desktop/licor_thermal_data/results_200731-103709/out_200731-103724_0-visible.png", "/Users/jose/Desktop/licor_thermal_data/6800_side")

function [image_thermal_representative] = analyze_image(file_thermal_mat, file_visible, folder_out)
    scalefactor = 2;

    if (nargin < 3)
        folder_out = 'output';
    end
    
    if ~exist(folder_out, 'dir')
        mkdir(folder_out);
    end
    
    data_thermal = matfile(file_thermal_mat); % load slice-by-slice
    image_visible_lores_registered = data_thermal.image_visible_lores_registered;
    stats = data_thermal.finalstats;
    time_elapsed = stats.time_elapsed;
    num_files = size(data_thermal.finalstats, 1);
    
    % preprocess thermal image
    image_thermal_representative = double(data_thermal.Tkelvin_aligned_calibrated(:,:,floor(num_files/3)))/100;
    image_thermal_representative = imresize(image_thermal_representative, scalefactor);
    image_thermal_representative = rescale_image_quantile(image_thermal_representative, 0.05, 0.95);
    %image_thermal_mean = adapthisteq(image_thermal_mean, 'NumTiles',[8 8]);
    image_thermal_representative = imsharpen(image_thermal_representative,'Radius',2,'Amount',1.5);
    image_thermal_representative = ind2rgb(floor(255*image_thermal_representative),hot(255));
    fprintf('\n')

    points_thermal = [];
    points_visible = [];
    
    image_visible_hires = imread(file_visible);

    done = false;
    doalign = true;
    while (~done)
        if (doalign)
            [image_fused, image_visible_registered, points_thermal, points_visible] = image_align(image_thermal_representative, image_visible_hires, points_thermal, points_visible);
        end
        
        doanotherroi = true;
        
        while (doanotherroi==true)
            [bw finishroi] = chooseroi(image_visible_lores_registered, image_visible_registered, image_thermal_representative);

            if (strcmp(finishroi, 'No - realign'))
                doalign = true;
            else % if we did get a good ROI
                doalign = false;

                bw = imresize(bw, [480 640]);
                pixels_keep = bw>0;

                % calculate stats in this region
                thermal_stats = zeros([num_files 5]);
                for j=1:num_files
                    image_visible_hires = double(data_thermal.Tkelvin_aligned_calibrated(:,:,j))/100;
                    pixels_this = image_visible_hires(pixels_keep);
                    pixels_this = pixels_this(pixels_this>0);
                    thermal_stats(j,1) = mean(pixels_this);
                    thermal_stats(j,2) = std(pixels_this);
                    thermal_stats(j,3) = quantile(pixels_this,0.05);
                    thermal_stats(j,4) = quantile(pixels_this,0.5);
                    thermal_stats(j,5) = quantile(pixels_this,0.95);
                    fprintf('%.4f\n',j/num_files);
                end
                fprintf('\n')
                % show plot
                f2 = figure;
                subplot(1,2,1); 
                plot(time_elapsed, thermal_stats(:,3),'-g'); hold on;
                plot(time_elapsed, thermal_stats(:,4),'-r');
                plot(time_elapsed, thermal_stats(:,5),'-b');
                plot(time_elapsed, stats.temp_atm,'-k');
                xlabel('Time (seconds');
                ylabel('Temperature');
                subplot(1,2,2); 
                plot(stats.temp_atm, thermal_stats(:,3),'.r'); hold on;
                xlabel('Air temp'); ylabel('Leaf temp');
                xlim([263 343]);
                ylim([263 343]);
                line([0 500], [0 500]);

                ans_roi = MFquestdlg([0.5 0.5], 'Keep this ROI?','Prompt','yes','no','finished','yes');
                if (strcmp(ans_roi,'yes')==1)
                    [~, fp, ~] = fileparts(file_visible);
                    disp(fp);
                    filename_output = inputdlg('Enter sample name','Input',1,fp);
                    filename_output = filename_output{1};

                    if (~isempty(filename_output))
                        table_out = table;
                        table_out.file = repmat(file_thermal_mat, [num_files 1]);
                        table_out.thermal_mean = thermal_stats(:,1);
                        table_out.thermal_sd = thermal_stats(:,2);
                        table_out.thermal_q05 = thermal_stats(:,3);
                        table_out.thermal_q50 = thermal_stats(:,4);
                        table_out.thermal_q95 = thermal_stats(:,5);

                        table_out = horzcat(table_out, stats);

                        filename_output_final = sprintf('%s%s.csv',folder_out,filename_output);
                        writetable(table_out, filename_output_final)
                        fprintf('wrote file %s\n', filename_output_final);

                        imwrite(bw, sprintf('%s/%s-mask.png',folder_out,filename_output));
                        imwrite(image_visible_registered, sprintf('%s/%s-visible-registered.jpg',folder_out,filename_output));
                        done = true;
                    end
                elseif (strcmp(ans_roi,'finished')==1)
                    done = true;
                end

                close(f2);
            end
            
            ans_roi = questdlg('Do another ROI in this aligned image?','Query','Yes','No','Yes');
            doanotherroi = strcmp(ans_roi, 'Yes');
        end
    end

end