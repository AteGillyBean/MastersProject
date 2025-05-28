clc; % command window
clear all; % workspace
close all; % close all pictures

%% Load Files
segmentationFile = 'KJ-SOL.nrrd'; % what muscle segmentation we're looking at
[segmentation, metaSeg] = nrrdread(segmentationFile);

fatDICOMFolder = '/Users/gillianaloyola/Documents/Research/Segs/KJ Final/t1_vibe_dixon_tra_2echo_F_22';
waterDICOMFolder = '/Users/gillianaloyola/Documents/Research/Segs/KJ Final/t1_vibe_dixon_tra_2echo_W_23';

% Read DICOM volumes
fatImage = squeeze(double(dicomreadVolume(fatDICOMFolder)));
waterImage = squeeze(double(dicomreadVolume(waterDICOMFolder)));
segmentation = double(segmentation);

if ~isequal(size(segmentation), size(fatImage), size(waterImage))
    error('Mismatch in image dimensions. Resample segmentation to match DICOM images.');
end

%% Erode Segmentation Mask
se = strel('square', 3);
erodedSegmentation = imerode(segmentation, se);

% %% TO GET BACKGROUND First (run only once)
% % Extract background intensities from selected points
% for i = 1:length(x)
%     bgIntensity(i) = fatImage(round(y(i)), round(x(i)), sliceIdx);
% end
% 
% disp('Background Intensities:');
% disp(bgIntensity);
% 
% % Compute statistics for background
% backgroundVoxels = fatImage(segmentation == 0); % Pixels where segmentation is 0
% meanBackground = mean(backgroundVoxels);
% stdBackground = std(backgroundVoxels);
% 
% fprintf('Mean Background Intensity: %.2f\n', meanBackground);
% fprintf('Standard Deviation of Background: %.2f\n', stdBackground);
%%
clc;
% equation: fat img / (fat + water imgs)

% Compute Fat Fraction
threshold = 53; %background - take background mean from earlier
denominator = fatImage + waterImage;
validMask = (denominator > threshold & denominator > 0);
fatFraction = zeros(size(fatImage));
fatFraction(validMask) = fatImage(validMask) ./ denominator(validMask);

% Remove values greater than 0.3
fatFraction(fatFraction > 0.3) = 0; % Set values above 0.3 to zero

% Compute Muscle Fat Fraction
muscleFatFraction = fatFraction .* erodedSegmentation;

muscleVoxels = muscleFatFraction(erodedSegmentation == 1);
muscleVoxels = muscleVoxels(muscleVoxels > 0);

fprintf('Mean Fat Fraction in Eroded Muscle: %.4f (%.2f%%)\n', mean(muscleVoxels), mean(muscleVoxels) * 100);
fprintf('Standard Deviation of Fat Fraction in Eroded Muscle: %.4f\n', std(muscleVoxels));


%% Compute Fat Fraction (Original Segmentation)
originalMuscleVoxels = fatFraction(segmentation == 1);
originalMuscleVoxels = originalMuscleVoxels(originalMuscleVoxels > 0);

fprintf('Mean Fat Fraction in Original Muscle: %.4f\n', mean(originalMuscleVoxels));
fprintf('Standard Deviation of Fat Fraction in Original Muscle: %.4f\n', std(originalMuscleVoxels));

% %% Compare to Data from paper
% %clc;
% 
% % Define the normal fat fraction range (soleus)
% normal_min = 2.2 / 100; normal_max = 4.4 / 100;
% 
% %normal_min = 0.7 / 100; normal_max = 1.9 / 100; %TA Values
% 
% %normal_min = 2.4 / 100; normal_max = 5.8 / 100; %MG Values
% 
% %normal_min = 1.0 / 100; normal_max = 5.7 / 100; %LG Values
% 
% % Compute the mean fat fraction
% mean_fat_fraction = mean(muscleVoxels);
% 
% % Calculate the differences
% diff_min = mean_fat_fraction - normal_min;
% diff_max = mean_fat_fraction - normal_max;
% 
% % Display results
% fprintf('Mean Fat Fraction in Eroded Muscle: %.4f (%.2f%%)\n', mean_fat_fraction, mean_fat_fraction * 100);
% fprintf('Difference from Normal Min (2.2%%): %.4f (%.2f%%)\n', diff_min, diff_min * 100);
% fprintf('Difference from Normal Max (4.4%%): %.4f (%.2f%%)\n', diff_max, diff_max * 100);

%% For JH - LG %segment check

nonzeroSlices = find(squeeze(sum(sum(erodedSegmentation, 1), 2)) > 0);
if ~isempty(nonzeroSlices)
    sliceIdx = nonzeroSlices(round(end/2)); % Choose a slice with valid segmentation
else
    error('No valid slices found in eroded segmentation.');
end

figure;
imshow(erodedSegmentation(:,:,sliceIdx), []);
title('Eroded Segmentation Check');
%% Visualization
%sliceIdx = round(size(fatImage, 3) / 2); %everything but JH LG
sliceIdx = nonzeroSlices(round(end/2)); % Choose a slice with valid segmentation

% fat & water images with both segmentations
figure('Name', 'Fat and Water Images', 'NumberTitle', 'off', 'WindowState', 'maximized');
tiledlayout(2,2);
nexttile;
imshow(mat2gray(fatImage(:,:,sliceIdx))); colormap jet; colorbar; title('Fat Image');
nexttile;
imshow(mat2gray(waterImage(:,:,sliceIdx))); colormap jet; colorbar; title('Water Image');
nexttile;
imshow(segmentation(:,:,sliceIdx)); colormap gray; title('Original Segmentation');
nexttile;
imshow(erodedSegmentation(:,:,sliceIdx)); colormap gray; title('Eroded Segmentation');

% fat and muscle fractions
figure('Name', 'Fat Fraction Maps', 'NumberTitle', 'off', 'WindowState', 'maximized');
tiledlayout(1,2);
nexttile;
imshow(mat2gray(fatFraction(:,:,sliceIdx))); colormap jet; colorbar; title('Fat Fraction');
nexttile;
imshow(mat2gray(muscleFatFraction(:,:,sliceIdx))); colormap jet; colorbar; title('Muscle Fat Fraction');

% Overlay of Eroded Segmentation
figure('Name', 'Segmentation Overlay on Water Image', 'NumberTitle', 'off', 'WindowState', 'maximized');
imshow(waterImage(:,:,sliceIdx), []);
hold on;
contour(erodedSegmentation(:,:,sliceIdx), 'r', 'LineWidth', 2);
title('Eroded Segmentation Overlay on Water Image');

%% Comparing Muscle Fat Fraction B/W Segmentations
% Compute Fat Fraction for Original and Eroded Segmentations
originalMuscleFatFraction = fatFraction .* segmentation;
erodedMuscleFatFraction = fatFraction .* erodedSegmentation;

% Visualize the comparison
figure('Name', 'Fat Fraction in Muscle - Original vs Eroded Segmentation', 'NumberTitle', 'off', 'WindowState', 'maximized');
tiledlayout(1,2);

% Original Muscle Fat Fraction
nexttile;
imshow(mat2gray(originalMuscleFatFraction(:,:,sliceIdx))); 
colormap jet; colorbar;
title('Original Segmentation - Muscle Fat Fraction');

% Eroded Muscle Fat Fraction
nexttile;
imshow(mat2gray(erodedMuscleFatFraction(:,:,sliceIdx))); 
colormap jet; colorbar;
title('Eroded Segmentation - Muscle Fat Fraction');

%% Comparing Muscle Fat Fraction B/W Segmentations
% Compute Fat Fraction for Original and Eroded Segmentations
originalMuscleFatFraction = fatFraction .* segmentation;
erodedMuscleFatFraction = fatFraction .* erodedSegmentation;

% Visualize the comparison
figure('Name', 'Fat Fraction in Muscle - Original vs Eroded Segmentation', 'NumberTitle', 'off', 'WindowState', 'maximized');
tiledlayout(1,2);

% Original Muscle Fat Fraction
nexttile;
imshow(mat2gray(originalMuscleFatFraction(:,:,sliceIdx))); 
colormap jet; colorbar;
title('Original Segmentation - Muscle Fat Fraction');

% Eroded Muscle Fat Fraction
nexttile;
imshow(mat2gray(erodedMuscleFatFraction(:,:,sliceIdx))); 
colormap jet; colorbar;
title('Eroded Segmentation - Muscle Fat Fraction');

%% Save Results in Structured Folder
% Base folder
baseFolder = 'Results';
if ~exist(baseFolder, 'dir')
    mkdir(baseFolder);
end

customName = 'KJ-SOL'; % <<< Change this to whatever you want
subFolder = fullfile(baseFolder, customName);
if ~exist(subFolder, 'dir')
    mkdir(subFolder);
end

% --- Save Figures ---25 ;9
figHandles = findall(groot, 'Type', 'figure');
for i = 1:length(figHandles)
    figName = get(figHandles(i), 'Name');
    if isempty(figName)
        figName = ['Figure_' num2str(figHandles(i).Number)];
    end
    saveas(figHandles(i), fullfile(subFolder, [figName '.png']));
end

% --- Compute and Save Statistics ---
meanEroded = mean(muscleVoxels);
stdEroded = std(muscleVoxels);
meanOriginal = mean(originalMuscleVoxels);
stdOriginal = std(originalMuscleVoxels);

% Save as .txt
statsText = sprintf([ ...
    'Mean Fat Fraction in Eroded Muscle: %.4f (%.2f%%)\n', ...
    'Standard Deviation in Eroded Muscle: %.4f\n', ...
    'Mean Fat Fraction in Original Muscle: %.4f\n', ...
    'Standard Deviation in Original Muscle: %.4f\n'], ...
    meanEroded, meanEroded * 100, ...
    stdEroded, ...
    meanOriginal, ...
    stdOriginal);

txtFile = fullfile(subFolder, 'FatFractionStats.txt');
fid = fopen(txtFile, 'w');
fprintf(fid, '%s', statsText);
fclose(fid);

% Save as .csv
statsCell = {
    'Metric', 'Value';
    'Mean Fat Fraction in Eroded Muscle', meanEroded;
    'Std Dev in Eroded Muscl', stdEroded;
    'Mean Fat Fraction in Original Muscle', meanOriginal;
    'Std Dev in Original Muscle', stdOriginal;
};

csvFile = fullfile(subFolder, 'FatFractionStats.csv');
writecell(statsCell, csvFile);

close all;