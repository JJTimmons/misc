function FLAIR_path = find_flair()
% Need to work through a dir full of MRI dirs
% Generate the lesion area using the lesion prediction algorithm


D = dir;
D = D(~ismember({D.name}, {'.', '..'}));

% Loop thru all subdirectories
lesions = {};
% for k = 1:numel(D)
for k = 1:2
    currD = D(k).name;
    cd(currD);
    
    % load all images
    images = dir;
    images = images(~ismember({images.name}, {'.', '..'}));
    images = {images.name};
    images = arrayfun(@(x) fullfile(pwd, filesep, x), images);
    images = char(images);
    
    % convert to nii
    header_array = spm_dicom_headers(images); 
    nifti_image = spm_dicom_convert(header_array, 'all', 'flat', 'nii');
    
    % run LST algo on nii
    nifti = char(nifti_image(1).files);
    ps_LST_lpa(nifti, '');
    
    cd ..
end
end