function FLAIR_path = find_flair(patient_path)
% Need to work through a dir full of MRI dirs and return the path
% corresponding to the one with T2 FLAIR

files = dir(patient_path);
directoryNames = {files([files.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

% Loop thru all subdirectories
for sub_dir = directoryNames
   dir_path = fullfile(patient_path, sub_dir);
   test_loc = dir(dir_path{1});
   test_img_path = fullfile(patient_path, sub_dir, test_loc(4).name);
   test_img = spm_dicom_headers(test_img_path{1});
   
   if strcmp(test_img{1,1}.SeriesDescription, 'AX FLAIR')
        FLAIR_path = dir_path;
        break;
   end
end
end