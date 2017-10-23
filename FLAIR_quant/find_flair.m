function FLAIR_path = find_flair(patient_path)
% Need to work through a dir full of MRI dirs and return the path
% corresponding to the one with T2 FLAIR

files = dir(patient_path);
directoryNames = {files([files.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

% Loop thru all subdirectories
for sub_dir = directoryNames
   dir_path = fullfile(patient_dir, sub_dir);
   test_loc = dir(dir_path);
   test_img = spm_dicom_headers(fullfile(patient_dir, sub_dir, test_loc(4).name));
   
   if strcmp(test_img{1,1}.SeriesDescription, 'AX FLAIR')
        FLAIR_path = dir_path;
        break;
   end
end
end