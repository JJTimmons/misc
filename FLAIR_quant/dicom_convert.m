function nifti_image = dicom_convert(img_path)
% Uses SPM's DICOMImport functions to convert a DICOM folder of images
% into a single nii

header_array = spm_dicom_headers(img_path); 
nifti_image = spm_dicom_convert(header_array, 'all', 'flat', 'nii');

end