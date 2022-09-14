% addpath 'E:\GitHub_Design\cccp'
% folder = 'E:\DesignData\ligands\CoiledCoil\C4_all_CYS\output_x40_y20\'
% filepaths = strcat(folder, '*.txt')
% cd 'E:\DesignData\ligands\CoiledCoil\C4_all_CYS\output_x40_y20\'

addpath 'E:\GitHub_Design\cccp'
folder = 'E:\DesignData\ligands\LigandBB\_reverse_design\c2_coil\_z_fix_c2\'
filepaths = strcat(folder, '*.txt')
cd 'E:\DesignData\ligands\LigandBB\_reverse_design\c2_coil\_z_fix_c2\'

listing = dir(filepaths)

parpool

parfor i=1:length(listing)
    baseFileName = listing(i).name
    fullFileName = fullfile(folder, baseFileName)
    %fprintf(1, fullFileName)
    try
        fcrick(fullFileName, 4, 'GENERAL', 0, baseFileName, [], [], [], []) 
    catch
        warning('Crashing...');
    end
end

