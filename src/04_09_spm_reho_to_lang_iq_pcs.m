%% prep

% Clear workspace and set SPM defaults
clear; close all; clc;
spm('defaults', 'FMRI');spm_jobman('initcfg');


%% data
cd('/wdata/msmuhammad/projects/RPOE/mri/data/derivatives/func/SPM/lang-cog-pcs/features_lang_cog/ReHo');
data = readtable('spm-table.txt');

% Extract and prepare
participant_ids = data.Subj;
age = data.MRI_age;
sex = data.sex;
shared_CC1 = data.sharedLanguage_cognitionCC1;
lang_specific_PC1 = data.language_specific_PC1;
lang_specific_PC2 = data.language_specific_PC2;
cog_specific_PC1 = data.cognition_specific_PC1;
nifti_paths = data.InputFile;


% Convert to cell array if needed
if ~iscell(nifti_paths)
    nifti_paths = cellstr(nifti_paths);
end

% Make sure paths are absolute
for i = 1:length(nifti_paths)
    if ~startsWith(nifti_paths{i}, '/') && ~startsWith(nifti_paths{i}, 'C:')
        nifti_paths{i} = fullfile(pwd, nifti_paths{i});
    end
end

% Check and uncompress files
fprintf('Checking and uncompressing files...\n');
for i = 1:length(nifti_paths)
    current_file = nifti_paths{i};
    
    % Check if file exists
    if ~exist(current_file, 'file')
        error('File not found: %s', current_file);
    end
    
    % Check if compressed
    if endsWith(current_file, '.gz')
        fprintf('Uncompressing: %s\n', current_file);
        
        % Define uncompressed path
        if endsWith(current_file, '.nii.gz')
            uncompressed_file = strrep(current_file, '.nii.gz', '.nii');
        else
            uncompressed_file = strrep(current_file, '.gz', '');
        end
        
        % Check if already uncompressed
        if exist(uncompressed_file, 'file')
            fprintf('  Already uncompressed: %s\n', uncompressed_file);
        else
            % Uncompress
            try
                gunzip(current_file);
                fprintf('  Created: %s\n', uncompressed_file);
            catch
                error('Failed to uncompress: %s', current_file);
            end
        end
        
        % Update path to uncompressed version
        nifti_paths{i} = uncompressed_file;
    end
end


% Verify all files are uncompressed
fprintf('\nVerifying files...\n');
all_good = true;
for i = 1:length(nifti_paths)
    if endsWith(nifti_paths{i}, '.gz')
        fprintf('ERROR: File still compressed: %s\n', nifti_paths{i});
        all_good = false;
    elseif ~exist(nifti_paths{i}, 'file')
        fprintf('ERROR: File not found: %s\n', nifti_paths{i});
        all_good = false;
    else
        fprintf('✓ OK: %s\n', nifti_paths{i});
    end
end


if ~all_good
    error('Some files are still compressed or missing.');
end

% Prepare file paths for SPM
% SPM needs ',1' appended for 3D NIfTI files
for i = 1:length(nifti_paths)
    if ~contains(nifti_paths{i}, ',')
        nifti_paths{i} = [nifti_paths{i} ',1'];
    end
end

% Create output directory
output_dir = fullfile(pwd, 'spm_results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    cd(output_dir);
end



%% Run SPM analysis
fprintf('\nRunning SPM analysis...\n');

% Design specification
matlabbatch{1}.spm.stats.factorial_design.dir = {output_dir};
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = nifti_paths;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;

% Covariates
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).c = age;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).iCC = 5;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).c = sex;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).iCC = 5;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).c = shared_CC1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).cname = 'shared_CC1';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).iCC = 5;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(4).c = lang_specific_PC1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(4).cname = 'lang_specific_PC1';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(4).iCC = 5;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(5).c = lang_specific_PC2;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(5).cname = 'lang_specific_PC2';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(5).iCC = 5;

matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(6).c = cog_specific_PC1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(6).cname = 'cog_specific_PC1';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(6).iCC = 5;


% Other settings
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCF', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCF', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Run design specification
try
    spm_jobman('run', matlabbatch);
    fprintf('Design specification completed.\n');
catch ME
    fprintf('Error in design specification: %s\n', ME.message);
    return;
end

% Model estimation
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(output_dir, 'SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;


try
    spm_jobman('run', matlabbatch);
    fprintf('Model estimation completed.\n');
catch ME
    fprintf('Error in model estimation: %s\n', ME.message);
    return;
end


%% Diagnostic: Check regressor order
load(fullfile(output_dir, 'SPM.mat'), 'SPM');

fprintf('\n=== REGRESSOR ORDER DIAGNOSTIC ===\n');
for i = 1:length(SPM.xX.name)
    fprintf('Regressor %d: %s\n', i, SPM.xX.name{i});
end

fprintf('\nExpected order:\n');
fprintf('1. Constant/Intercept\n');
fprintf('2. age\n');
fprintf('3. sex\n');
fprintf('4. shared_CC1\n');
fprintf('5. lang_specific_PC1\n');
fprintf('6. lang_specific_PC2\n');
fprintf('7. cog_specific_PC1\n');


%% Create contrasts 
fprintf('\nCreating contrasts...\n');

matlabbatch = {};
matlabbatch{1}.spm.stats.con.spmmat = {fullfile(output_dir, 'SPM.mat')};

% T-contrasts for each variable (positive and negative)
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'shared_CC1_pos';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [0 0 0 1 0 0 0];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'shared_CC1_neg';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [0 0 0 -1 0 0 0];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'lang_specific_PC1_pos';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 0 0 1 0 0];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'lang_specific_PC1_neg';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 0 -1 0 0];
matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'lang_specific_PC2_pos';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 0 0 1 0];
matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'lang_specific_PC2_neg';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 0 -1 0];
matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'cog_specific_PC1_pos';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 0 1];
matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'cog_specific_PC1_neg';
matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 0 0 0 -1];
matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

% F-contrast for all variables
matlabbatch{1}.spm.stats.con.consess{9}.fcon.name = 'All_variables';
matlabbatch{1}.spm.stats.con.consess{9}.fcon.weights = [0 0 0 1 0 0 0; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1];
matlabbatch{1}.spm.stats.con.consess{9}.fcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 0;

try
    spm_jobman('run', matlabbatch);
    fprintf('✓ Contrasts created successfully\n');
catch ME
    fprintf('✗ Error creating contrasts: %s\n', ME.message);
end

 

%% pvalue calculation manually cuz SPM doesn't save it

fprintf('\nConverting t-maps to p-value maps...\n');

% Load SPM.mat to get degrees of freedom
load(fullfile(output_dir, 'SPM.mat'), 'SPM');
df = SPM.xX.erdf;  % Error degrees of freedom

% Process each contrast
for contrast_num = 1:8
    fprintf('Processing contrast %d/6...\n', contrast_num);
    
    % Define input/output filenames
    tmap_file = fullfile(output_dir, sprintf('spmT_%04d.nii', contrast_num));
    pmap_file = fullfile(output_dir, sprintf('p_uncorr_%04d.nii', contrast_num));
    pFDR_file = fullfile(output_dir, sprintf('p_FDR_%04d.nii', contrast_num));
    
    % Check if t-map exists
    if ~exist(tmap_file, 'file')
        fprintf('  Warning: %s not found. Skipping.\n', tmap_file);
        continue;
    end
    
    % Read t-map
    t_vol = spm_vol(tmap_file);
    t_data = spm_read_vols(t_vol);
    
    % Calculate two-tailed p-values (uncorrected)
    p_data = 2 * (1 - tcdf(abs(t_data), df));
    p_data(t_data == 0) = 1;  % Set background to p=1
    
    % Save uncorrected p-map
    p_vol = t_vol;
    p_vol.fname = pmap_file;
    p_vol.dt(1) = 16;  % Set data type to float (32-bit)
    p_vol.descrip = 'Uncorrected p-values (two-tailed)';
    spm_write_vol(p_vol, p_data);
    fprintf('  ✓ Saved: %s\n', pmap_file);
end

fprintf('\nConversion complete.\n');
fprintf('Note: These are CONTINUOUS p-value maps (0-1), not thresholded.\n');

%% Rename beta images with parameter names
load(fullfile(output_dir,'SPM.mat'),'SPM');

Vbeta = SPM.Vbeta;              % struct array of beta image handles [web:2]
reg_names = SPM.xX.name;        % cell array of regressor names [web:1]

for k = 1:numel(Vbeta)
    old_path = fullfile(output_dir, Vbeta(k).fname);

    % Clean up regressor name to make it filename‑safe
    this_name = reg_names{k};
    this_name = regexprep(this_name,'\s+','_');      % spaces -> _
    this_name = regexprep(this_name,'[^\w\-]','');   % remove weird chars

    % Build new filename, keeping the beta index if you like
    new_fname = sprintf('beta_%04d_%s.nii', k, this_name);
    new_path = fullfile(output_dir, new_fname);

    movefile(old_path, new_path);

    % Update SPM structure so future SPM calls see the new names
    Vbeta(k).fname = new_fname;
end

SPM.Vbeta = Vbeta;
save(fullfile(output_dir,'SPM.mat'),'SPM','-v7.3');

%% Rename contrast and T images with contrast names
load(fullfile(output_dir,'SPM.mat'),'SPM');

for i = 1:numel(SPM.xCon)
    cname = SPM.xCon(i).name;              % contrast name [web:11]
    cname_clean = regexprep(cname,'\s+','_');
    cname_clean = regexprep(cname_clean,'[^\w\-]','');

    % Rename con image, if present (T and F contrasts use Vcon) [web:2][web:19]
    if isfield(SPM.xCon(i),'Vcon') && ~isempty(SPM.xCon(i).Vcon)
        old_con = fullfile(output_dir, SPM.xCon(i).Vcon.fname);
        if exist(old_con,'file')
            new_con = fullfile(output_dir, sprintf('con_%04d_%s.nii', i, cname_clean));
            movefile(old_con,new_con);
            SPM.xCon(i).Vcon.fname = sprintf('con_%04d_%s.nii', i, cname_clean);
        end
    end

    % Rename statistic image, if present (T → spmT, F → spmF) [web:2][web:19]
    if isfield(SPM.xCon(i),'Vspm') && ~isempty(SPM.xCon(i).Vspm)
        old_spm = fullfile(output_dir, SPM.xCon(i).Vspm.fname);
        if exist(old_spm,'file')
            % keep prefix spmT/spmF based on STAT [web:2][web:19]
            if SPM.xCon(i).STAT == 'T'
                prefix = 'spmT';
            elseif SPM.xCon(i).STAT == 'F'
                prefix = 'spmF';
            else
                prefix = 'spmX';
            end
            new_spm = fullfile(output_dir, sprintf('%s_%04d_%s.nii', prefix, i, cname_clean));
            movefile(old_spm,new_spm);
            SPM.xCon(i).Vspm.fname = sprintf('%s_%04d_%s.nii', prefix, i, cname_clean);
        end
    end
end

save(fullfile(output_dir,'SPM.mat'),'SPM','-v7.3');

%% p-value maps with contrast names
fprintf('\nConverting t-maps to p-value maps...\n');
load(fullfile(output_dir,'SPM.mat'),'SPM');
df = SPM.xX.erdf;  % error degrees of freedom [web:18]

for contrast_num = 1:8  % only T-contrasts here
    cname = SPM.xCon(contrast_num).name;
    cname_clean = regexprep(cname,'\s+','_');
    cname_clean = regexprep(cname_clean,'[^\w\-]','');

    fprintf('Processing contrast %d (%s)...\n', contrast_num, cname);

    % If you ran the renaming above, use the new spmT prefix with name
    tmap_file = fullfile(output_dir, ...
        sprintf('spmT_%04d_%s.nii', contrast_num, cname_clean));

    if ~exist(tmap_file,'file')
        % fall back to plain spmT_XXXX if not yet renamed
        tmap_file = fullfile(output_dir, sprintf('spmT_%04d.nii', contrast_num));
    end

    if ~exist(tmap_file,'file')
        fprintf('  Warning: %s not found. Skipping.\n', tmap_file);
        continue;
    end

    pmap_file  = fullfile(output_dir, ...
        sprintf('p_uncorr_%04d_%s.nii', contrast_num, cname_clean));
    pFDR_file  = fullfile(output_dir, ...
        sprintf('p_FDR_%04d_%s.nii',   contrast_num, cname_clean)); %#ok<NASGU>

    t_vol  = spm_vol(tmap_file);
    t_data = spm_read_vols(t_vol);

    p_data = 2 * (1 - tcdf(abs(t_data), df));
    p_data(t_data == 0) = 1;

    p_vol        = t_vol;
    p_vol.fname  = pmap_file;
    p_vol.dt(1)  = 16;
    p_vol.descrip = sprintf('Uncorrected p-values (two-tailed), %s', cname);
    spm_write_vol(p_vol, p_data);
    fprintf('  ✓ Saved: %s\n', pmap_file);
end

fprintf('\nConversion complete.\n');
