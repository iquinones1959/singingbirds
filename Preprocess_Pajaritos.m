% Manuscript Title: Humans in love are singing birds: Socially-mediated brain activity in language production
% Abbreviated Title: Humans in love are singing birds
% Authors: Clara Martin1,2*, Ileana Quiñones1* & Manuel Carreiras1,2,3
% Affiliations:
% 1Basque Center on Cognition, Brain and Language, BCBL, Donostia-San Sebastian, Spain
% 2IKERBASQUE, Basque Foundation for Science, Bilbao, Spain
% 3University of the Basque Country, UPV/EHU, Bilbao, Spain
% *equal contribution
% Corresponding authors: c.martin@bcbl.eu; i.quinones@bcbl.eu. 

% Developer: Ileana Quiñones

% Preprocessing functional MRI data 
% fMRI experiment folder "Pajaritos"
% April 2020
% Before runing you need to add SPM12 to the Matlab setpath
% as follows, addpath(genpath('/bcbl/home/public/Ileana/spm12'));
% set up SPM defaults, defaults = spm_get_defaults;
% You need to add the folder where you have the functions you are going to use
% as follows, addpath('/bcbl/home/public/Ileana/PAJARITOS');

function Preprocess_Pajaritos
path = '/bcbl/home/public/Ileana/PAJARITOS';
cd(path);
sub = dir('*'); % I run all the participants within the working directory '/bcbl/home/public/Ileana/PAJARITOS'
Slice_T = 1;
x = 256-([245,233,228,234,231,229,224,244,249,245,244,237,234,244,240,240,245,247,248,243,243,248,244,241,245,247,246,252,244,249,224,234,238,248,246,246,236,245,252,240]); % values required by the SPM function trim_img
y = 256-([78,71,72,63,60,62,53,72,72,78,76,82,69,74,81,61,76,77,75,72,75,76,66,72,88,84,75,67,90,88,59,53,73,72,73,68,66,72,63,77]);
Step1 = 1; % this variable indicates whether the segmentation of the T1 image will be done (1) or not (0)

%% Parallelizing the participants
if length(sub)>32
    parobj = parpool('ips_base',32);
else
    parobj = parpool('ips_base',length(sub));
end

%% Loop for participant
parfor nsub = 1:length(sub)
    %% Preprocessing Step 1
    sub(nsub).name
    cd([path,filesep,sub(nsub).name]);
    %-----------Trim neck (T1-T2)---------------
    t1folder = dir([path,filesep,sub(nsub).name,filesep,'*t1_mprage_sag_p2_1iso_MGH*']);
    t2folder = dir([path,filesep,sub(nsub).name,filesep,'*t2_space_sag_p2_1iso_MGH*']);
    data_name = dir([path,filesep,sub(nsub).name,filesep,'*Pajaritos_R*']);
    if Step1 == 1
        imageT1 = spm_select('FPList',[path,filesep,sub(nsub).name,filesep,t1folder(1).name],'^0.*\.nii$');
        imageT2 = spm_select('FPList',[path,filesep,sub(nsub).name,filesep,t2folder.name],'^0.*\.nii$');
        trim_img(imageT1,2,x(nsub),y(nsub));
        trim_img(imageT2,2,x(nsub),y(nsub));
        %----------Adjust commisures T1-----------
        cd([path,filesep,sub(nsub).name,filesep,t1folder.name]);
        prov = dir('tb*.nii');
        if isempty(prov)
            display([sub(nsub).name,' - Not exist']);
            continue
        else
            spm_auto_reorient(prov.name,'T1');
        end
        %----------Adjust commisures T2-----------
        cd([path,filesep,sub(nsub).name,filesep,t2folder.name]);
        prov = dir('tb*.nii');
        if isempty(prov)
            continue
        else
            spm_auto_reorient(prov.name,'T2');
        end
        %------------Preprocessing Step 2 - Batch Corregister (T1-T2)
        matlabbatch = load([path,filesep,'Template_Coregister_Reslice.mat']);
        matlabbatch = matlabbatch.matlabbatch;
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,t1folder.name],'^tb.*\.nii$'));
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,t2folder.name],'^tb.*\.nii$'));
        spm_jobman('run', matlabbatch);
        display(['Coregister Done - ' sub(nsub).name]);
        % -------------- Batch Segment (anatomical) ---------------
        matlabbatch = load ([path,filesep,'Template_SegmentT1T2_spm12.mat']);
        matlabbatch = matlabbatch.matlabbatch;
        matlabbatch{1}.spm.spatial.preproc.channel(1).vols =  cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,t1folder.name], '^tb.*.nii$'));
        matlabbatch{1}.spm.spatial.preproc.channel(2).vols =  cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,t2folder.name], '^rtb.*.nii$'));
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1]; % to save the grey matter output in native space and modulated
        matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1]; % to save the white matter output in native space and modulated
        matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1]; % to save the csf output in native space and modulated
        matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 1];
        matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1]; % to save deformations matrixes
        spm_jobman('run', matlabbatch);
        display(['Done - Suj_' sub(nsub).name]);
        % -------------- Batch Smooth (anatomical segmentations) ---------------
        matlabbatch = load ([path,filesep,'Template_Smooth.mat']);
        matlabbatch = matlabbatch.matlabbatch;
        matlabbatch{1}.spm.spatial.smooth.data = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,t1folder.name], '^wc.*tb.*.nii$'))
        matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
        spm_jobman('run', matlabbatch);
        display(['Done - Suj_' sub(nsub).name]);
        % -------------- Batch Smooth (anatomical segmentations) ---------------
        matlabbatch = load ([path,filesep,'Template_Smooth.mat']);
        matlabbatch = matlabbatch.matlabbatch;
        matlabbatch{1}.spm.spatial.smooth.data = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,t1folder.name], '^wc.*tb.*.nii$'))
        matlabbatch{1}.spm.spatial.smooth.fwhm = [10 10 10];
        matlabbatch{1}.spm.spatial.smooth.prefix = 's10';
        spm_jobman('run', matlabbatch);
        display(['Done - Suj_' sub(nsub).name]);
    end
    
    %% Preprocessing Step 2 --- Slice Timing
    if Slice_T == 1
        cd(path);
        TR = 1.8;
        Nslices = 72;
        matlabbatch = load ('Template_SliceTiming.mat');
        matlabbatch = matlabbatch.matlabbatch;
        for exp = 1:size(data_name,1)
            images = spm_select('FPList',[path,filesep,sub(nsub).name,filesep,data_name(exp).name],'^.*\.nii$');
            matlabbatch{1}.spm.temporal.st.scans{1,exp} = cellstr(images(:,:));
            matlabbatch{1}.spm.temporal.st.nslices = Nslices;
            matlabbatch{1}.spm.temporal.st.tr = TR;
            matlabbatch{1}.spm.temporal.st.ta = TR - TR/Nslices;
            matlabbatch{1}.spm.temporal.st.so = [1:2:Nslices 2:2:Nslices];
            matlabbatch{1}.spm.temporal.st.refslice = ceil(Nslices/2);
        end
        spm_jobman('run',matlabbatch);
    end
    
    %% Ajuste automatico de comisuras de las funcionales
    spmDir = which('spm'); % Definici�n de par�metros
    spmDir = spmDir(1:end-5);
    tmpl = [spmDir 'tpm' filesep 'EPI.nii'];
    vg = spm_vol(tmpl);
    flags = struct('regtype','rigid');
    for k = 1:size(data_name,1)
        vf = struct;
        cd([path,filesep,sub(nsub).name,filesep,data_name(k).name]);
        if Slice_T == 1
            prov = dir('a*.nii');
        else
            prov = dir('*.nii');
        end
        for zz = 1:size(prov,1)
            p = [prov(zz).name];
            f = strtrim(p(1,:));
            spm_smooth(f,'temp.nii',[12 12 12]);
            vf = spm_vol('temp.nii');
            [M,~] = spm_affreg(vg,vf,flags);
            M3 = M(1:3,1:3);
            [u, ~, v] = svd(M3);
            M3 = u*v';
            M(1:3,1:3) = M3;
            N = nifti(f);
            N.mat = M*N.mat;
            create(N);
        end
        delete(vf.fname);
    end
    
    %% Batch modification Realign_Unwarp
    cd(path);
    matlabbatch = load('Template_RealigUnwarp.mat');
    matlabbatch = matlabbatch.matlabbatch;
    matlabbatch = repmat(matlabbatch,1,size(data_name,1));
    for exp = 1:size(data_name,1)
        if Slice_T == 1
            matlabbatch{1,exp}.spm.spatial.realignunwarp.data(1,1).scans = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,data_name(exp).name],'^a.*\.nii$'));
        else
            matlabbatch{1,exp}.spm.spatial.realignunwarp.data(1,1).scans = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,data_name(exp).name],'.*\.nii$'));
        end
    end
    spm_jobman('run',matlabbatch);
    
    %% Batch modification Coregister
    matlabbatch = load('Template_Coregister.mat');
    matlabbatch = matlabbatch.matlabbatch;
    matlabbatch = repmat(matlabbatch,1,size(data_name,1));
    for exp = 1:size(data_name,1)
        matlabbatch{1,exp}.spm.spatial.coreg.estimate.ref = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,t1folder(1).name],'^tb.*\.nii$'));
        matlabbatch{1,exp}.spm.spatial.coreg.estimate.source = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,data_name(exp).name],'^meanu.*\.nii$'));
        matlabbatch{1,exp}.spm.spatial.coreg.estimate.other = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,data_name(exp).name],'^u.*\.nii$'));
    end
    spm_jobman('run', matlabbatch);
    
    %% Batch modification Normalise
    matlabbatch = load ('Template_NormaliseWriteT1_12.mat');
    matlabbatch = matlabbatch.matlabbatch;
    for exp = 1:size(data_name,1)
        matlabbatch{1}.spm.spatial.normalise.write.subj(1,exp).def = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,t1folder(1).name], '^y_tb.*.nii$'));
        matlabbatch{1}.spm.spatial.normalise.write.subj(1,exp).resample = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,data_name(exp).name],'^u.*\.nii$'));
        
    end
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [90,-126,-72;-90,90,108];
    spm_jobman('run', matlabbatch);
    
    %% Batch modification Detrending
    Images_det = cell(1,size(data_name,1));
    cd([path,filesep,sub(nsub).name,filesep]);
    for exp = 1:size(data_name,1)
        Images = spm_select('FPList',[path,filesep,sub(nsub).name,filesep,data_name(exp).name],'^wu.*\.nii$');
        Images_det{1,exp} = char(Images);
    end
    cspm_lmgs_2010b(Images_det);
    
    %% Batch modification Smooth
    cd(path);
    matlabbatch = load ('Template_Smooth.mat');
    matlabbatch = matlabbatch.matlabbatch;
    matlabbatch = repmat(matlabbatch,1,size(data_name,1));
    for exp = 1:size(data_name,1)
        matlabbatch{1,exp}.spm.spatial.smooth.data = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,data_name(exp).name],'^dwu.*\.nii$'));
        matlabbatch{1,exp}.spm.spatial.smooth.prefix = 's8';
    end
    spm_jobman('run', matlabbatch);
    
    %% Step 3 - Batch modification Smooth 6*6*6
    matlabbatch = load ('Template_Smooth.mat');
    matlabbatch = matlabbatch.matlabbatch;
    matlabbatch = repmat(matlabbatch,1,size(data_name,1));
    for exp = 1:size(data_name,1)
        matlabbatch{1,exp}.spm.spatial.smooth.data = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,data_name(exp).name],'^dwu.*\.nii$'));
        matlabbatch{1,exp}.spm.spatial.smooth.prefix = 's6';
        matlabbatch{1,exp}.spm.spatial.smooth.fwhm = [6 6 6];
    end
    spm_jobman('run', matlabbatch);
end

%% Close and clean
delete(parobj);
end
