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

function STAT_Pajaritos
%% Defining working directory
close all; clc;
path ='/bcbl/home/public/Ileana/PAJARITOS';
addpath(path);
cd(path);
sub = dir('0_*');

%% Batch for statistical design definition
for nsub = 1:length(sub)
    TR = 1.8;
    cd(path);
    matlabbatch = load ('Template_Modelo1Nivel.mat');
    matlabbatch = matlabbatch.matlabbatch;
    cd([path,filesep,sub(nsub).name,filesep]);
    data_name = dir('Pajaritos*');
    mkdir([path,filesep,sub(nsub).name,filesep,'ModeloRobusto_sm8_condfix_det_1Ses_scans_withoutM']);
    matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(spm_select('CPath',[path,filesep,sub(nsub).name,filesep,'ModeloRobusto_sm8_condfix_det_1Ses_scans_withoutM'],'.'));
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    images = {};
    for exp = 1:length(data_name)
        if exp == 1
            images = spm_select('FPList',[path,filesep,sub(nsub).name,filesep,data_name(exp).name],'^s8dwuav.*\.nii$');
        else
            images = cat(1,images,spm_select('FPList',[path,filesep,sub(nsub).name,filesep,data_name(exp).name],'^s8dwuav.*\.nii$'));
        end
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess(1,1).scans = cellstr(images);
    matlabbatch{1}.spm.stats.fmri_spec.sess(1,1).multi = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep],'.*1Ses_condfix\.mat$'));
    matlabbatch{1}.spm.stats.fmri_spec.sess(1,1).hpf = 656;
    spm_jobman('run', matlabbatch);
end

%% Batch for model estimation
for nsub = 1:length(sub)
    cd(path);
    matlabbatch = load ('Template_Estimation.mat');
    matlabbatch = matlabbatch.matlabbatch;
    matlabbatch{1,1}.spm.stats.fmri_est.spmmat = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,'ModeloRobusto_sm8_condfix_det_1Ses_scans_withoutM'],'^SPM.*\.mat$'));
    spm_jobman('run', matlabbatch);
end

%% Batch for statistical contrast definition and estimation
cd(path);
clear matlabbatch;
contrastes = load('Contrastes_Factorial.txt');
contrastes_name = textread('Contrastes_Factorial_names_add.txt','%s');
for nsub = 1:length(sub)
    matlabbatch = load ('Template_Contrast.mat');
    matlabbatch = matlabbatch.matlabbatch;
    matlabbatch{1,1}.spm.stats.con.spmmat = cellstr(spm_select('FPList',[path,filesep,sub(nsub).name,filesep,'ModeloRobusto_sm8_condfix_det_4Ses'],'^SPM.*\.mat$'));
    matlabbatch{1,1}.spm.stats.con.delete = 1;
    for cont = 1:length(contrastes_name)
        matlabbatch{1,1}.spm.stats.con.consess{cont}.tcon.name = contrastes_name{cont};
        matlabbatch{1,1}.spm.stats.con.consess{cont}.tcon.convec = contrastes(cont,:);
    end
    spm_jobman('run', matlabbatch);
end

%% Batch for 2nd level statistical design definition 
clear matlabbatch contrastes;
load ('Template_Modelo2level_Norm_none.mat');
%contrastes_name = textread('Contrastes_names.txt','%s');
load([path,filesep,sub(1).name,filesep,'ModeloRobusto_sm8_4cond_det',filesep,'SPM.mat']);
for i = 1:length(SPM.xCon)
    contrastes_name{i} = SPM.xCon(i).name;
end
matlabbatch = repmat(matlabbatch,1,length(contrastes_name));

for i = 1:length(contrastes_name)
    mkdir([path,filesep,'Modelo_Robusto_Stat_Grupo_4Ses_21suj_4cond_det_4Ses_NormNone',filesep,contrastes_name{i}]);
    direct{i} = [path,filesep,'Modelo_Robusto_Stat_Grupo_4Ses_21suj_4cond_det_4Ses_NormNone',filesep,contrastes_name{i}];
end
for ncon = 1:length(contrastes_name)
    clear contrastes
    cont_suj = 0;
    for nsub = 1:length(sub)
        clear cont
        load([path,filesep,sub(nsub).name,filesep,'ModeloRobusto_sm8_4cond_det',filesep,'SPM.mat']);
        for i = 1:length(SPM.xCon)
            cont(i) = strcmp(SPM.xCon(i).name,contrastes_name{ncon});
        end
        contrast_corresp = find(cont == 1);
        if contrast_corresp ~= 0
            cont_suj = cont_suj +1;
            if contrast_corresp < 10
                contrastes{cont_suj,1} = spm_select('FPList',[path,filesep,sub(nsub).name,filesep,'ModeloRobusto_sm8_4cond_det'],['^con_000',num2str(contrast_corresp),'.*\.img$']); %#ok<*SAGROW>
            elseif contrast_corresp >= 10
                contrastes{cont_suj,1} = spm_select('FPList',[path,filesep,sub(nsub).name,filesep,'ModeloRobusto_sm8_4cond_det'],['^con_00',num2str(contrast_corresp),'.*\.img$']);
            end
        end
    end
    matlabbatch{1,ncon}.spm.stats.factorial_design.dir = cellstr(direct{ncon});
    matlabbatch{1,ncon}.spm.stats.factorial_design.des.t1.scans = cellstr(contrastes);
end
save Modelo2level_sm8_4Ses_21suj_4cond_det_4Ses_NormNone.mat matlabbatch
spm_jobman('run', matlabbatch);

%% Batch for model estimation
clear matlabbatch;
load ('Template_est_2nivel.mat');
load([path,filesep,sub(1).name,filesep,'ModeloRobusto_sm8_4cond_det',filesep,'SPM.mat']);
for i = 1:length(SPM.xCon)
    contrastes_name{i} = SPM.xCon(i).name;
end
matlabbatch = repmat(matlabbatch,1,length(contrastes_name));
for ncon = 1: length(contrastes_name)
    matlabbatch{1,ncon}.spm.stats.fmri_est.spmmat = cellstr(spm_select('FPList',[path,filesep,'Modelo_Robusto_Stat_Grupo_4Ses_21suj_4cond_det_4Ses_NormNone',filesep,contrastes_name{ncon}],'^SPM.*\.mat$'));
end
save ModeloRobusto_est_2level_21suj_4cond_det_4Ses_NormNone.mat matlabbatch
spm_jobman('run', matlabbatch);

%% Batch for statistical contrast definition and estimation
clear matlabbatch;
load ('Template_Contrast.mat');
load([path,filesep,sub(1).name,filesep,'ModeloRobusto_sm8_4cond_det',filesep,'SPM.mat']);
for i = 1:length(SPM.xCon)
    contrastes_name{i} = SPM.xCon(i).name;
end
matlabbatch = repmat(matlabbatch,1,length(contrastes_name));
for ncon = 1: length(contrastes_name)
    matlabbatch{ncon}.spm.stats.con.spmmat = cellstr(spm_select('FPList',[path,filesep,'Modelo_Robusto_Stat_Grupo_4Ses_21suj_4cond_det_4Ses_NormNone',filesep,contrastes_name{ncon}],'^SPM.*\.mat$'));
    matlabbatch{ncon}.spm.stats.con.consess{1}.tcon.name = 'allsuj';
    matlabbatch{ncon}.spm.stats.con.consess{1}.tcon.convec = 1;
end
save ModeloRobusto_Contrastes_2level_21suj_det_4Ses_NormNone.mat matlabbatch
spm_jobman('run', matlabbatch);
