close all;
clear;
clc;
addpath(genpath('./'));
warning off;

params.filtVar = 1;
params.norm = 1;
params.log2 = 1;

dataDir = '/Users/eapple/Documents/MATLAB/MRGC/MRGC/Data/BRCA';
%dataDir = '/Users/eapple/Documents/testdata/';
outDir = './res';

subDir = sprintf('filtVar_%d_norm_%d_log2_%d', params.filtVar, params.norm, params.log2);
outDir = sprintf('%s/%s', outDir, subDir);
dataDir = sprintf('%s/%s', dataDir, subDir);

if ~exist(outDir, 'dir') 
    mkdir(outDir);
end

files = dir(fullfile(dataDir, '*.mat'));
%files = dir(fullfile(dataDir, '*.csv'));
files = {files.name}';
% rand('twister',5489);
for i = 1 : numel(files)
    fname = fullfile(dataDir, files{i});
	dataStr = strsplit(files{i}, '.');

	load(fname);
	%data{1} = struct2array(exp);
	%data{2} = struct2array(methy);
	%data{3} = struct2array(mirna);
    data{1} = exp;
	data{2} = methy;
	data{3} = mirna;
    
    fprintf('size(A) is %s\n', mat2str(size(data{1})))
	par.num_bases = 10 * ones(size(data,2), 1);
	par.alpha = 0.01 * ones(size(data,2), 1);
	par.beta = 0.001 * ones(size(data,2), 1);
	par.num_iters = 50;
	[idx_eg, idx_rc] = process_TCGA_datasets(data, par);
	sampleNames = fieldnames(fexp);
    
	
	if ~exist(outDir, 'dir') 
		mkdir(outDir);
	end

	outFile = sprintf('%s/%s.mat', outDir, dataStr{1});
	save(outFile, 'sampleNames', 'idx_eg', 'idx_rc');
    
    % p-values on rotational cost rc and eigen value eg
   
    X = idx_rc;
    [b,logl,H,stats] = coxphfit(X, surv,'censoring', censored);
    
    X1 = idx_eg;
    [b1,logl1,H1,stats1] = coxphfit(X1, surv,'censoring', censored);
    fprintf("%s: rc p-val= %f, eg p-val= %f",dataStr{1}, stats.p, stats1.p)
   
end
