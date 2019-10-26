%% this should work to divide the results adn to get SampleType

analysis_folder = pwd
figures_folder = 'Figures_sigmoids';
imaging_string = '20X';

%%%% the best way is to vary thighs at one, so decide a sampletype and vary the timespoints (the two are indipendednt) 


sampletypes = {'_40RPM_','_40RPM2_','_CONTROL_'}   %%%% divide for mucus or not mucus
%sampletypes = {'_22.3.19','_5.4.19'}%,'_12.4.19'}   %%%% divide for mucus or not mucus

%sampletypes = {'M'}
%timepoints = {'00h';'24h';'48h';'72h'};
%timepoints = {'_T25_';'_T30_';'_T37_'};  %% different temperature
timepoints = {'_BE_','_AF_'};
%donors = {'d1';'d2';'d3'};
donors={''};
inserts = {''};
positions = {''};

% use default quick one, which is in this case the correct one
% boxsizes_vector = [16 32 48 64 96 128 160 192 224 256 340 512 1024];
boxsizes_vector = [16 32 64 128 256 512 1024];

% use the default one, works with this magnification/camera combo
q_limits_1oum = [];

flag_dryrun = false;
flag_recalculate_goodboxes = 'true';

[SampleType] = populate_DDM_Sigmoids_struct( analysis_folder, figures_folder,...
    imaging_string, sampletypes, timepoints, donors, inserts, positions,  boxsizes_vector, q_limits_1oum,...
    flag_dryrun, flag_recalculate_goodboxes);

save 'AccumData.mat';
mkdir(figures_folder)
%%  plot sigmoids
%cd '/media/np451/Seagate Backup Plus Drive1/DATA/ependymal_paper/March2019/DDM/allAnalysis'

load('AccumData_40X_12.4.19.mat',...
    'SampleType',...
    'figures_folder',...
    'sampletypes',...
    'timepoints',...
    'donors',...
    'inserts',...
    'positions');


% control first
for i = 1:numel(sampletypes)
    for j = 1:numel(timepoints)
        %MergedData(i,j) = merge_SampleType_data(SampleType, i,j,donors, '','');
        MergedData(i,j) = merge_SampleType_data(SampleType, i,j,'', '','');  %%%% i don t have different donors
    end
end

%MergeData=reshape(MergedData,[1,2]);
bsi = find((sqrt(MergedData(1).window_area_um2) / 0.146) == 64);

hf = plot_CBF_histograms_compare_sampletypes(MergedData, bsi, 0:0.5:45, 'lines',1);
 %print2svg(hf,'Figures_sigmoids\CBF_histograms(time)',1,1)



[hf,ha] = plot_sigmoid_all_errorbar_compare_sampletypes(MergedData, 'jet');
 %print2svg(hf,'Figures_sigmoids\sigmoids_compare(time)',1,1)


[hf, ha, he, hpf, hpv, hpp] = plot_sigmoid_errorbar_compare(MergedData,'lines',[],true,true,'left');
set(ha(1,:),'YLim', [0 10]);
% print2svg(hf,'Figures_sigmoids\sigmoids(time)',1,1)

%%  rename files fro movies
str2add='_10.8_v2_5_';
d=dir('*.movie');
for jj=1:numel(d)
filename=d(jj).name;
newname= strcat(filename(1:4),str2add,filename(5:end));
movefile(filename,newname)
%movefile(filename,strcat(filename(1:end-6),str2add,'.movie'));
end

