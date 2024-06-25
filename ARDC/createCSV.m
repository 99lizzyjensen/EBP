%Created By: Andrew Sivaprakasam
%Last Updated: October 2023
%Description: Compile all measures into a CSV (one for Left, one for Right)
%for statistical analysis in R

addpath(pwd);
cwd = pwd();
%% Directory setup

data_dir = '/media/asivapr/AndrewNVME/Pitch_Study/All_ARDC';
output_dir = [data_dir,'/CSV_for_R'];

%%  
cd(data_dir);
filenames = dir('*.mat');
filenames = {filenames.name};

%TODO how to handle WBT??
ac_freqs = [250,500,1e3,2e3,3e3,4e3,6e3,8e3,10e3,11.2e3,12.5e3,14e3,16e3];
bc_freqs = [250,500,1e3,2e3,3e3,4e3];
dpoae_f2 = [1000,2344,3750,4781,6000,8000];
reflex_freqs = [500 1000 2000 4000];

headers = ["ID","Researcher","Time","Age",strcat("AC_",string(ac_freqs)),...
    strcat("BC_",string(bc_freqs)),strcat("DPF2_",string(dpoae_f2)),...
    strcat("REFLEX_IPSI_",string(reflex_freqs)),strcat("REFLEX_CONTRA_",string(reflex_freqs)),"QuickSIN"];

%TODO sort by date or ARDC ID...can do this in R too 

%%

%could be more efficient
csv_l = table;
csv_r = table;
for i = 1:length(filenames)
    load(filenames{i});
    id = visit.subjectID;
    res = visit.researcher;
    age = visit.Age;
    if ~isstring(res)
        res = "unknown";
    end
    
    tim = visit.time;
    
    %AUD
    %AC
    ac_r = matchAndExtract(visit.Audiogram.AC.R(:,1:2),ac_freqs);
    ac_l = matchAndExtract(visit.Audiogram.AC.L(:,1:2),ac_freqs);
    
    %BC
    bc_r = matchAndExtract(visit.Audiogram.BC.R(:,1:2),bc_freqs);
    bc_l = matchAndExtract(visit.Audiogram.BC.L(:,1:2),bc_freqs);

    %dpoae
    dp_r = matchAndExtract([visit.dpOAE.R.f2',visit.dpOAE.R.DP'],dpoae_f2);
    dp_l = matchAndExtract([visit.dpOAE.L.f2',visit.dpOAE.L.DP'],dpoae_f2);

    %reflexes (NaN for any NR or not collected)
    try
        temp = [visit.Reflexes.Frequencies',str2double(visit.Reflexes.ProbeL.Ipsi)'];
        ref_L_ipsi = matchAndExtract(temp,reflex_freqs);
    
        temp = [visit.Reflexes.Frequencies',str2double(visit.Reflexes.ProbeL.Contra)'];
        ref_L_contr = matchAndExtract(temp,reflex_freqs);
    
        temp = [visit.Reflexes.Frequencies',str2double(visit.Reflexes.ProbeR.Ipsi)'];
        ref_R_ipsi = matchAndExtract(temp,reflex_freqs);
    
        temp = [visit.Reflexes.Frequencies',str2double(visit.Reflexes.ProbeR.Contra)'];
        ref_R_contr = matchAndExtract(temp,reflex_freqs);
    catch
        ref_L_ips = nan(1,length(reflex_freqs));
        ref_R_ips = ref_L_ips;
        ref_L_contr = ref_L_ips;
        ref_R_contr = ref_L_ips;
    end

    quickSin_L = visit.QuickSIN.L;
    quickSin_R = visit.QuickSIN.R;

    dat_l = [id,res,{tim},age,num2cell(ac_l),num2cell(bc_l),num2cell(dp_l),num2cell(ref_L_ipsi),num2cell(ref_L_contr),num2cell(quickSin_L)];
    dat_r = [id,res,{tim},age,num2cell(ac_r),num2cell(bc_r),num2cell(dp_r),num2cell(ref_R_ipsi),num2cell(ref_R_contr),num2cell(quickSin_R)];
    csv_l(i,:) = cellstr(dat_l);
    csv_r(i,:) = cellstr(dat_r);

end

csv_l.Properties.VariableNames=headers;
csv_r.Properties.VariableNames=headers;

if ~isfolder(output_dir)
    mkdir(output_dir);
end

cd(output_dir);

writetable(csv_l,'all_ardc_l.csv');
writetable(csv_r,'all_ardc_r.csv');
%%
cd(cwd);