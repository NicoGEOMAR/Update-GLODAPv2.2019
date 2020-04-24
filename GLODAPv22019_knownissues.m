%%% This file executes all changes we agreed upon for GLODAPv2.2019
%%% If you want to run this code yourself, please make sure that you cloned the entire repository AND change all "cd..." occurances 
%%% in this code according to your folder structure. Note that changes for the 74EQ cruise are not yet incooperated, as we are still waiting on replies from the PIs.

clear all
cd 'E:\'

%% Apply "general" fixes and updates
load('GLODAPv2.2019_Merged_Master_File.mat')

% CFC113
ind_1046=find(G2cruise==1046);
ind_1046_station=find(G2station(ind_1046)==31);
G2cfc113(ind_1046(ind_1046_station))=NaN;
G2cfc113f(ind_1046(ind_1046_station))=9;
G2pcfc113(ind_1046(ind_1046_station))=NaN;
clear ind*

% pH
% Cruise 1027
ind_1027=find(G2cruise==1027);
ind_1027_station_196=find(G2station(ind_1027)==196);
ind_1027_station_200=find(G2station(ind_1027)==200);
ind_1027_station_203=find(G2station(ind_1027)==203);
G2phts25p0(ind_1027(ind_1027_station_196))=NaN;
G2phts25p0(ind_1027(ind_1027_station_200))=NaN;
G2phts25p0(ind_1027(ind_1027_station_203))=NaN;
G2phtsinsitutp(ind_1027(ind_1027_station_196))=NaN;
G2phtsinsitutp(ind_1027(ind_1027_station_200))=NaN;
G2phtsinsitutp(ind_1027(ind_1027_station_203))=NaN;
G2phts25p0f(ind_1027(ind_1027_station_196))=9;
G2phts25p0f(ind_1027(ind_1027_station_200))=9;
G2phts25p0f(ind_1027(ind_1027_station_203))=9;
G2phtsinsitutpf(ind_1027(ind_1027_station_196))=9;
G2phtsinsitutpf(ind_1027(ind_1027_station_200))=9;
G2phtsinsitutpf(ind_1027(ind_1027_station_203))=9;
clear ind*

% Cruise 1039
ind_1039=find(G2cruise==1039);
ind_1039_station=find(G2station(ind_1039)==133);
ind_cast1=find(G2cast(ind_1039(ind_1039_station))==1);
ind_cast2=find(G2cast(ind_1039(ind_1039_station))==2);
ind_bottle=find(ismember(G2bottle(ind_1039(ind_1039_station(ind_cast2))),[8,9,10,11]));
G2phts25p0(ind_1039(ind_1039_station(ind_cast1)))=NaN;
G2phts25p0(ind_1039(ind_1039_station(ind_cast2(ind_bottle))))=NaN;
G2phts25p0(ind_1039(ind_1039_station(ind_cast1)))=9;
G2phts25p0(ind_1039(ind_1039_station(ind_cast2(ind_bottle))))=9;
G2phtsinsitutp(ind_1039(ind_1039_station(ind_cast1)))=NaN;
G2phtsinsitutp(ind_1039(ind_1039_station(ind_cast2(ind_bottle))))=NaN;
G2phtsinsitutpf(ind_1039(ind_1039_station(ind_cast1)))=9;
G2phtsinsitutpf(ind_1039(ind_1039_station(ind_cast2(ind_bottle))))=9;
clear ind*

% C13 
%Cruise 350
ind_350=find(G2cruise==350);
ind_350_station=find(G2station(ind_350)==7);
ind_350_station_16=find(G2station(ind_350)==16);
ind_350_station_bottle_32=find(G2bottle(ind_350(ind_350_station))==32);
ind_350_station_bottle_22=find(G2bottle(ind_350(ind_350_station_16))==22);
G2c13(ind_350(ind_350_station(ind_350_station_bottle_32)))=NaN;
G2c13f(ind_350(ind_350_station(ind_350_station_bottle_32)))=9;
G2c13(ind_350(ind_350_station(ind_350_station_bottle_22)))=NaN;
G2c13f(ind_350(ind_350_station(ind_350_station_bottle_22)))=9;
clear ind*

%Cruise 370
ind_370=find(G2cruise==370);
ind_370_station=find(G2station(ind_370)==21);
ind_cast_3=find(G2cast(ind_370(ind_370_station))==3);
ind_370_station_bottle=find(G2bottle(ind_370(ind_370_station(ind_cast_3)))==4);
G2c13(ind_370(ind_370_station(ind_cast_3(ind_370_station_bottle))))=NaN;
G2c13f(ind_370(ind_370_station(ind_cast_3(ind_370_station_bottle))))=9;
clear ind*

%Cruise 703
ind_703=find(G2cruise==703);
ind_703_station=find(G2station(ind_703)==37);
ind_703_station_bottle=find(G2bottle(ind_703(ind_703_station))==7);
G2c13(ind_703(ind_703_station(ind_703_station_bottle)))=NaN;
G2c13f(ind_703(ind_703_station(ind_703_station_bottle)))=9;
clear ind*

% Cruise 1040: Flag all variables station 52 due to bad temperature
ind_final= 1078484:1078565;
variables=who;

% Delete general variables from list, e.g. G2bottle
generals={'G2bottle','G2bottomdepth','G2cast','G2cruise','G2day','G2depth','G2hour','G2latitude','G2longitude','G2maxsampdepth','G2minute','G2month','G2pressure','G2station','G2year','expocode','expocodeno'};
for i=1:length(generals)
    ii=find(strcmp(generals{i},variables)==1);
    variables(ii)=[];
    clear ii
end
clear ii generals
for i=1:length(variables)
    if strcmp(variables{i}(end),'f')==1
        str=cat(2,variables{i},'(ind_final)', '=9;');
        eval(str)
    elseif strcmp(variables{i}(end),'qc')==1
        str=cat(2,variables{i},'(ind_final)', '=0;');
        eval(str)    
    else
        str=cat(2,variables{i},'(ind_final)', '=NaN;');
        eval(str)   
    end
end

clear ind_final i variables

% Update 33RR20160208
% Load data as structure
cd 'E:\33RR20160208'
A=load('33RR20160208.mat');
A.DELC14(A.DELC14_FLAG_W~=2 & A.DELC14_FLAG_W~=6)=NaN;
A.DELC14_FLAG_W(isnan(A.DELC14))=9;
A.C14ERR(isnan(A.DELC14))=NaN;
A.DELC13(A.DELC13_FLAG_W~=2 & A.DELC13_FLAG_W~=6)=NaN;
A.DELC13_FLAG_W(isnan(A.DELC13))=9;
A.CFC_11_FLAG_W(A.CFC_11_FLAG_W~=2 & A.CFC_11_FLAG_W~=6)=9;
A.CFC_12_FLAG_W(A.CFC_12_FLAG_W~=2 & A.CFC_12_FLAG_W~=6)=9;
A.SF6_FLAG_W(A.SF6_FLAG_W~=2 & A.SF6_FLAG_W~=6)=9;
A.CCL4_FLAG_W(A.CCL4_FLAG_W~=2 & A.CCL4_FLAG_W~=6)=9;

% Find cruise number
exno=find(strcmp(expocode,'33RR20160208')==1);
crno=expocodeno(exno);
ind=find(G2cruise==crno);
  
% Find each sample individually as file and GLODAP do not fit 100%  
Amatcher=cellstr(strcat(num2str(G2station(ind)),num2str(G2cast(ind)),num2str(G2bottle(ind))));
Bmatcher=cellstr(strcat(num2str(A.STNNBR),num2str(A.CASTNO),num2str(A.BTLNBR)));
for l=1:min(length(Amatcher),length(Bmatcher))
    final{l}=find(strcmp(Amatcher{l},Bmatcher)==1);
end
final=final';

% Assign new data
for m=1:length(final) 
    al=final{m};
    if ~isempty(al)==1
        G2c13(ind(m))=A.DELC13(final{m}(1));
        G2c13f(ind(m))=A.DELC13_FLAG_W(final{m}(1));
        G2c14(ind(m))=A.DELC14(final{m}(1));
        G2c14f(ind(m))=A.DELC14_FLAG_W(final{m}(1));
        G2c14err(ind(m))=A.C14ERR(final{m}(1));  
        if A.CFC_11_FLAG_W(final{m}(1))==9
            G2cfc11(ind(m))=NaN;
            G2cfc11f(ind(m))=9;
            G2pcfc11(ind(m))=NaN;
        end
        if A.CFC_12_FLAG_W(final{m}(1))==9
            G2cfc12(ind(m))=NaN;
            G2cfc12f(ind(m))=9;
            G2pcfc12(ind(m))=NaN;
        end
        if A.SF6_FLAG_W(final{m}(1))==9
            G2sf6(ind(m))=NaN;
            G2sf6f(ind(m))=9;
            G2psf6(ind(m))=NaN;
        end
        if A.CCL4_FLAG_W(final{m}(1))==9
            G2ccl4(ind(m))=NaN;
            G2ccl4f(ind(m))=9;
            G2pccl4(ind(m))=NaN;
        end
    end
end
    
clear final file test A ans al crno exno Amatcher Bmatcher ind l m 

% Update 33RO20110926
% Load data as structure
cd 'E:\33RO20110926'
A=load('33RO20110926.mat');
A.DELC14(A.DELC14_FLAG_W~=2 & A.DELC14_FLAG_W~=6)=NaN;
A.DELC14_FLAG_W(isnan(A.DELC14))=9;
A.C14ERR(isnan(A.DELC14))=NaN;
A.DELC13(A.DELC13_FLAG_W~=2 & A.DELC13_FLAG_W~=6)=NaN;
A.DELC13_FLAG_W(isnan(A.DELC13))=9;

% Find cruise number
exno=find(strcmp(expocode,'33RO20110926')==1);
crno=expocodeno(exno);
ind=find(G2cruise==crno);
  
% Find each sample individually as file and GLODAP do not fit 100%  
Amatcher=cellstr(strcat(num2str(G2station(ind)),num2str(G2cast(ind)),num2str(G2bottle(ind))));
Bmatcher=cellstr(strcat(num2str(A.STNNBR),num2str(A.CASTNO),num2str(A.BTLNBR)));
for l=1:min(length(Amatcher),length(Bmatcher))
    final{l}=find(strcmp(Amatcher{l},Bmatcher)==1);
end
final=final';

% Assign new data
for m=1:length(final) 
    al=final{m};
    if ~isempty(al)==1
        G2c13(ind(m))=A.DELC13(final{m}(1));
        G2c13f(ind(m))=A.DELC13_FLAG_W(final{m}(1));
        G2c14(ind(m))=A.DELC14(final{m}(1));
        G2c14f(ind(m))=A.DELC14_FLAG_W(final{m}(1));
        G2c14err(ind(m))=A.C14ERR(final{m}(1));  
    end
end
    
clear final file test A ans al crno exno Amatcher Bmatcher ind l m 

% Update 33RO20150525
% Load data as structure
cd 'E:\33RO20150525'
A=load('33RO20150525.mat');
A.DELC14(A.DELC14_FLAG_W~=2 & A.DELC14_FLAG_W~=6)=NaN;
A.DELC14_FLAG_W(isnan(A.DELC14))=9;
A.C14ERR(isnan(A.DELC14))=NaN;
A.DELC13(A.DELC13_FLAG_W~=2 & A.DELC13_FLAG_W~=6)=NaN;
A.DELC13_FLAG_W(isnan(A.DELC13))=9;

% Find cruise number
exno=find(strcmp(expocode,'33RO20150525')==1);
crno=expocodeno(exno);
ind=find(G2cruise==crno);
  
% Find each sample individually as file and GLODAP do not fit 100%  
Amatcher=cellstr(strcat(num2str(G2station(ind)),num2str(G2cast(ind)),num2str(G2bottle(ind))));
Bmatcher=cellstr(strcat(num2str(A.STNNBR),num2str(A.CASTNO),num2str(A.BTLNBR)));
for l=1:min(length(Amatcher),length(Bmatcher))
    final{l}=find(strcmp(Amatcher{l},Bmatcher)==1);
end
final=final';

% Assign new data
for m=1:length(final) 
    al=final{m};
    if ~isempty(al)==1
        G2c13(ind(m))=A.DELC13(final{m}(1));
        G2c13f(ind(m))=A.DELC13_FLAG_W(final{m}(1));
        G2c14(ind(m))=A.DELC14(final{m}(1));
        G2c14f(ind(m))=A.DELC14_FLAG_W(final{m}(1));
        G2c14err(ind(m))=A.C14ERR(final{m}(1));  
    end
end
    
clear final file test A ans al crno exno Amatcher Bmatcher ind l m 

% Update 33RO20150410
% Load data as structure
cd 'E:\33RO20150410'
A=load('33RO20150410.mat');
A.DELC14(A.DELC14_FLAG_W~=2 & A.DELC14_FLAG_W~=6)=NaN;
A.DELC14_FLAG_W(isnan(A.DELC14))=9;
A.C14ERR(isnan(A.DELC14))=NaN;
A.DELC13(A.DELC13_FLAG_W~=2 & A.DELC13_FLAG_W~=6)=NaN;
A.DELC13_FLAG_W(isnan(A.DELC13))=9;

% Find cruise number
exno=find(strcmp(expocode,'33RO20150410')==1);
crno=expocodeno(exno);
ind=find(G2cruise==crno);
  
% Find each sample individually as file and GLODAP do not fit 100%  
Amatcher=cellstr(strcat(num2str(G2station(ind)),num2str(G2cast(ind)),num2str(G2bottle(ind))));
Bmatcher=cellstr(strcat(num2str(A.STNNBR),num2str(A.CASTNO),num2str(A.BTLNBR)));
for l=1:min(length(Amatcher),length(Bmatcher))
    final{l}=find(strcmp(Amatcher{l},Bmatcher)==1);
end
final=final';

% Assign new data
for m=1:length(final) 
    al=final{m};
    if ~isempty(al)==1
        G2c13(ind(m))=A.DELC13(final{m}(1));
        G2c13f(ind(m))=A.DELC13_FLAG_W(final{m}(1));
        G2c14(ind(m))=A.DELC14(final{m}(1));
        G2c14f(ind(m))=A.DELC14_FLAG_W(final{m}(1));
        G2c14err(ind(m))=A.C14ERR(final{m}(1));  
    end
end

% Change flag 6 to 2
G2c13f(G2c13f==6)=2;
G2c14f(G2c14f==6)=2;

clear final file test A ans al crno exno Amatcher Bmatcher ind l m 

G2hour(G2hour==24)=0;

% Fix Merian expocodes
for i=1:length(expocode)
    if strcmp(expocode{i}(1:4),'06MM')==1
        expocode{i}=strcat('06M2',expocode{i}(5:end));
    end
end
clear i 


%% Apply changes agreed during RG meeting to old cruises
% Sea of Japan
ind=find(strcmp('49SH20081021',expocode)==1);
ind=expocodeno(ind);
ind=find(G2cruise==ind);
G2salinityqc(ind)=1;
G2oxygenqc(ind)=1;
G2phosphateqc(ind)=1;
G2silicateqc(ind)=1;
G2nitrateqc(ind)=1;
G2tco2qc(ind)=1;
G2tco2(ind)=G2tco2(ind)+6;
clear ind

ind=find(strcmp('49UF20121024',expocode)==1);
ind=expocodeno(ind);
ind=find(G2cruise==ind);
G2talkqc(ind)=1;
G2phtsqc(ind)=1;
G2tco2qc(ind)=1;
G2talk(ind)=G2talk(ind)+6;
clear ind

% 49NZ20050525 TA adjustment
ind=find(strcmp('49NZ20050525',expocode)==1);
ind=expocodeno(ind);
ind=find(G2cruise==ind);
G2talk(ind)=G2talk(ind)-4;
clear ind

cd 'E:\'
save GLODAPv2.2019_updated_interim1.mat
clear all


%% Apply Bob's 1st QC to 2019 JMA cruises
cd 'E:\01 Backup\Bob_changes'

% Load changes_JMA.csv as column vectors
% Initialize variables.
filename = 'E:\01 Backup\Bob_changes\changes_JMA.csv';
delimiter = ',';
startRow = 3;
formatSpec = '%s%s%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
file = dataArray{:, 1};
parameter = dataArray{:, 2};
Bob_station = dataArray{:, 3};
Bob_bottle = dataArray{:, 4};
%old_station = dataArray{:, 5};
%old_bottle = dataArray{:, 6};
%old_flag = dataArray{:, 7};
Bob_flag = dataArray{:, 8};

clearvars filename delimiter startRow formatSpec fileID dataArray ans;
clear old*

for j=1:length(parameter)
    param{j}=parameter{j}(1:end-1);
end

% Apply changes 
% Load latest GLODAP update as structure
cd 'E:\'
G2=load ('GLODAPv2.2019_updated_interim1.mat');

% Find cruise numbers of cruises which receive an update
for i=1:length(file)
    expo=find(strcmp(file{i},G2.expocode)==1);
    crsno(i)=G2.expocodeno(expo);
end
crsno=crsno';
clear expo i j

% Apply changes; Note that there is no need to include cast seperately as all castno==1
for i=1:length(crsno)
    ind=find(crsno(i)==G2.G2cruise);
    ind_station=find(Bob_station(i)==G2.G2station(ind));
    ind_bottle=find(Bob_bottle(i)==G2.G2bottle(ind(ind_station)));
    G2.(parameter{i})(ind(ind_station(ind_bottle)))=9;
    G2.(param{i})(ind(ind_station(ind_bottle)))=NaN;
end

% Save the new data - using a weird way.....
clearvars -except G2
load ('GLODAPv2.2019_updated_interim1.mat')
G2silicate=G2.G2silicate;
G2silicatef=G2.G2silicatef;
G2phosphatef=G2.G2phosphatef;
G2phosphate=G2.G2phosphate;
G2nitrate=G2.G2nitrate;
G2nitratef=G2.G2nitratef;
G2oxygen=G2.G2oxygen;
G2oxygenf=G2.G2oxygenf;
G2talk=G2.G2talk;
G2talkf=G2.G2talkf;

clear G2
cd 'E:\'
save('GLODAPV2.2019_updated_interim2.mat')


%% Fix depth issues 
for j=1:length(expocodeno)
    ind=find(G2cruise==expocodeno(j));
    egal=unique(G2station(ind)); 

    % Check for each station of each cruise   
    for i=1:length(egal)
        ind_st=find(G2station(ind)==egal(i));

        % Set bad maxsampdepth to actual max pressure of station 
        if max(G2pressure(ind(ind_st)))>max(G2maxsampdepth(ind(ind_st)))
           G2maxsampdepth(ind(ind_st))=max(G2pressure(ind(ind_st)));
        end

        % Set bad bottomdepth to actual max pressure of station + 5m 
        if max(G2pressure(ind(ind_st)))>max(G2bottomdepth(ind(ind_st)))
           G2bottomdepth(ind(ind_st))=max(G2pressure(ind(ind_st)))+5;
        end

    end
    clear ind_st ind egal
end

% Double check for bottomdepth vs maxsampdepth
clear j i egal ind
for i=1:length(G2pressure)
    if G2bottomdepth(i)<G2maxsampdepth(i)
        G2bottomdepth(i)=G2maxsampdepth(i)+5;
    end
end

clear i j
clear ans
save('GLODAPV2.2019_updated_interim3.mat')


%% Add PCO2 data!
%Create new empty columns 
G2pco=G2c13;
G2pco(:)=NaN;
G2pcof=G2pco;
G2pcotemp=G2pco;
G2pcof(:)=9;

cd 'E:\'

%Read original cruise data for cruises whcih contain discrete pco2 samples with temp
cd 'E:\01 Backup\Nico\GLODAPv2.2020\GLODAPv2_cruises\PCO2'
inputFiles = dir( fullfile('*') );
fileNames = { inputFiles.name };
fileNames(1:2)=[];

clear inputFiles

% Fill Pco2 data and flag columns with life
for kk=1:length(fileNames)
    %Load cruise data of each cruise individually and only use "good" pco2
    file=fileNames{kk}(1:12); %Cruise with Pco2 data
    cd (file)
    A=load(strcat(file,'.mat'));
    A.PCO2(A.PCO2_FLAG_W~=2 & A.PCO2_FLAG_W~=6)=NaN;
    A.PCO2_FLAG_W(isnan(A.PCO2))=9;
    A.PCO2TMP(isnan(A.PCO2))=NaN;
    
    % Cruises merged into one cruise for GLODAPv2
    if strcmp(file,'316N19920530')==1
        file='316N19920502';
    elseif strcmp(file,'316N19920713')==1
        file='316N19920502';
    end
    
    % Find the corresponding cruise number in the GLODAP file 
    exno=find(strcmp(expocode,file)==1);
    crno=expocodeno(exno);
    ind=find(G2cruise==crno);

    % Necessary excemptions as cruise files are not all in correct format
    if length(find(isnan(A.BTLNBR)))==length(A.BTLNBR) & length(find(~isnan(A.SAMPNO)))==length(A.BTLNBR)
        A.BTLNBR=A.SAMPNO;
    elseif ismember(kk,[3 9 10])==1
        A.BTLNBR=A.SAMPNO;
    elseif ismember(kk,[8])==1
        A.BTLNBR=A.SAMPNO;
        A.STNNBR(end)=100;
    elseif kk==16
        A.BTLNBR(357)=90;
    end

    % Find index of each PCO2 sample
    Amatcher=cellstr(strcat(num2str(G2station(ind)),num2str(G2cast(ind)),num2str(G2bottle(ind))));
    Bmatcher=cellstr(strcat(num2str(A.STNNBR),num2str(A.CASTNO),num2str(A.BTLNBR)));
    for l=1:min(length(Amatcher),length(Bmatcher))
        final{l}=find(strcmp(Amatcher{l},Bmatcher)==1);
    end
    final=final';

    % Assign pco2 value and flag to each sample
    for m=1:length(final) 
        al=final{m};
        if ~isempty(al)==1
            G2pco(ind(m))=A.PCO2(final{m}(1));
            G2pcotemp(ind(m))=A.PCO2TMP(final{m}(1));
            G2pcof(ind(m))=A.PCO2_FLAG_W(final{m}(1));

            % Through out values where Cruise file pressure and GLODAP pressure values aren't the same
            test(m)=G2pressure(ind(m))-A.CTDPRS(final{m}(1));
            if test(m)>0
                G2pco(ind(m))=NaN;
                G2pcotemp(ind(m))=NaN;
                G2pcof(ind(m))=9;
            end

        end
    end

    % Show files with "bad" pressure deviations
    if ~any(test==0)
        file
    end

    clear A exno crno ind file final al m l Amatcher Bmatcher
    cd ../
end

% Assign pco2 temperatures (for the cruises where the data doesn't have the needed column, i.e. are set to NaN until now) 
% First cruises with reported 20°C
cruises_20={'06MT19900123','06MT19910210','06MT19940329','316N19871123','316N19970717','316N19970815','31WT19910716','33MW19930704'};
for i=1:length(cruises_20)
    ind=find(strcmp(cruises_20{i},expocode)==1);
    ind=expocodeno(ind);
    ind=find(G2cruise==ind);

    nans=find(~isnan(G2pco(ind))==1);
    G2pcotemp(ind(nans))=20;

    clear nans ind
end
clear cruises_20 i nans ind

% Cruise with 4°C
ind=find(strcmp('320619960503',expocode)==1);
ind=expocodeno(ind);
ind=find(G2cruise==ind);

nans=find(~isnan(G2pco(ind))==1);
G2pcotemp(ind(nans))=4;

clear nans ind

%Cruises with inSitu(?) temperature
ind=find(strcmp('33LK19960415',expocode)==1);
ind=expocodeno(ind);
ind=find(G2cruise==ind);

nans=find(~isnan(G2pco(ind))==1);
G2pcotemp(ind(nans))=G2temperature(ind(nans));

clear nans ind

ind=find(strcmp('35LU19950909',expocode)==1);
ind=expocodeno(ind);
ind=find(G2cruise==ind);

nans=find(~isnan(G2pco(ind))==1);
G2pcotemp(ind(nans))=G2temperature(ind(nans));

clear nans ind ans


%% Add (and convert) fco2 cruise data
cd 'E:\01 Backup\Nico\GLODAPv2.2020\GLODAPv2_cruises\FCO2'
cruises={'31DS19940126','325020060213','32L919930718','33LG20060321','33LG20090916','33MW19920224','33MW19950922','33RO19980123','33RO20030604','33RO20070710','33RO20071215','31DS19960105','33RO20100308','33RO20131223'};

% Cruise-by-cruise
for kk=1:length(cruises)
    file=cruises{kk};

    % Clean up and load data as structure 
    clearvars -except cruises kk file G2* expocode expocodeno
    cd (file)
    A=load(strcat(file,'.mat')); 
    A.FCO2(A.FCO2_FLAG_W~=2 & A.FCO2_FLAG_W~=6)=NaN;
    A.FCO2(isnan(A.FCO2_TMP))=NaN;
    A.FCO2_FLAG_W(isnan(A.FCO2))=9;

    % Fill "missing" SiO2 and PO4 values with canyonb
    SILCAT_ca=zeros(length(A.SILCAT),1);
    PHSPHT_ca=zeros(length(A.PHSPHT),1);
    A.OXYGEN(isnan(A.OXYGEN))=A.CTDOXY(isnan(A.OXYGEN));
    A.SALNTY(isnan(A.SALNTY))=A.CTDSAL(isnan(A.SALNTY));
    data=[A.YEAR,A.LATITUDE,A.LONGITUDE,A.CTDPRS,A.CTDTMP,A.SALNTY,A.OXYGEN];
    for ii=1:length(A.SILCAT)
        if isnan(A.SILCAT(ii))==1 && isnan(A.FCO2(ii))==0 
            SILCAT_ca(ii)=silcat_nncanyonb_bit18(data(ii,:));
        else
            SILCAT_ca(ii)=A.SILCAT(ii);
        end
        if isnan(A.PHSPHT(ii))==1 && isnan(A.FCO2(ii))==0 
            PHSPHT_ca(ii)=phspht_nncanyonb_bit18(data(ii,:));
        else
            PHSPHT_ca(ii)=A.PHSPHT(ii);          
        end
    end
      
    SILCAT_ca=SILCAT_ca';
    PHSPHT_ca=PHSPHT_ca';
    clear data ii 

    % Caluclate "missing" TA values using S*67 approx.
    for k=1:length(A.ALKALI)
        if isnan(A.ALKALI(k))==1 
           TA_ca(k)=A.SALNTY(k)*67;
        else
            TA_ca(k)=A.ALKALI(k);
        end
    end

    % Convert fco2 to pco2 using CO2SYS with "GLODAP" constants and assuming 0dbar
    if length(find(~isnan(A.FCO2)))~=0 && length(find(~isnan(TA_ca)))~=0 %Only if fco2 and S data is actually present... 
       [DATA,~,~]=CO2SYS(TA_ca,A.FCO2,1,5,A.SALNTY,A.FCO2_TMP,A.FCO2_TMP,0,0,SILCAT_ca,PHSPHT_ca,1,10,1);
       A.PCO2=DATA(:,19);
       A.PCO2_FLAG_W=A.FCO2_FLAG_W;
       A.PCO2TMP=A.FCO2_TMP;
    end
    clear DATA

    % Find G2 cruise index
    exno=find(strcmp(expocode,file)==1);
    crno=expocodeno(exno);
    ind=find(G2cruise==crno);

    % Some special conversions due to bad file format
    if kk==1
        for l=1:length(A.BTLNBR)
            aholder=num2str(A.SAMPNO(l));
            A.BTLNBR(l)=str2num(aholder(2:3));
            clear aholder
        end
    end

    % Assign converted data of cruise to GLODAP
    Amatcher=cellstr(strcat(num2str(G2station(ind)),num2str(G2cast(ind)),num2str(G2bottle(ind))));
    Bmatcher=cellstr(strcat(num2str(A.STNNBR),num2str(A.CASTNO),num2str(A.BTLNBR)));
    for l=1:min(length(Amatcher),length(Bmatcher))
        final{l}=find(strcmp(Amatcher{l},Bmatcher)==1);
    end
    final=final';
    for m=1:length(final) 
        al=final{m};
        if ~isempty(al)==1
            G2pco(ind(m))=A.PCO2(final{m}(1));
            G2pcotemp(ind(m))=A.FCO2_TMP(final{m}(1));
            G2pcof(ind(m))=A.PCO2_FLAG_W(final{m}(1));
            test(m)=G2pressure(ind(m))-A.CTDPRS(final{m}(1));
            if test(m)>0
                G2pco(ind(m))=NaN;
                G2pcotemp(ind(m))=NaN;
                G2pcof(ind(m))=9;
            end
        end
    end
    clear SILCAT* PHSPHT* TA* A exno crno ind file final al m l Amatcher Bmatcher
    cd ../
end

% Change flags from 6 to 2 and 2 to 9 for values which couldn't be
% calculated due to e.g. not existing salinity values (neither bottle nor
% ctd)
G2pcof(G2pcof==6)=2;
G2pcof(isnan(G2pco))=9;

% Throw out pco2 values without given temperature
G2pcof(isnan(G2pcotemp))=9;
G2pco(isnan(G2pcotemp))=NaN;
G2pcof(isnan(G2pco))=9;

clearvars -except G2* expocode expocodeno
cd 'E:\'
save('GLODAPv2.2019_updated_interim4.mat')
clear all


%% Caluclate pco2 on constant 20°C
% Load GLODAP data as structure 
G2=load('GLODAPv2.2019_updated_interim4.mat');

% Calculate "missing" SiO2 and PO4 using CANYONB
silicate=G2.G2silicate;
phosphate=G2.G2phosphate;
Si=CANYONB(G2.G2year(isnan(silicate)),G2.G2latitude(isnan(silicate)),G2.G2longitude(isnan(silicate)),G2.G2pressure(isnan(silicate)),G2.G2temperature(isnan(silicate)),G2.G2salinity(isnan(silicate)),G2.G2oxygen(isnan(silicate)),'SiOH4');
silicate(isnan(silicate))=Si.SiOH4;
Phos=CANYONB(G2.G2year(isnan(phosphate)),G2.G2latitude(isnan(phosphate)),G2.G2longitude(isnan(phosphate)),G2.G2pressure(isnan(phosphate)),G2.G2temperature(isnan(phosphate)),G2.G2salinity(isnan(phosphate)),G2.G2oxygen(isnan(phosphate)),'PO4');
phosphate(isnan(phosphate))=Phos.PO4;
clear Phos Si

% Calculate missing TA using S*67 approx.
talk=G2.G2talk;
talk(isnan(talk))=G2.G2salinity(isnan(talk)).*67;

% Find pco2 data with "wrong" Temperature 
ind=find(~isnan(G2.G2pcotemp));
ind_temp=find(G2.G2pcotemp(ind)~=20);

% Calculate values for 20°C using TA
[DATA,~,~]=CO2SYS(talk(ind(ind_temp)),G2.G2pco(ind(ind_temp)),1,4,G2.G2salinity(ind(ind_temp)),G2.G2pcotemp(ind(ind_temp)),20,0,0,silicate(ind(ind_temp)),phosphate(ind(ind_temp)),1,10,1);
G2.G2pco(ind(ind_temp))=DATA(:,19);
G2.G2temp(ind(ind_temp))=20;

% Set flag to 9 for missing pco2 values which could not be calculated due
% to missing silicate and phosphate, i.e. missing oxygen, temperature
% etc...
G2.G2pcof(isnan(G2.G2pco))=9;

clear ind*


%% Calculate missing carbon parameters
% First search for missing carbon parameters using the "old" method, i.e. only using pH, DIC and TA
ind_misstco=find(isnan(G2.G2tco2)==1 & G2.G2phts25p0f==2 & G2.G2talkf==2 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);
ind_misstalk=find(G2.G2tco2f==2 & G2.G2phts25p0f==2 & isnan(G2.G2talk)==1 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);
ind_misspHtot=find(G2.G2tco2f==2 & isnan(G2.G2phts25p0)==1 & G2.G2talkf==2 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);

% Calculate missing parameters using CO2SYS
misstco=CO2SYS(G2.G2talk(ind_misstco),G2.G2phts25p0(ind_misstco),1,3,G2.G2salinity(ind_misstco),25,G2.G2temperature(ind_misstco),0,G2.G2pressure(ind_misstco),silicate(ind_misstco),phosphate(ind_misstco),1,10,1);
G2.G2tco2(ind_misstco)=misstco(:,2);
misstalk=CO2SYS(G2.G2tco2(ind_misstalk),G2.G2phts25p0(ind_misstalk),2,3,G2.G2salinity(ind_misstalk),25,G2.G2temperature(ind_misstalk),0,G2.G2pressure(ind_misstalk),silicate(ind_misstalk),phosphate(ind_misstalk),1,10,1);
G2.G2talk(ind_misstalk)=misstalk(:,1);
misspHtot=CO2SYS(G2.G2talk(ind_misspHtot),G2.G2tco2(ind_misspHtot),1,2,G2.G2salinity(ind_misspHtot),G2.G2temperature(ind_misspHtot),25,G2.G2pressure(ind_misspHtot),0,silicate(ind_misspHtot),phosphate(ind_misspHtot),1,10,1);
G2.G2phts25p0(ind_misspHtot)=misspHtot(:,18);

% Set flags for calculated values to 0 
G2.G2tco2f(ind_misstco)=0;
G2.G2talkf(ind_misstalk)=0;
flag0=find(~isnan(G2.G2phts25p0(ind_misspHtot)));
G2.G2phts25p0f(ind_misspHtot(flag0))=0;

% Missing pHinsitu at last as now calculated pHtot values can be used too
ind_misspHinsitu=find(isnan(G2.G2phts25p0)==0 & isnan(G2.G2phtsinsitutp)==1 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);
misspHinsitu=CO2SYS(talk(ind_misspHinsitu),G2.G2phts25p0(ind_misspHinsitu),1,3,G2.G2salinity(ind_misspHinsitu),25,G2.G2temperature(ind_misspHinsitu),0,G2.G2pressure(ind_misspHinsitu),silicate(ind_misspHinsitu),phosphate(ind_misspHinsitu),1,10,1);
G2.G2phtsinsitutp(ind_misspHinsitu)=misspHinsitu(:,18);

% Set pHinsitu flag to either 0 or 2 depending on pHtot
flag0=find(~isnan(G2.G2phtsinsitutp(ind_misspHinsitu)) & G2.G2phts25p0f(ind_misspHinsitu)~=2);
flag2=find(~isnan(G2.G2phtsinsitutp(ind_misspHinsitu)) & G2.G2phts25p0f(ind_misspHinsitu)==2);
flag9=find(G2.G2phtsinsitutpf~=9 & isnan(G2.G2phtsinsitutp)==1);
G2.G2phtsinsitutpf(ind_misspHinsitu(flag2))=2;
G2.G2phtsinsitutpf(ind_misspHinsitu(flag0))=0;
G2.G2phtsinsitutpf(flag9)=9;

clearvars -except G2 talk silicate phosphate

% Now also include new pco2 data 
% Find data where pco2 and only one further carbon parameter is present 
ind_misstco_ta=find(isnan(G2.G2tco2)==1 & G2.G2phts25p0f~=2 & G2.G2talkf==2 & G2.G2pcof==2 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);
ind_misstco_ph=find(isnan(G2.G2tco2)==1 & G2.G2phts25p0f==2 & G2.G2talkf~=2 & G2.G2pcof==2 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);

ind_missta_tco=find(isnan(G2.G2talk)==1 & G2.G2phts25p0f~=2 & G2.G2tco2f==2 & G2.G2pcof==2 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);
ind_missta_ph=find(isnan(G2.G2talk)==1 & G2.G2phts25p0f==2 & G2.G2tco2f~=2 & G2.G2pcof==2 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);

ind_missph_ta=find(isnan(G2.G2phts25p0)==1 & G2.G2tco2f~=2 & G2.G2talkf==2 & G2.G2pcof==2 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);
ind_missph_tco=find(isnan(G2.G2phts25p0)==1 & G2.G2tco2f==2 & G2.G2talkf~=2 & G2.G2pcof==2 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);

% Also find missing pco2 values, i.e samples where at least two other carbon parameters are present; Always use DIC and TA if possible
ind_misspco=find(isnan(G2.G2pco)==1 & G2.G2talkf==2 & G2.G2tco2f==2 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);
ind_misspco_ta=find(isnan(G2.G2pco)==1 & G2.G2phts25p0f==2 & G2.G2talkf==2 & G2.G2tco2f~=2 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);
ind_misspco_tco=find(isnan(G2.G2pco)==1 & G2.G2phts25p0f==2 & G2.G2talkf~=2 & G2.G2tco2f==2 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);

% Calculate missing values using CO2SYSmisstco_ta=CO2SYS(G2.G2talk(ind_misstco_ta),G2.G2pco(ind_misstco_ta),1,4,G2.G2salinity(ind_misstco_ta),20,20,0,0,silicate(ind_misstco_ta),phosphate(ind_misstco_ta),1,10,1);
% DIC & Ta not using pH
misstco_ta=CO2SYS(G2.G2talk(ind_misstco_ta),G2.G2pco(ind_misstco_ta),1,4,G2.G2salinity(ind_misstco_ta),20,G2.G2temperature(ind_misstco_ta),0,G2.G2pressure(ind_misstco_ta),silicate(ind_misstco_ta),phosphate(ind_misstco_ta),1,10,1);
G2.G2tco2(ind_misstco_ta)=misstco_ta(:,2);
missta_tco=CO2SYS(G2.G2tco2(ind_missta_tco),G2.G2pco(ind_missta_tco),2,4,G2.G2salinity(ind_missta_tco),20,G2.G2temperature(ind_missta_tco),0,G2.G2pressure(ind_missta_tco),silicate(ind_missta_tco),phosphate(ind_missta_tco),1,10,1);
G2.G2talk(ind_missta_tco)=missta_tco(:,1);

% pH
misspH_ta=CO2SYS(G2.G2talk(ind_missph_ta),G2.G2pco(ind_missph_ta),1,4,G2.G2salinity(ind_missph_ta),20,25,0,0,silicate(ind_missph_ta),phosphate(ind_missph_ta),1,10,1);
G2.G2phts25p0(ind_missph_ta)=misspH_ta(:,18);
misspH_tco=CO2SYS(G2.G2pco(ind_missph_tco),G2.G2tco2(ind_missph_tco),4,2,G2.G2salinity(ind_missph_tco),20,25,0,0,silicate(ind_missph_tco),phosphate(ind_missph_tco),1,10,1);
G2.G2phts25p0(ind_missph_tco)=misspH_tco(:,18);

% PCO2
misspco=CO2SYS(G2.G2talk(ind_misspco),G2.G2tco2(ind_misspco),1,2,G2.G2salinity(ind_misspco),G2.G2temperature(ind_misspco),20,G2.G2pressure(ind_misspco),0,silicate(ind_misspco),phosphate(ind_misspco),1,10,1);
G2.G2pco(ind_misspco)=misspco(:,19);
misspco_ta=CO2SYS(G2.G2talk(ind_misspco_ta),G2.G2phts25p0(ind_misspco_ta),1,3,G2.G2salinity(ind_misspco_ta),25,20,0,0,silicate(ind_misspco_ta),phosphate(ind_misspco_ta),1,10,1);
G2.G2pco(ind_misspco_ta)=misspco_ta(:,19);
misspco_tco=CO2SYS(G2.G2phts25p0(ind_misspco_tco),G2.G2tco2(ind_misspco_tco),3,2,G2.G2salinity(ind_misspco_tco),25,20,0,0,silicate(ind_misspco_tco),phosphate(ind_misspco_tco),1,10,1);
G2.G2pco(ind_misspco_tco)=misspco_tco(:,19);

% DIC and TA using pH
% First calculate pco2 for 25°C to agree with pH input temperature and then do regular CO2SYS caluclation
holder=CO2SYS(talk,G2.G2pco,1,4,G2.G2salinity,20,25,0,0,silicate,phosphate,1,10,1);
holder_pco=holder(:,19);
misstco_ph=CO2SYS(holder_pco(ind_misstco_ph),G2.G2phts25p0(ind_misstco_ph),4,3,G2.G2salinity(ind_misstco_ph),25,G2.G2temperature(ind_misstco_ph),0,G2.G2pressure(ind_misstco_ph),silicate(ind_misstco_ph),phosphate(ind_misstco_ph),1,10,1);
G2.G2tco2(ind_misstco_ph)=misstco_ph(:,2);    
missta_ph=CO2SYS(holder_pco(ind_missta_ph),G2.G2phts25p0(ind_missta_ph),4,3,G2.G2salinity(ind_missta_ph),25,G2.G2temperature(ind_missta_ph),0,G2.G2pressure(ind_missta_ph),silicate(ind_missta_ph),phosphate(ind_missta_ph),1,10,1);
G2.G2talk(ind_missta_ph)=missta_ph(:,1);
clear holder*

% Last but not least calculate pHinsitu from "new" pHtot
ind_misspHinsitu=find(isnan(G2.G2phts25p0)==0 & isnan(G2.G2phtsinsitutp)==1 & isnan(G2.G2salinity)==0 & isnan(G2.G2oxygen)==0);
misspHinsitu=CO2SYS(talk(ind_misspHinsitu),G2.G2phts25p0(ind_misspHinsitu),1,3,G2.G2salinity(ind_misspHinsitu),25,G2.G2temperature(ind_misspHinsitu),0,G2.G2pressure(ind_misspHinsitu),silicate(ind_misspHinsitu),phosphate(ind_misspHinsitu),1,10,1);
G2.G2phtsinsitutp(ind_misspHinsitu)=misspHinsitu(:,18);

% Assign flag 0 to calculated data
G2.G2tco2f(ind_misstco_ta)=0;
G2.G2tco2f(ind_misstco_ph)=0;
G2.G2talkf(ind_missta_tco)=0;
G2.G2talkf(ind_missta_ph)=0;
G2.G2phts25p0f(ind_missph_ta)=0;
G2.G2phts25p0f(ind_missph_tco)=0;
G2.G2pcof(ind_misspco)=0;
G2.G2pcof(ind_misspco_ta)=0;
G2.G2pcof(ind_misspco_tco)=0;
G2.G2phtsinsitutpf(ind_misspHinsitu)=0;

% Set flag to 9 for missing pco2 values which could not be calculated due
% to missing silicate and phosphate, i.e. missing temperature
% etc...
G2.G2pcof(isnan(G2.G2pco))=9;

% Save
clearvars -except G2 
load('GLODAPv2.2019_updated_interim4.mat')
G2tco2=G2.G2tco2;
G2tco2f=G2.G2tco2f;
G2talk=G2.G2talk;
G2talkf=G2.G2talkf;
G2phts25p0=G2.G2phts25p0;
G2phts25p0f=G2.G2phts25p0f;
G2phtsinsitutp=G2.G2phtsinsitutp;
G2phtsinsitutpf=G2.G2phtsinsitutpf;
G2pco=G2.G2pco;
G2pcof=G2.G2pcof;
G2pcotemp=G2.G2pcotemp;

clear G2

save('GLODAPv2.2019_updated_interim5.mat')

%% Sort product according to: 1)Cruise 2)STNNBR 3)CTDPRS 4)BOTTLE
variables={'G2cruise','G2station','G2cast','G2year','G2month','G2day','G2hour','G2minute','G2latitude','G2longitude','G2bottomdepth','G2maxsampdepth','G2bottle','G2pressure','G2depth','G2temperature','G2theta','G2salinity','G2salinityf','G2salinityqc','G2sigma0','G2sigma1','G2sigma2','G2sigma3','G2sigma4','G2gamma','G2oxygen','G2oxygenf','G2oxygenqc','G2aou','G2aouf','G2nitrate','G2nitratef','G2nitrateqc','G2nitrite','G2nitritef','G2silicate','G2silicatef','G2silicateqc','G2phosphate','G2phosphatef','G2phosphateqc','G2tco2','G2tco2f','G2tco2qc','G2talk','G2talkf','G2talkqc','G2pco','G2pcof','G2pcotemp','G2phts25p0','G2phts25p0f','G2phtsinsitutp','G2phtsinsitutpf','G2phtsqc','G2cfc11','G2pcfc11','G2cfc11f','G2cfc11qc','G2cfc12','G2pcfc12','G2cfc12f','G2cfc12qc','G2cfc113','G2pcfc113','G2cfc113f','G2cfc113qc','G2ccl4','G2pccl4','G2ccl4f','G2ccl4qc','G2sf6','G2psf6','G2sf6f','G2c13','G2c13f','G2c13qc','G2c14','G2c14f','G2c14err','G2h3','G2h3f','G2h3err','G2he3','G2he3f','G2he3err','G2he','G2hef','G2heerr','G2neon','G2neonf','G2neonerr','G2o18','G2o18f','G2toc','G2tocf','G2doc','G2docf','G2don','G2donf','G2tdn','G2tdnf','G2chla','G2chlaf'};
Matrix=[G2cruise,G2station,G2cast,G2year,G2month,G2day,G2hour,G2minute,G2latitude,G2longitude,G2bottomdepth,G2maxsampdepth,G2bottle,G2pressure,G2depth,G2temperature,G2theta,G2salinity,G2salinityf,G2salinityqc,G2sigma0,G2sigma1,G2sigma2,G2sigma3,G2sigma4,G2gamma,G2oxygen,G2oxygenf,G2oxygenqc,G2aou,G2aouf,G2nitrate,G2nitratef,G2nitrateqc,G2nitrite,G2nitritef,G2silicate,G2silicatef,G2silicateqc,G2phosphate,G2phosphatef,G2phosphateqc,G2tco2,G2tco2f,G2tco2qc,G2talk,G2talkf,G2talkqc,G2pco,G2pcof,G2pcotemp,G2phts25p0,G2phts25p0f,G2phtsinsitutp,G2phtsinsitutpf,G2phtsqc,G2cfc11,G2pcfc11,G2cfc11f,G2cfc11qc,G2cfc12,G2pcfc12,G2cfc12f,G2cfc12qc,G2cfc113,G2pcfc113,G2cfc113f,G2cfc113qc,G2ccl4,G2pccl4,G2ccl4f,G2ccl4qc,G2sf6,G2psf6,G2sf6f,G2c13,G2c13f,G2c13qc,G2c14,G2c14f,G2c14err,G2h3,G2h3f,G2h3err,G2he3,G2he3f,G2he3err,G2he,G2hef,G2heerr,G2neon,G2neonf,G2neonerr,G2o18,G2o18f,G2toc,G2tocf,G2doc,G2docf,G2don,G2donf,G2tdn,G2tdnf,G2chla,G2chlaf];
Matrix=sortrows(Matrix,[1 2 14 13]);

for i=1:length(variables)
     str=strcat(variables{i}, '=Matrix(:,i);');
     eval(str)
     clear str
end

 clear i variables Matrix str
 
 save('GLODAPv2.2019_updated_interim6.mat')

    
