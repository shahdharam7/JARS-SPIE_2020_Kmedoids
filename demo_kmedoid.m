clc;
close all;
clear all;

%% Add Path
addpath('C:\Users\Dharam\Downloads\My_Algorithms\Image');
addpath('C:\Users\Dharam\Downloads\My_Algorithms\Other Function');
addpath('C:\Users\Dharam\Downloads\My_Algorithms\M');
addpath('C:\Users\Dharam\Downloads\My_Algorithms\gt_cuprite');

%% Image Read
s=load('Cuprite.mat');
p=s.nRow;
q=s.nCol;
Bands=188;
Y=s.Y;
x=hyperConvert3d(Y,p,q,Bands);

%% Hysime
n_Hysime=12;

%VCA
[U, endmemberindex_2] = hyperVca(Y,n_Hysime);
endmemberindex_VCA=change_index(endmemberindex_2,p,q);

%kmedoid_algorithm
[endmemberindex_1] = kmedoid_algorithm(Y,n_Hysime);
endmemberindex_kmedoid_algorithm=change_index(endmemberindex_1,p,q);

%gt compare
t1=load('groundTruth_Cuprite_nEnd12.mat');
gt=t1.M;
tit=t1.cood;
n1=gt(3:103,:);
n2=gt(114:147,:);
n3=gt(168:220,:);
gt=[n1;n2;n3];

[gt_m,gt_n]=size(gt);

for i=1:gt_n
    for j=1:Bands
        extracted_kmedoid_algorithm(j,i)=x(endmemberindex_kmedoid_algorithm(i,1),endmemberindex_kmedoid_algorithm(i,2),j);
        extracted_VCA(j,i)=x(endmemberindex_VCA(i,1),endmemberindex_VCA(i,2),j);
    end
end

%% SAM Calculation
ex_n=gt_n;

for i=1:gt_n
    for j=1:ex_n
        %gt vs extraceted
        Mat_SAM_kmedoid_algorithm(i,j)=real(acos(dot(gt(:,i),extracted_kmedoid_algorithm(:,j))/(norm(gt(:,i)*norm(extracted_kmedoid_algorithm(:,j))))));
        Mat_SAM_VCA(i,j)=real(acos(dot(gt(:,i),extracted_VCA(:,j))/(norm(gt(:,i)*norm(extracted_VCA(:,j))))));
    end
end

store_kmedoid_algorithm=[0,0];
store_VCA=[0,0];

sam_kmedoid_algorithm=0;
sam_VCA=0;

sam_total_kmedoid_algorithm=0;
sam_total_VCA=0;

for i=1:gt_n
    %kmedoid_algorithm
    [max_value1,mrow]=min(Mat_SAM_kmedoid_algorithm);
    [max_value,col_kmedoid_algorithm]=min(max_value1);
    sam_total_kmedoid_algorithm=sam_total_kmedoid_algorithm+max_value;
    sam_kmedoid_algorithm=[sam_kmedoid_algorithm;max_value];
    row_kmedoid_algorithm=mrow(col_kmedoid_algorithm);
    s1=[row_kmedoid_algorithm,col_kmedoid_algorithm];
    store_kmedoid_algorithm=[store_kmedoid_algorithm;s1];
    save_kmedoid_algorithm(row_kmedoid_algorithm)=max_value;
    Mat_SAM_kmedoid_algorithm(row_kmedoid_algorithm,:)=[100*ones];
    Mat_SAM_kmedoid_algorithm(:,col_kmedoid_algorithm)=[100*ones];
    e_kmedoid_algorithm(:,row_kmedoid_algorithm)=extracted_kmedoid_algorithm(:,col_kmedoid_algorithm);
    %VCA
    [max_value1,mrow]=min(Mat_SAM_VCA);
    [max_value,col_VCA]=min(max_value1);
    sam_total_VCA=sam_total_VCA+max_value;
    sam_VCA=[sam_VCA;max_value];
    row_VCA=mrow(col_VCA);
    s1=[row_VCA,col_VCA];
    store_VCA=[store_VCA;s1];
    save_VCA(row_VCA)=max_value;
    Mat_SAM_VCA(row_VCA,:)=[100*ones];
    Mat_SAM_VCA(:,col_VCA)=[100*ones];
    e_VCA(:,row_VCA)=extracted_VCA(:,col_VCA);
end

s_E=store_kmedoid_algorithm(2:end,:);
s_CM=store_VCA(2:end,:);

for i=1:gt_n
    s_E(s_E(i,1),3)=s_E(i,2);
    s_CM(s_CM(i,1),3)=s_CM(i,2);
end

rms_sae=[rms(save_kmedoid_algorithm);
    rms(save_VCA)];

rms_sae = radtodeg(rms_sae);
disp('RMSSAE of VCA');
disp(radtodeg(rms(save_VCA)));

disp('RMSSAE of kmedoid_algorithm');
disp(radtodeg(rms(save_kmedoid_algorithm)));