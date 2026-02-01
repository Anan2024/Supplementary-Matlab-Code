clear all
clc
warning off
%%For more information on the dataset, please refer to the website: http://archive.ics.uci.edu/ml/datasets/SECOM
addpath(genpath(pwd));
load('.\secom.data');
 load('Labels.txt');
x=secom;
x=[x,Labels];
data=x;

col_std = std(data, 0, 1, 'omitnan');  

constant_cols = col_std == 0;

data_step1 = data(:, ~constant_cols);

nan_counts = sum(isnan(data_step1 ), 1);

columns_to_remove = nan_counts > 5;

data_step2 = data_step1(:, ~columns_to_remove);

rows_with_nan = any(isnan(data_step2), 2); 

data_step3 =data_step2 (~rows_with_nan, :);

A=data_step3;

colsToKeep = true(1, size(A, 2));

B = A(:, colsToKeep);

Phase_I=B(:,249)==-1;
Phase_II=B(:,249)==1;

Phase_I_data11=B(Phase_I,1:248);
Phase_II_data11=B(Phase_II,1:248);

var_ref = var(Phase_I_data11);
tol = 1e-10;
scaling_factors = 1 ./ sqrt(max(var_ref, tol)); 
U = diag(scaling_factors);
Phase_I_data=Phase_I_data11*U;
Phase_II_data=Phase_II_data11*U;

save('Phase_I_data.mat', 'Phase_I_data');
save('Phase_II_data.mat','Phase_II_data');












