close all;
clear all;


fprintf('********************************\n');
filename='total_dtw_paired_matrices_phonemes_negative.mat';
fprintf('Loading file: %s\n',filename);

load('total_dtw_paired_matrices_phonemes_negative.mat')



alfa = 0.42;
gamma = 0;

fftsize=640;

index_test=round(rand(1)*size(X,1));


fprintf('Phoneme: <%s>\n', cell2mat(PH(index_test)));


mcep_x=X(index_test,:);
mcep_y=Y(index_test,:);
mcep_diff=mcep_y-mcep_x;

lenvelop_x=real(mgc2sp(mcep_x,alfa, gamma, fftsize));
lenvelop_y=real(mgc2sp(mcep_y,alfa, gamma, fftsize));

figure;
hold on;
plot(lenvelop_x,'b');
plot(lenvelop_y,'r');

%dyx_env1=lenvelop_y-lenvelop_x;

dyx_env=real(mgc2sp(mcep_diff,alfa, gamma, fftsize));

%plot(dyx_env1,'g');

plot(dyx_env,'k--');


grid on;