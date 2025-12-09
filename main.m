%% Clustering on the emotional speech material from CNR/Loquendo
%
% 2010, Giampiero Salvi, <giampisalvi@gmail.com>

%% load data
load data/total_dtw_paired_matrices_phonemes_negative.mat
[N, D] = size(X);
coeffs = 2:27;
D = length(coeffs);
%% get gaussian per phoneme
PHsimp = regexprep(PH, ':', '');
for h=1:length(PHsimp)
    if strcmp('M', PHsimp(h)), PHsimp{h} = 'm'; end
 %   if strcmp('Dz:', PHsimp(h)), PHsimp{h} = 'Dz'; end
 %   if strcmp('b:', PHsimp(h)), PHsimp{h} = 'b'; end
 %   if strcmp('dd:', PHsimp(h)), PHsimp{h} = 'dd'; end
 %   if strcmp('g:', PHsimp(h)), PHsimp{h} = 'g'; end
 %   if strcmp('m:', PHsimp(h)), PHsimp{h} = 'm'; end
 %   if strcmp('v:', PHsimp(h)), PHsimp{h} = 'v'; end
end
phlist = unique(PHsimp);
% eliminate silence
phlist(strcmp('#', phlist)) = [];
phlist(strcmp('@', phlist)) = [];
NPH = length(phlist);

% preallocate for speed
ns = zeros(NPH,1);
means = zeros(NPH, D);
covars = zeros(NPH, D, D);
for h=1:NPH
    smps = find(strcmp(phlist(h), PHsimp));
    ns(h) = length(smps);
    data = X(smps,coeffs)-Y(smps,coeffs);
    means(h,:) = mean(data);
    covars(h,:,:) = cov(data);
end

%% create distance matrix (compatible with pdist)
bdist = zeros(1, NPH*(NPH-1)/2);
idx = 1;
for h=1:(NPH-1)
    for k=(h+1):NPH
        bdist(idx) = bhattacharyya(means(h,:), means(k,:), squeeze(covars(h,:,:)), squeeze(covars(k,:,:)));
        idx = idx+1;
    end
end

%% perform clustering
link = linkage(bdist, 'average');
figure(1)
[h,t,perm] = dendrogram(link, 0, 'labels', phlist);
%set(gca, 'ylim', [0 15], 'fontsize', 8)
%print -depsc2 dendro_complete.eps
figure(2)
bdistsq = squareform(bdist);
imagesc(bdistsq(perm,perm))
set(gca, 'xtick', 1:length(phlist), 'xticklabel', phlist(perm))
set(gca, 'ytick', 1:length(phlist), 'yticklabel', phlist(perm))


