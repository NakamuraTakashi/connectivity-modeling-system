close all
%% parameters
start_year = 1993;
end_year = 2018;
% n = number of polygons
n=30;

%% Load connectivity matrix data

settle2 =zeros(n,n);

for i = 1:(end_year-start_year+1)
    YEAR=num2str(start_year+(i-1));
    load(['connectivity_matrix_' YEAR '.mat'],'settle');
    settle2 = settle2 + settle;
end

settle = settle2;

%% plot the matrix

% number of particles transiting form node i to node j
% clf;
% [xi,yi]=meshgrid(1:1:n, 1:1:n);
load('MyColormaps')

f1=figure;
f1.Color=[1 1 1]; f1.Position=[0 0 580 500];
% im=image(setimg,'CDataMapping','direct');
clims = [0 2000];
im=imagesc(flipud(settle),clims);
colormap(colmap1);
hold on
% set(gca,'YDir','normal');

ylabel('Source Node');
xlabel('Sink Node');
xticks(1:n)
xticklabels(1:n)
yticks(1:n)
yticklabels(n:-1:1)
colorbar
hold off

drawnow
% hgexport(figure(1), 'output/cm_all.png', hgexport('factorystyle'),'Format','png');

%% Markov chain
% Q = settle.'; %未完成!!! 行と列を入れ替え、各列が合計1になるようにノーマライズする。
% Q = Q./sum(Q);
% 
% 定常状態になる場合？？？（未完成）
% p0 = ones(1,n);
% [X,EV]=eig(Q);
% Y=inv(X);
% W=zeros(n,n);
% for i=1:1:n
%     W(:,i) = p0*X(:,i)*Y(i,:);
% end
% Q = settle.';
Q = settle;
Q = Q./sum(Q);
Q = Q.';
Q2=Q^1000;  %!!!!!!!!!!!!!!!!! Q^n n: n-th generation
% for i=1:2
%     Q2 = Q2.'*Q;
% end
% Q2=Q2.';
f2=figure;
f2.Color=[1 1 1]; f2.Position=[0 0 580 500];
% im=image(setimg,'CDataMapping','direct');
clims = [0 0.5];
im=imagesc(flipud(Q2),clims);
colormap(colmap1);
hold on

ylabel('Site Node');
xlabel('Mixed ratio');
xticks(1:n)
xticklabels(1:n)
yticks(1:n)
yticklabels(n:-1:1)
colorbar
hold off
%% 
colors=jet(n);
f3 = figure;
axes1 = axes('Parent',f3,...
    'Position',[0.05 0.138571428571429 0.775 0.815000000000001]);
hold(axes1,'on');
b = bar(Q2,'stacked');

for k=1:n
  set(b(k),'facecolor',colors(k,:));
end
xticks(1:n)
ylim([0 1]);
legend1 = legend;
set(legend1,...
    'Position',[0.84619625740835 0.127954630229281 0.127187501132488 0.816068158404384]);

