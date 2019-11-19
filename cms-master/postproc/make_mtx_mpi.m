close all
%% parameters
np = 10; % Number of CPU
YEAR = '1993'
% output_dir = '../expt/expt_cots_2011/output';
output_dir = ['../expt/expt_cots_' YEAR '/output'];
% n = number of polygons
n=30;

%% 

for j = 1:np
    if np >= 100      
        str_file_num = num2str(j,'%03d');
    elseif np >= 10
        str_file_num = num2str(j,'%02d');
    else
        str_file_num = num2str(j,'%01d');
    end
    
    con_filename=[output_dir,'/con_file_',str_file_num];
    
    % con_file = CMS connectivity output
    con_file_part=load(con_filename);
    if j == 1
        con_file = con_file_part;
    else
        con_file = vertcat(con_file, con_file_part);
    end
end

%% compute a square connectivity matrix(ixj)

settle =zeros(n,n);

% Build the sparse matrix
for i = 1:n
    
    %recruitment polygon in column #2
    pol =con_file(con_file(:,2)==i,:);
    if isempty(pol)
        a = zeros(n,1);
    end
    
    for j=1:n
        %source polygon in column #1
        source=pol(pol(:,1)==j,1);
        a (j,:) =numel(source);
    end
    
    settle(:,i)=a;
end

% settle(settle ==0)=nan;
% save('connectivity_matrix.mat','settle');
save(['connectivity_matrix_' YEAR '.mat'],'settle');

%% plot the matrix

% number of particles transiting form node i to node j
% clf;
% [xi,yi]=meshgrid(1:1:n, 1:1:n);
load('MyColormaps')

f1=figure;
f1.Color=[1 1 1]; f1.Position=[0 0 580 500];
clims = [0 129];
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
% hgexport(figure(1), ['output/cm_' YEAR '.png'], hgexport('factorystyle'),'Format','png');

%% Markov chain
Q = settle.'; %未完成!!! 行と列を入れ替え、各列が合計1になるようにノーマライズする。
% Q = settle; %未完成!!! 行と列を入れ替え、各列が合計1になるようにノーマライズする。
Q = Q./sum(Q);

p0 = ones(1,12);
[X,EV]=eig(Q);
Y=inv(X);
W=zeros(n,n);
for i=1:1:n
    W(:,i) = p0*X(:,i)*Y(i,:);
end

figure;
im=imagesc(W);
colormap(colmap5);
hold on
% set(gca,'YDir','normal');

ylabel('Site Node');
xlabel('Mixed ratio');
colorbar
hold off
%% 
colors=jet(n);
figure;
b = bar(W,'stacked');
for k=1:n
  set(b(k),'facecolor',colors(k,:));
end
legend
