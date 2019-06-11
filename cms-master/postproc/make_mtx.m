%% parameters
np = 10; % Number of CPU
output_dir = '../expt/expt_japan/output';

% con_file = CMS connectivity output
con_file=load('../expt/expt_japan/output/con_file_09');
% n = number of polygons
n=5;
%% 

% con_file = CMS connectivity output
con_file=load('../expt/expt_japan/output/con_file_09');
% n = number of polygons
n=5;

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

%% plot the matrix

% number of particles transiting form node i to node j
clf;
[xi,yi]=meshgrid(1:1:n, 1:1:n);
pcolor(xi,yi,settle);hold on
axis square
ylabel('Source Node');
xlabel('Receiving Node');
colorbar