clear all;
close all;
clc;

filepath = 'mat/';
filenames = {'human'        % Average human
             'cat'          % Mixed species cat
             'roundworm'    % C elegans, the only complete connectome
             'fly'          % Medulla of drosophila melanogaster
             'macaque'      % Rhesus ape
             'mouse'        % mouse
             'rat'};        % Rattus norvegicus
filesuffix = '.mat';

numfiles = numel(filenames);

stats = zeros(numfiles,6);

for file_index=1:numfiles
    % Read adjacency matrix from file
    file = join({filepath,filenames{file_index},filesuffix},'');
    A = dlmread(file);
    
    % Make A symmetric
    A = max(A,A');
    A = min(A,1);
    
    % Number of nodes
    N = size(A,1);
    
    % Compute eigenpairs (unordered)
    [eigs_vec,eigs] = eig(A);
    
    % Extract eigenvalues
    eigs = diag(eigs);
    
    % Sort by eigenvalues (min -> max)
    [eigs_sorted,ind_sorted] = sort(eigs);
    eigs_vec = eigs_vec(:,ind_sorted);
    
    % Spectral gap (largest - 2nd largest eigenvalues)
    spectral_gap = eigs_sorted(N) - eigs_sorted(N-1); 
    
    % Get principal eigenvector (<=> max eigenvalue)
    principal_eigenvector = eigs_vec(:,N);
    
    % Subgraph centrality (sum of principal_eigenvector(j)^2)
    SC = zeros(N,1);
    for i = 1:N
        SC_i = 0;
        for j = 1:N
            lambda_j = eigs_sorted(j);
            u_ij = eigs_vec(i,j);
            SC_i = SC_i + u_ij.^2.*sinh(lambda_j);
        end
        SC(i) = SC_i;
    end
    
    A_ = sinh(eigs_sorted(N)).^(-1/2);
    xi_A = log(A_);
    diff_vector = log(abs(principal_eigenvector)) - ... 
                        (xi_A + 0.5*log(SC));
    xi_vector = ( diff_vector ).^2;
    xi = sqrt(mean(xi_vector));
    
    % Make plots
    fig = figure;
    
    subplot(2,2,1);
    G = graph(A);
    plot(G);
    title('Network')
   
    subplot(2,2,2);
    stem(eigs_sorted);
    ylabel('Eigenvalue');
    title(sprintf('Spectral gap: %d',spectral_gap))
    
    subplot(2,2,3);
    spy(A);
    title('Adjacency matrix')
    
    subplot(2,2,4);
    loglog(SC,abs(principal_eigenvector),'*');
    xlabel('Subgraph centrality')
    ylabel('Principal eigenvector')
    title(sprintf('Xi: %d',xi));
    
    % Add reference line with slope 1/2
    hold on;
    startpoint = min(SC);
    endpoint = max(SC);
    loglog([startpoint, endpoint], ...
           [sqrt(startpoint),sqrt(endpoint)] ./ sqrt(mean(SC)));
    
    title_name = [upper(filenames{file_index}(1)), ...
                  lower(filenames{file_index}(2:end))];
    p = mtit(fig,title_name,'fontsize',14);
    
    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits',...
            'Inches','PaperSize',[pos(3), pos(4)]);
    print(fig,title_name,'-dpdf','-r0');
    
    stats(file_index,:) = ...
        [N min(sum(A)) mean(sum(A)) max(sum(A)) spectral_gap, xi];
end

% Print stats and make plots to compare expansion and
% subgraph centrality between species
stats

%% Unused plots below
%figure;
%stem(stats(:,1:3))

%figure;
%stem(log(stats(:,4:5)))
