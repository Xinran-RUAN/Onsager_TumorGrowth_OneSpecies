% summarize_example3_newton_iter.m
% Summarize Newton iteration counts for Example 3.

clear; clc;

outdir = 'example3';

gammaList = [10, 20, 40, 80, 160, 320];
Nfixed = 6400;

avgNewton = zeros(length(gammaList), 1);
maxNewton = zeros(length(gammaList), 1);
minNewton = zeros(length(gammaList), 1);
TuseList  = zeros(length(gammaList), 1);

for ig = 1:length(gammaList)

    gamma = gammaList(ig);

    filePattern = sprintf('example3_gamma%d_N%d*.mat', gamma, Nfixed);
    files = dir(fullfile(outdir, filePattern));

    if isempty(files)
        error('No data file found for gamma = %d and N = %d.', gamma, Nfixed);
    end

    [~, idx] = max([files.datenum]);
    matfileName = fullfile(outdir, files(idx).name);

    S = load(matfileName);

    if ~isfield(S, 'newton_iter')
        error('File %s does not contain newton_iter.', files(idx).name);
    end

    iter = S.newton_iter(:);

    avgNewton(ig) = mean(iter);
    maxNewton(ig) = max(iter);
    minNewton(ig) = min(iter);

    if isfield(S, 'time')
        TuseList(ig) = S.time(end);
    elseif isfield(S, 'T')
        TuseList(ig) = S.T;
    else
        TuseList(ig) = NaN;
    end
end

fprintf('\nNewton iteration summary, N = %d\n', Nfixed);
fprintf(' gamma      avgNewton      maxNewton      minNewton\n');
for ig = 1:length(gammaList)
    fprintf('%6d      %8.3f      %8.0f      %8.0f\n', ...
        gammaList(ig), avgNewton(ig), maxNewton(ig), minNewton(ig));
end

save(fullfile(outdir, sprintf('example3_newton_iter_summary_N%d.mat', Nfixed)), ...
    'gammaList', 'Nfixed', 'TuseList', ...
    'avgNewton', 'maxNewton', 'minNewton');