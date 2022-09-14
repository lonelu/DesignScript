function M = readCA(file, coorType)

if (coorType == 1)
    [status, result] = system(sprintf('grep " CA " %s | cut -c32-54 | gawk ''{print $1"\t"$2"\t"$3}''', file));
%    [status, result] = system(sprintf('grep " CA " %s | gawk ''{print $7"\t"$8"\t"$9}''', file));
    M = strread(result);
    if (size(M, 2) == 1)
        M = reshape(M, 3, size(M, 1)/3)';
    end
elseif (coorType == 0)
    M = load(file);
else
    M = file;
end
