function m2bin_eemd(F, filename)

fid = fopen(filename, 'w');
dim = length(size(F));
fwrite(fid, dim, 'int');
fwrite(fid, size(F), 'int');
fwrite(fid, F(:), 'double');
fclose(fid);

end