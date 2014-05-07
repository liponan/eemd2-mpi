function data = bin2m(filename)

	fid = fopen(filename, 'r');
	dim = fread(fid, 1, 'int');
	
	lg = fread(fid, dim, 'int').';

	data = fread(fid, prod(lg), 'double');
	fclose(fid);

	data = reshape(data, lg);
