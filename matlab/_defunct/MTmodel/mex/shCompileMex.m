function shCompileMex()

startDir = cd;

mexDir = which('validCorrDn3.c');
if ( isempty(dir) )
  error('Can''t find shModel MEX files in current path');
end

pos = findstr(mexDir, '/');
if ( isempty(pos) )
  error('Subdirectories appear not to be separated by slashes');
end
mexDir = mexDir(1:pos(end));


cd(mexDir);
eval(['mex -outdir ', mexDir, ' ', mexDir, 'validCorrDn3.c ', mexDir, 'svalidconvolve.c']);
eval(['mex -outdir ', mexDir, ' ', mexDir, 'dsqr.c']);
eval(['mex -outdir ', mexDir, ' ', mexDir, 'destructiveMatrixWriteAtIndices.c']);


cd(startDir);
return