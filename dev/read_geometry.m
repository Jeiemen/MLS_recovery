function geometry = read_geometry(filename)
f = fopen(filename);
raw = fread(f,inf);
strinfo = char(raw');
geometry = jsondecode(strinfo);
fclose(f);