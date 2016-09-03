%提取FEKO三角形剖分文件*.nas数据
%输入：*.nas中以GRID*开头数据——GRID.txt，*.nas中以CTRIA3开头数据——CTRIA.txt
%输出：三角形顶点坐标数据——NODE.txt，三角面片对应三顶点数据——FACE.txt
%By klingy@ahu.edu.cn

%提取FEKO剖分三角面片顶点坐标数据，格式：node x y z
fid = fopen('GRID.txt','r');
a = fscanf(fid, '%s %d %E %E %d %c %d %E', [12 inf]);
fclose(fid);
b = a([6,7,8,12],:);
fid = fopen('NODE.txt','wt');
count = fprintf(fid,'%8d % 1.9E % 1.9E % 1.9E\n',b);
fclose(fid);
%提取FEKO剖分三角面片对应顶点数据，格式：face node1 node2 node3
fid = fopen('CTRIA.txt','r');
a = fscanf(fid, '%s %d %d %d %d %d', [11 inf]);
fclose(fid);
b = a([7,9,10,11],:);
fid = fopen('FACE.txt','wt');
count = fprintf(fid,'%8d %8d %8d %8d\n',b);
fclose(fid);
