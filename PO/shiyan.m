X0 = load('NODE.txt');
tri = load('FACE.txt');
X0(:,1)=[];
tri(:,1)=[];
save('E:\\pm\\X0.mat','X0');
save('E:\\pm\\tri.mat','tri');
%
m = size(tri,1);
for iii=1:m;
    xx = [X0(tri(iii,1),1);X0(tri(iii,2),1);X0(tri(iii,3),1);X0(tri(iii,1),1)];
    yy = [X0(tri(iii,1),2);X0(tri(iii,2),2);X0(tri(iii,3),2);X0(tri(iii,1),2)];
    zz = [X0(tri(iii,1),3);X0(tri(iii,2),3);X0(tri(iii,3),3);X0(tri(iii,1),3)];
    plot3(xx,yy,zz);
    hold on 
end