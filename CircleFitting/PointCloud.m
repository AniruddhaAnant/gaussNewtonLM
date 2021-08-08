clc;
choice=input('Enter your choice: 1 for cube 2 for sphere');
if choice == 1
    x=0;
    y=0;
    z=0;
    NumberOfPointsForCube=input('Enter number of points to be contained in and on cube');
    SideOfCube=input('Enter value of side of cube in which the points are to be contained');
    rng(0,'twister'); 
    cx = (SideOfCube-x).*rand(NumberOfPointsForCube,1) + x;
    cy = (SideOfCube-y).*rand(NumberOfPointsForCube,1) + y;
    cz = (SideOfCube-z).*rand(NumberOfPointsForCube,1) + z;
    scatter3(cx,cy,cz);
end

if choice == 2
    RadiusOfSphere=input('Enter radius of sphere in which the points are to be contained');
NumberOfPoints=input('Enter number of points to be contained in and on sphere');
rng(0,'twister');   

a = 0;
b = 2*pi;
theta = (b-a).*rand(NumberOfPoints,1) + a;
phi = (b-a).*rand(NumberOfPoints,1) + a;
j=sin(theta);
h=sin(phi);
k=cos(phi);
m=cos(theta);

hold on
for g=1:NumberOfPoints
    rx= RadiusOfSphere * j(g) * k(g);
    ry= RadiusOfSphere * j(g) * h(g);
    rz= RadiusOfSphere * m(g);
    
    scatter3(rx,ry,rz);
end
hold off
end

axis([-10 10 -10 10 -10 10])
pbaspect([1 1 1]);


 