function [I1,m1]=gauss3(r,winsz,O,T)

m1=[cos(T(1))  -sin(T(1)) 0
    sin(T(1))  cos(T(1))  0
    0           0           1];
m2=[cos(T(2)) -sin(T(2)) 0
    0               0       1
    sin(T(2))   cos(T(2)) 0
    ];
m3=[  0               0       1
    cos(T(3)) -sin(T(3)) 0
    sin(T(3))   cos(T(3)) 0
    ];
m=m1*m2*m3;
[X, Y,Z] = meshgrid(-winsz:winsz, -winsz:winsz,-winsz:winsz);
% I1 = exp( -1/2* (( X-O(1)).^2/r(1)^2 + (Y-O(2)).^2/r(2)^2+(Z-O(3)).^2/r(3)^2)) ;
I1 = exp( -1/2* (( (X-O(1))*m(1,1)+(Y-O(2))*m(1,2)+(Z-O(3))*m(1,3)).^2/r(1)^2 + ((X-O(1))*m(2,1)+(Y-O(2))*m(2,2)+(Z-O(3))*m(2,3)).^2/r(2)^2+((X-O(1))*m(3,1)+(Y-O(2))*m(3,2)+(Z-O(3))*m(3,3)).^2/r(3)^2)) ;
m1=m(1,:);
if m1(1)<0
    m1 = -m1;
end
end