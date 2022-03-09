figure(1),
fsurf(@(x,y) 5.*x.^2+4.*x.*y+2.*y.^2);
figure(2),
fcontour(@(x,y) 5.*x.^2+4.*x.*y+2.*y.^2);

figure(3),
fsurf(@(x,y) -1.*x.^2+4.*x.*y+2.*y.^2);
figure(4),
fcontour(@(x,y) -1.*x.^2+4.*x.*y+2.*y.^2);

figure(5),
fsurf(@(x,y) 2.*x.^2+4.*x.*y+2.*y.^2);
figure(6),
fcontour(@(x,y) 2.*x.^2+4.*x.*y+2.*y.^2);