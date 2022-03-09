b = [1,2,1,1.414, 1, 2];
b1 = [1,2,1,1.4142, 1, 2];
x1 = invhilb(6)*b';
x2 = invhilb(6)*b1';
c = cond(hilb(6));