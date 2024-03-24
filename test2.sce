// un exemple simple
a=[1,7;2,8];b=[10;11];
deff('[y]=fsol1(x)','y=a*x+b');
deff('[y]=fsolj1(x)','y=a');
[xres]=fsolve([100;100],fsol1);
a*xres+b
[xres]=fsolve([100;100],fsol1,fsolj1);
disp(a*xres+b);

// voir SCI/modules/optimization/sci_gateway/fortran/Ex-fsolve.f
[xres]=fsolve([100;100],'fsol1','fsolj1',1.e-7);
disp(a*xres+b);