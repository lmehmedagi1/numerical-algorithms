function x = Newton_Raphson(f, x0, eps, maxiters)
	
	if (eps<=0 )
		error('Invalid parameters!');
		return;
	end

	x = x0(:);
	for i = 0 : maxiters

		[J, H] = numderivative(f, x, [], [], 'hypermat');		

		if (norm(J) < eps)
			return;
		end;

		x = x - H\J.';
	end	
endfunction
















function[x] = BFGS(f, x0, eps, max_iters)
	
	if (eps<=0 )
		error('Invalid parameters!');
		return;
	end

	n = length(x0);	
	H = eye(n, n);
	x = x0(:);	
				
	dx = numderivative(f, x);

	for i = 1 : max_iters
				
		if (norm(dx) <= eps)
			return;
		end

		p = - H * dx.';
		
		h = strong_wolfe_powell(f, x, f(x), numderivative(f, x), p);
			
		x_new = x + h*p;
		
		dx_old = dx;
		dx = numderivative(f, x_new);

		if(isnan(dx))
			mprintf('NaN value detected');
			return;
		end
		
		s_i = x_new - x;
		y_i = dx.' - dx_old.';	
		rho_i = 1/(y_i.'*s_i);
		
		if(i==1)	
		   H = (y_i.'*s_i)/(y_i.'*y_i) * eye(n,n);
		end

		H = (eye(n,n) - rho_i * s_i * y_i')*H*(eye(n,n) - rho_i * y_i * s_i.') + rho_i * s_i * s_i.';

		x = x_new;
	end

	disp('Minimum cannot be found within given number of iterations or it does not exist!');
endfunction









function [h] = strong_wolfe_powell(f, x0, f0, g0, p)
c1 = 1e-4;
c2 = 0.9;
h_min = 0.25;
h_max = 2.5;
h = 1;
h_im1 = 0;
h_i = h;
f_im1 = f0;
dphi0 = g0*p;
i = 1;
max_iters = 20;

    while 1
  	x = x0 + h_i*p;
      
	f_i = f(x);
      	g_i = numderivative(f, x);
	
	if (f_i > f0 + c1*h_i*dphi0) || ((f_i >= f_im1) && (i>1))
    		h = zoom(f, x0, f0, g0, p, h_im1, h_i);
    		break;
      	end

        dphi = g_i*p;
 	
	if ( abs(dphi) <= -c2*dphi0 )
    		h = h_i;
    		break;
  	end
      
  	if ( dphi >= 0 )
    		h = zoom(f, x0, f0, g0, p, h_i, h_im1);
    		break;
  	end

	h_im1 = h_i;
  	f_im1 = f_i;

  	h_i = h_i + 0.8*(h_max - h_i);

	if (i > max_iters)
    		h = h_i;
    		break;
  	end
  
  	i = i + 1;
    end
endfunction










function [h] = zoom(f, x0, f0, g0, p, h_low, h_high)
c1 = 1e-4;
c2 = 0.9;
i = 1;
max_iters = 20;
dphi0 = g0*p;

    while 1
	h_i = (h_low + h_high)/2;	
	
	x = x0 + h_i*p;
       
	f_i = f(x);
       	g_i = numderivative(f, x);
	
  	x_low = x0 + h_low*p;
  	f_low = f(x_low); 

  	if ( (f_i > f0 + c1*h_i*dphi0) || ( f_i >= f_low))
    		h_high = h_i;
  	else
		dphi = g_i*p;
    	
		if (abs(dphi) <= -c2*dphi0)
      			h = h_i;
      			return;
    		end
    
		if (dphi * (h_high - h_low) >= 0)
      			h_high = h_low;
    		end

            	h_low = h_i;
  	end

	if (i > max_iters)
		h = h_i;
    		break;
  	end
	
	i = i + 1;
    end
endfunction














function y = parabola1(x)
	y = x^2 + 4*x - 5;
endfunction;

function y = f1(x)
	y = sin(x);
endfunction;

function y = f2(x)
	y = cos(3*%pi*x)/x;
endfunction;

function y = f3(x)
	y = x^3 + 3*x^2 - 2*x + 1;
endfunction;



function y = rosenbrock(x)
	y = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
endfunction; 

function y = rosenbrock2(x1, x2)
	y = 100*(x2 - x1.^2).^2 + (1 - x1).^2;
endfunction;


 

function y = rastrigin(x)
	y = 20 + x(1).^2 + x(2).^2 -10*(cos(2*%pi*x(1)) + cos(2*%pi*x(2)));
endfunction;

function y = rastrigin2(x1, x2)
	y = 20 + x1.^2 + x2.^2 -10*(cos(2*%pi*x1) + cos(2*%pi*x2));
endfunction;





function [y] = levy(xx)
	d = length(xx);

	for ii = 1:d
		w(ii) = 1 + (xx(ii) - 1)/4;
	end

	term1 = (sin(%pi*w(1)))^2;
	term3 = (w(d)-1)^2 * (1+(sin(2*%pi*w(d)))^2);

	sum = 0;
	for ii = 1:(d-1)
		wi = w(ii);
	        new = (wi-1)^2 * (1+10*(sin(%pi*wi+1))^2);
		sum = sum + new;
	end

	y = term1 + sum + term3;
endfunction

function [y] = levy2(x1, x2)

	w1 = 1 + (x1 - 1)./4;
	w2 = 1 + (x2 - 1)./4;

	term1 = (sin(%pi*w1)).^2;
	term3 = (w2-1).^2 .* (1+(sin(2*%pi*w2)).^2);

	wi = w1;
	new = (wi-1).^2 .* (1+10*(sin(%pi*wi+1)).^2);

	y = term1 + new + term3;
endfunction





function [y] = ackley(xx, a, b, c)

	d = length(xx);

	if (nargin < 4)
  	  	c = 2*%pi;
	end
	if (nargin < 3)
   		b = 0.2;
	end
	if (nargin < 2)
    		a = 20;
	end

	sum1 = 0;
	sum2 = 0;
	for ii = 1:d
		xi = xx(ii);
		sum1 = sum1 + xi^2;
		sum2 = sum2 + cos(c*xi);
	end

	term1 = -a * exp(-b*sqrt(sum1/d));
	term2 = -exp(sum2/d);

	y = term1 + term2 + a + exp(1);	
endfunction

function [y] = ackley2(x1, x2)

	a = 20;
	b = 0.2;
	c = 2*%pi;

	sum1 = 0;
	sum2 = 0;

	xi = x1;
	sum1 = sum1 + xi.^2;
	sum2 = sum2 + cos(c*xi);

	xi = x2;
	sum1 = sum1 + xi.^2;
	sum2 = sum2 + cos(c*xi);

	term1 = -a * exp(-b*sqrt(sum1/2));
	term2 = -exp(sum2/2);

	y = term1 + term2 + a + exp(1);
endfunction




function [y] = spherefmod(xx)
	d = length(xx);
	sum = 0;
	for ii = 1:d
		xi = xx(ii);
		sum = sum + xi^2;
	end
	
	y = sum;	
endfunction

function [y] = spherefmod2(x1, x2)
	y = x1.^2 + x2.^2;
endfunction





function [y] = himmelblau(xx)
	x1 = xx(1);
	x2 = xx(2);

	y = (x1^2 + x2 -11)^2 + (x1 + x2^2 -7)^2;
endfunction

function [y] = himmelblau2(x1, x2)
	y = (x1.^2 + x2 -11).^2 + (x1 + x2.^2 -7).^2;
endfunction






function [y] = hosaki(x)
	y = (1-8*x(1) + 7*x(1)^2 - 7/3*x(1)^3 + 1/4*x(1)^4) * x(2)^2*exp(-x(1));
endfunction;

function [y] = hosaki2(x1, x2)
	y = (1-8*x1 + 7*x1.^2 - 7/3*x1.^3 + 1/4*x1.^4) .* x2.^2.*exp(-x2);
endfunction;







function [y] = holdertable(x)
	y = -abs(exp(abs(1-sqrt(x(1)^2+x(2)^2)/%pi))*sin(x(1))*cos(x(2)));
endfunction;

function [y] = holdertable2(x1, x2)
	y = -abs(exp(abs(1-sqrt(x1.^2+x2.^2)/%pi)).*sin(x1).*cos(x2));
endfunction;






function [y] = crossit(xx)
	x1 = xx(1);
	x2 = xx(2);

	fact1 = sin(x1)*sin(x2);
	fact2 = exp(abs(100 - sqrt(x1^2+x2^2)/%pi));

	y = -0.00001 * (abs(fact1*fact2)+1)^0.1;
endfunction


function [y] = crossit2(x1, x2)
	fact1 = sin(x1).*sin(x2);
	fact2 = exp(abs(100 - sqrt(x1.^2+x2.^2)/%pi));

	y = -0.00001 * (abs(fact1.*fact2)+1).^0.1;
endfunction






function [y] = camelback(x)
	y = (4-2.1*x(1)^2+x(1)^4/3)*x(1)^2 + x(1)*x(2) + (-4+4*x(2)^2)*x(2)^2;
endfunction;

function [y] = camelback2(x1, x2)
	y = (4-2.1*x1.^2+x1.^4/3).*x1.^2 + x1.*x2 + (-4+4*x2.^2).*x2.^2;
endfunction;






function [y] = easom(xx)
	x1 = xx(1);
	x2 = xx(2);

	fact1 = -cos(x1)*cos(x2);
	fact2 = exp(-(x1-%pi)^2-(x2-%pi)^2);

	y = fact1*fact2;
endfunction

function [y] = easom2(x1, x2)

	fact1 = -cos(x1).*cos(x2);
	fact2 = exp(-(x1-%pi).^2-(x2-%pi).^2);

	y = fact1.*fact2;
endfunction






function [y] = langer(xx)
	d = 2;
	m = 5;
	c = [1, 2, 5, 2, 3];
	A = [3, 5; 5, 2; 2, 1; 1, 4; 7, 9]

	outer = 0;
	for ii = 1:m
	    inner = 0;
	    for jj = 1:d
	        xj = xx(jj);
	        Aij = A(ii,jj);
	        inner = inner + (xj-Aij)^2;
	    end
	    new = c(ii) * exp(-inner/%pi) * cos(%pi*inner);
	    outer = outer + new;
	end
	y = outer;
endfunction

function [y] = langer2(x1, x2)
	d = 2;
	m = 5;
	c = [1, 2, 5, 2, 3];
	A = [3, 5; 5, 2; 2, 1; 1, 4; 7, 9]

	outer = 0;
	for ii = 1:m
	    inner = 0;
	
	    Aij = A(ii,1);
	    inner = inner + (x1-Aij).^2;
	
	    Aij = A(ii,2);
	    inner = inner + (x2-Aij).^2;		

	    new = c(ii) .* exp(-inner/%pi) .* cos(%pi.*inner);
	    outer = outer + new;
	end
	y = outer;
endfunction







function test(f, x0, eps, maxiters, xmin, xmax, tacke, stringf)
	x = BFGS(f, x0, eps, maxiters);
	interval = linspace(xmin, xmax, tacke);
	plot(interval,f);
	plot(x, f(x), 'r*');
	l2 = 'Nadjeni minimum: ' + string(x);
	title('Funkcija: ' + stringf)
	xlabel(l2);
endfunction



function test3d(f, x0, eps, maxiters, xmin, xmax, tacke, naslov, f2, ugao1, ugao2)
	x = BFGS(f, x0, eps, maxiters);
	
	x1 = linspace(xmin, xmax, tacke);
	x2 = linspace(xmin, xmax, tacke);
	[XX1, XX2] = meshgrid(x1, x2);
	
	Z = f2(XX1, XX2);
	
	f = gcf();
	f.color_map = jetcolormap(25);
	plot3d1(x1, x2, Z);
	
	a = get("current_axes");
	a.box = "back_half";
	a.grid = [-1, -1, 1];
	a.grid_position = "background";
	a.grid_thickness = [6, 6, 6];
	a.grid_style = [7, 7, 7];
	a.rotation_angles = [ugao1, ugao2];

	t = a.title;
	t.text = naslov;

	l = a.x_label;
	l.text = "Nadjeni minimum: " + string(x);
endfunction



subplot(221);
test(parabola1, 13, 0.00001, 100, -6, 2, 100, 'x^2 + 4*x - 5');

subplot(222);
test(f1, -2, 0.00001, 100, -6, 2, 100, 'sin(x)');

subplot(223);
test(f2, -2, 0.00001, 100, -6, -1, 1500, 'cos(3*pi*x)/x');

subplot(224);
test(f3, -2, 0.00001, 100, -6, 2, 100, 'x^3 + 3*x^2 - 2*x + 1');



scf;



subplot(231);
test3d(rastrigin, [-1; 2], 0.00001, 100, -5.12, 5.12, 100, "Rastrigin funkcija", rastrigin2, 10, 130);

subplot(233);
test3d(ackley, [1;3], 0.00001, 100, -20, 20, 100, "Ackley funkcija", ackley2, 70, 80);

subplot(232);
test3d(rosenbrock, [10;7], 0.00001, 100, -5, 5, 100, "Rosenbrock funkcija", rosenbrock2, 15, 50);

subplot(234);
test3d(levy, [3;-12], 0.00001, 100, -10, 10, 100, "Levy funkcija", levy2, 10, 130);

subplot(235);
test3d(spherefmod, [1;3], 0.00001, 100, -5.12, 5.12, 100, "Sphere funkcija", spherefmod2, 20, 130);

subplot(236);
test3d(himmelblau, [-10;88], 0.00001, 100, -4, 4, 100, "Himmelblau funkcija", himmelblau2, 45, 45);




scf;




subplot(231);
test3d(hosaki, [3;0], 0.00001, 100, 0, 5, 100, "Hosaki funkcija", hosaki2, 75, 25);

subplot(232);
test3d(holdertable, [2;3], 0.00001, 100, -10, 10, 100, "Holdertable funkcija", holdertable2, 45, 45);

subplot(233);
test3d(crossit, [-9;7], 0.0000001, 100, -10, 10, 100, "Cross-in-Tray funkcija", crossit2, 75, 70);

subplot(234);
test3d(camelback, [-1;-0.5], 0.00001, 100, -1, 1, 100, "Camelback funkcija", camelback2, 40, 50);

subplot(235);
test3d(easom, [3.5;2.1], 0.00001, 100, 0, 6.28, 100, "Easom funkcija", easom2, 80, 140);

subplot(236);
test3d(langer, [7;6], 0.00001, 100, 0, 10, 100, "Langermann funkcija", langer2, 50, 40);


