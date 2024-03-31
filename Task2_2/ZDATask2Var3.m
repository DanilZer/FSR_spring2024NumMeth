function ZDATask2Var3
  clear;clc;
  N1 = 10;
  N2 = 100;
  N3 = 1000;
  Tspan = [0,2];
  x_plot = linspace(Tspan(1), Tspan(2),100);
  True_sol = @(x) exp(sin(x));
  y_plot = True_sol(x_plot);
  [x1, y1]= Solution(N1,Tspan);
  [x2, y2]= Solution(N2,Tspan);
  [x3, y3]= Solution(N3,Tspan);
  plot(x1,y1,"r",x2,y2,"y",x3,y3,"b",x_plot,y_plot,"g")
  legend("Решение полученное методом Ньютона для N = 10","Решение полученное методом Ньютона для N = 100","Решение полученное методом Ньютона для N = 1000", "Точное решение",'Location', 'northwest')
endfunction

function [x,y] = Solution(N,Tspan)
  h = (Tspan(2) - Tspan(1))/N;
  x = linspace(Tspan(1), Tspan(2),N+1);
  y = zeros(size(x));
  dy = zeros(size(x));
  y(1) = 1;
  dy(1) = 1;

  for i = 1 : N
    [y(i + 1), dy(i + 1)] = Euler_im_Newton(y(i), dy(i), x(i), h);
  end

endfunction

function [y_next, dy_next] = Euler_im_Newton(y, dy, x, h)
  y_next = y;
  dy_next = dy;
  tol = 10e-16;
  delta = [1,1];
  while (norm(delta) > tol)
    matrix = [y_next - y - h*dy_next;
              dy_next - dy - h *(dy_next*cos(x+h) - y_next*log(y_next))];
    Jmatrix = [1 , -h;
              h*(log(y_next)+1), 1 - h*cos(x + h)];
    delta = - Jmatrix\matrix;
    y_next = y_next + delta(1);
    dy_next = dy_next + delta(2);
  endwhile
endfunction
