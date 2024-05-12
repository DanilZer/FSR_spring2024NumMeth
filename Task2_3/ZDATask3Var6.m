function ZDATask3Var6
    C = 36;
    D = 1;
    N = 50;
    K = 1200;  % Графики в отсчете были получены при K~10000чтобы проверяющие смогли запустить решение за адекватное время в программе K = 1200, графики получаются менее наглядными, но общие тенденции на них все равно видны
    figure(1);
    Solve(1, C, D, 19, N, K);
    figure(2);
    Solve(1, C, D, 18, N, K);
    figure(3);
    Solve(2, C, D, 19, N, K);
    figure(4);
    Solve(2, C, D, 18, N, K);
    figure(5);
    Solve(4, C, D, 17, N, K);
    figure(6);
    Solve(4, C, D, 16, N, K);
    figure(7);
    Solve(10, C, D, 16, N, K);
    figure(8);
    Solve(10, C, D, 15, N, K);
endfunction

function Solve(a, C, D, L, N, K)
    dx = a*L/N;
    dt = min(dx^2/4/D,dx^2/4/C);
    u = zeros(K+1, N+1);
    x = linspace(-a*L/2, a*L/2, N+1);
    u(1, abs(x) <= L/2) = 1;
    C_values = C * (abs(x) <= L/2);

    g_C = C_values * dt/dx^2;

    for k = 1:K
        for i = 2:N
            u(k+1, i) = g_C(i) * u(k, i-1) + (1-2*g_C(i)+D*dt) * u(k, i) + g_C(i) * u(k, i+1);
        end
    end

    t = (0:dt:dt*K)';
    u_bar = sum(u, 2) / N;
    plot(t, u_bar);
    xlabel('t');
    ylabel('Средняя плотность');
    title(sprintf('Средняя плотность от времени (L = %f, a = %f)', L, a));
end
