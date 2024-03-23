function ZDA_TASK1VAR3
  clc,clear;
  A_true = 2;
  B_true = -3;
  eps = 0.2;
  f = @(x) 1 ./ (A_true*x + B_true);
  n = 10;
  range = 2:20;
  perm_indices = randperm(length(range), n);
  x_data = range(perm_indices);
  x_data_2 = [2,3,4,5,6,7,8,9,10,11];
  y_data = f(x_data);
  y_data_2 = f(x_data_2);
  noise = -eps*min(abs(y_data)) + 2*eps*min(abs(y_data))*rand(1, n);
  y_data_noise = y_data + noise;
  A = [x_data'.*y_data_noise',y_data_noise'];
  noise_2 = -eps*min(abs(y_data_2)) + 2*eps*min(abs(y_data_2))*rand(1, n);
  y_data_noise_2 = y_data_2 + noise_2;
  A_2 = [x_data_2'.*y_data_noise_2',y_data_noise_2'];
  b = ones(n,1);
  [Q,R] = qr(A,0);
  [Q2,R2] = qr(A_2,0);
  paramsMNK = inv(R)*Q'*b;
  paramsMNK_2 = inv(R2)*Q2'*b;
  A_MNK = paramsMNK(1)
  B_MNK = paramsMNK'(2)
  A_MNK_2 = paramsMNK_2(1)
  B_MNK_2 = paramsMNK_2'(2)
  g = @(x) 1 ./ (A_MNK*x + B_MNK);
  g2 = @(x) 1 ./ (A_MNK_2*x + B_MNK_2);
  x_plot = linspace(2,25,1000);
  figure(1);
  plot(x_data,y_data_noise,"o",x_plot,g(x_plot));
  grid on;
  xlabel('X'); ylabel('Y');
  legend('Data with noise', 'Plot based on MNK parameters');
  figure(2);
  plot(x_data_2,y_data_noise_2,"o",x_plot,g2(x_plot));
  grid on;
  xlabel('X'); ylabel('Y');
  legend('Data with noise', 'Plot based on MNK parameters');
endfunction
