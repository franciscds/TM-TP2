%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Trabalho de técnicas de modelagem     
%%%          Francis Santos               
%%%         2012022167                    
%       Mínimos Quadrados                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                               EXERCÍCIO 1                               %
%% ITEM 1 
% Calcula Coeficientes

matricula = 2012022167;
[b1,b2,a1,a2] = gera_coeficientes(matricula);
r = ones(1,101);
% A função é y_1(n) = -a1y(n-1)-a2y(n-2) + b1x(n-1)+b2x(n-2)
y_1(1) = 0; y_1(2) = b1;
for i=3:101
    y_1(i) = -a1*y_1(i-1)-a2*y_1(i-2) + b1*r(i-1) + b2*r(i-2);
end

figure,
subplot(2,1,1)
stem(r);
ylim([0 2])
title('Ruído Branco de entrada')
xlabel('Amostras')
ylabel('Valores das amostras')

subplot(2,1,2)
stem(1:101,y_1,'r-','linewidth',2);
xlabel('amostras')
ylabel('y[n]')
title('Resposta do sistema ao degrau')

%% ITEM 2 - Simular resposta para ruido branco
% Entrada x(n) como um ruído branco

s = rng;
L=401; %(-200 < n < 200)
x=rand(1,L);
%calculo da saída y[n]
y_2(1) = 0; y_2(2) = b1*x(1);
for i=3:length(x)
    y_2(i) = -a1*y_2(i-1)-a2*y_2(i-2) + b1*x(i-1) + b2*x(i-2);
end

%print do ruido
figure,
subplot(2,1,1)
stem(-200:200,x);
title('Ruído Branco de entrada')
xlabel('Amostras')
ylabel('Valores das amostras')

subplot(2,1,2)
stem(-200:200,y_2,'g-','linewidth',2);
title('Resposta y[n] ao ruído branco do sistema')
xlabel('Amostras')
ylabel('y[n]')

%% ITEM 3-4 - Formular problema de minimos quadrados
dados_y = y_2(51:80);
dados_r = x(51:80);
row=[2 1];
col=2:29;

% matriz de regressores
% psi = [y(2)y(1)...u(2)u(1)]
psi = [-dados_y(toeplitz(col,row)) dados_r(toeplitz(col,row))]; 

theta = pinv(psi'*psi)*psi'*dados_y(3:30)'; % parametros estimados 

disp('       Valor Estimado    Valor Real ')
disp(sprintf('  a1    %.4f           %.4f   ', theta(1), a1));
disp(sprintf('  a2     %.4f            %.4f   ', theta(2), a2));
disp(sprintf('  b1     %.4f            %.4f   ', theta(3), b1));
disp(sprintf('  b2     %.4f            %.4f   ', theta(4), b2));

%% ITEM 5 Considerando Ruido aleatorio
%SNR = mi/sigma
sigma_y = std(y_2); % desvio padrao de y

sigma_e = sigma_y/sqrt(10); %  (SNR = 10dB)

e = sigma_e.*randn(1,401); % ruido aleatorio

ye = y_2 + e; % saida ruidosa

dados_ye = ye(51:80); % alguns dados de ye

% nova matriz de regressores
psi_e = [-dados_ye(toeplitz(col,row)) dados_r(toeplitz(col,row))];

theta_e = pinv(psi_e'*psi_e)*psi_e'*dados_ye(3:30)'; % parametros estimados

disp('       Valor Estimado    Valor Real ')
disp(sprintf('  a1    %.4f           %.4f   ', theta_e(1), a1));
disp(sprintf('  a2     %.4f            %.4f   ', theta_e(2), a2));
disp(sprintf('  b1     %.4f            %.4f   ', theta_e(3), b1));
disp(sprintf('  b2     %.4f            %.4f   ', theta_e(4), b2));

%% ITEM 6 - Resposta y e curvas ajustadas
% as curvas y e y_e ajustadas pelos parametros theta e theta_e,
% respectivamente
y_ajust  = ([y_2(toeplitz(2:401,[2 1])) x(toeplitz(2:401,[2 1]))]*theta)';
ye_ajust = ([ye(toeplitz(2:401,[2 1])) x(toeplitz(2:401,[2 1]))]*theta_e)';

figure,
plot(-40:40,y_2(3:83),'b-','linewidth',2); box off;
hold on
plot(-40:40,ye(3:83),'r-','linewidth',2);
plot(-40:40,y_ajust(3:83),'g--','linewidth',2);
plot(-40:40,ye_ajust(3:83),'m--','linewidth',2);
xlim([-10 10])
xlabel('amostras')
legend('y','y_e','y_{ajustado}','y_{e_{ajustado}}')
title('Exercício 1: item 6')

%% ITEM 7 - SNR e torno de 1
s = rng;
sigma_e2 = sigma_y; % std (SNR = 1 ou 0dB)

e_2 = sigma_e2.*rand(1,401); % ruido aleatorio

yerro_2 = y_2 + e_2; % saida ruidosa

dados_ye2 = yerro_2(51:80); % alguns dados de ye_2

% nova matriz de regressores
psi_e_2 = [-dados_ye2(toeplitz(col,row)) dados_r(toeplitz(col,row))];

theta_e_2 = pinv(psi_e_2'*psi_e_2)*psi_e_2'*dados_ye2(3:30)'; % parametros estimados

% a curva y_e ajustada pelos parametros theta_e
ye_ajust_2 = [-yerro_2(toeplitz(2:401,[2 1])) x(toeplitz(2:401,[2 1]))]*theta_e_2;

disp('       Valor Estimado    Valor Real ')
disp(sprintf('  a1    %.4f           %.4f   ', theta_e_2(1), a1));
disp(sprintf('  a2     %.4f            %.4f   ', theta_e_2(2), a2));
disp(sprintf('  b1     %.4f            %.4f   ', theta_e_2(3), b1));
disp(sprintf('  b2     %.4f            %.4f   ', theta_e_2(4), b2));


figure,
plot(-40:40,y_2(3:83),'b-','linewidth',2); box off;
hold on
plot(-40:40,yerro_2(3:83),'r-','linewidth',2);
plot(-40:40,y_ajust(3:83),'b--','linewidth',2);
plot(-40:40,ye_ajust_2(3:83),'r--','linewidth',2);
xlim([-10 10])
xlabel('amostras')
legend('y','y_e','y_{ajustado}','y_{e_{ajustado}}')
title('Exercício 1: item 7')

%% ITEM 8 - tuido branco 100 pontos U2(k)
rng
k=100; %
u2=rand(1,k);
%saida do sistema a u2
y_u2(1)=0; y_u2(2)=b1*u2(1);
for i=3:length(u2)
    y_u2(i) = -a1*y_u2(i-1)-a2*y_u2(i-2) + b1*u2(i-1) + b2*u2(i-2);
end
dados_u2 = u2(51:80);
dados_yu2 = y_u2(51:80);

y_u2_estimado = [-y_u2(toeplitz(2:99,[2 1])) u2(toeplitz(2:99,[2 1]))]*theta;

% saida do sistema a entrada u_2 baseado no theta_e calculado
y_u2_estimado_2  = [-y_u2(toeplitz(2:99,[2 1])) u2(toeplitz(2:99,[2 1]))]*theta_e;

figure
plot(0:97,y_u2(3:end),'k--','linewidth',2); box off;
hold on
plot(0:97,y_u2_estimado,'b-','linewidth',2);
plot(0:97,y_u2_estimado_2,'r-','linewidth',2);
xlabel('amostras')
legend('y baseado nos coeficientes do item 1)','y baseado nos coeficientes \theta','y baseado nos coeficientes \theta_e')
title('Exercício 1: item 8')

%% ITEM 9 estimar parametros 1 e 3 ordem
% Possibilidade 1 - b2=a2=0

% estimando os parametros 
index = toeplitz(2:19,2);
index3 = dados_ye(3:20)';
% matriz de regressores
psi_9_1(:,1) = -dados_ye(index);
psi_9_1(:,2) = dados_r(index);
theta_9_1 = inv(psi_9_1'*psi_9_1)*psi_9_1'*index3; % parametros estimados 

disp('       Valor Estimado (1a ordem)  ')
disp(sprintf('  a1     %.4f               ', theta_9_1(1)));
disp(sprintf('  b1     %.4f               ', theta_9_1(2)));

% 3 ordem ( definir b3 e a3)
% estimando os parametros do modelo de 2a ordem
index2 = toeplitz(3:29,[3 2 1]);

% matriz de regressores
psi_9_2 = [-dados_ye(index2) dados_r(index2)];

theta_9_2 = pinv(psi_9_2'*psi_9_2)*psi_9_2'*dados_ye(4:30)'; % parametros estimados

fprintf('\n');

disp('       Valor Estimado (3a ordem)  ')
disp(sprintf('  a1     %.4f                   ', theta_9_2(1)));
disp(sprintf('  a2     %.4f                   ', theta_9_2(2)));
disp(sprintf('  a3     %.4f                   ', theta_9_2(3)));
disp(sprintf('  b1     %.4f                   ', theta_9_2(4)));
disp(sprintf('  b2     %.4f                     ', theta_9_2(5)));
disp(sprintf('  b3     %.4f                   ', theta_9_2(6)));

%% ITEM 10
% resposta ao degrau a partir dos parametros estimados no item (5) 
% condicoes iniciais 
y_10_1(1) = 0; y_10_1(2) = theta_e(2); % y[0] = h[0] = 0 e y[1] = h[1]u[0] = b1

% y[n] = -a1*y[n-1]-a2*y[n-2]+b1*u[n-1]+b2*u[n-2]
for i = 3:length(r)
    y_10_1(i) = -theta_e(1)*y_10_1(i-1)-theta_e(2)*y_10_1(i-2)+...
        theta_e(3)*r(i-1)+theta_e(4)*r(i-2);
end

% resposta ao degrau a partir dos parametros estimados p/ modelo de 1a
% ordem

% condicoes iniciais
 % y[0] = h[0] = 0 e y[1] = h[1]u[0] = b1
y_10_2(1) = 0; y_10_2(2) = theta_9_1(2);

% y[n] = -a1*y[n-1]+b1*u[n-1]
for i = 3:length(r)
    y_10_2(i) = -theta_9_1(1)*y_10_2(i-1)+theta_9_1(2)*r(i-1);
end

% resposta ao degrau a partir dos parametros estimados p/ modelo de 3a
% ordem
% condicoes iniciais 
y_10_3(1) = 0; y_10_3(2) = theta_9_2(4); 
y_10_3(3) = theta_9_2(4)-theta_9_2(1)*theta_9_2(4)+theta_9_2(5); % y[2] = b1 - a1b1 + b2

% y[n] = -a1*y[n-1]-a2*y[n-2]-a3*y[n-3]+b1*u[n-1]+b2*u[n-2]+b3*u[n-3]
for i = 4:length(r)
    y_10_3(i) = -theta_9_2(1)*y_10_3(i-1)-theta_9_2(2)*y_10_3(i-2)...
        -theta_9_2(3)*y_10_3(i-3)+theta_9_2(4)*r(i-1)+...
        theta_9_2(5)*r(i-2)+theta_9_2(6)*r(i-3);
end

figure
plot(0:length(r)-1,r,'b-','linewidth',2);
hold on
plot(0:length(r)-1,y_1,'k-','linewidth',2);
plot(0:length(r)-1,y_10_1,'r-','linewidth',2);
plot(0:length(r)-1,y_10_2,'m-','linewidth',2);
plot(0:length(r)-1,y_10_3,'g-','linewidth',2);
legend('Degrau de entrada', 'Resposta ao degrau do sistema verdadeiro'...
,'Resposta ao degrau do sistema utilizando os parametros do item (5)'...
,'Resposta ao degrau do sistema utilizando os parametros do modelo de 1a ordem'...
,'Resposta ao degrau do sistema utilizando os parametros do modelo de 3a ordem')
title('Exercício 1: Item 10')

%%
