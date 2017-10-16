%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                               EXERCÍCIO 3                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Trabalho de técnicas de modelagem     
%%%         Francis Santos               
%%%         2012022167                    
%       Mínimos Quadrados                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IITE 1 
clear all;
[x, y] = simula_sistema(2012022167);
%plot sinais de entrada e saída
figure
subplot(2,1,1)
plot(x,'b-','linewidth',2);
box off;
ylabel('x(t)')
title(sprintf('Dados de entrada e saída do sistema\n(a)'))
subplot(2,1,2)
plot(y,'b-','linewidth',2);
box off;
xlabel('t(s)')
ylabel('y(t)')
title('(b)')
%% ITEM A
% calculo da FAC de y 
[acf,lags,bounds] = autocorr(y,200);
[acf_min, idmin] = min(acf); % primeiro minimo

figure;
plot(lags,acf,'b-','linewidth',2);
hold on
plot(lags,bounds(1)*ones(length(lags)),'k--','linewidth',0.05)
plot(lags,bounds(2)*ones(length(lags)),'k--','linewidth',0.05)
xlabel('Atraso')
ylabel('Correlação')
title('Autocorrelação do sinal de saída')

%% ITEM B
[FCC,lags] = xcov(x,y);
FCC_media = mean(FCC); % media da FCC entre x e y
FCC_sd = std(FCC); % desvio padrao da FCC entre x e y

figure 
plot(lags,FCC,'b-','linewidth',2);
hold on 
plot(lags,(FCC_media-2*FCC_sd)*ones(1,length(lags)),'k--');
plot(lags,(FCC_media+2*FCC_sd)*ones(1,length(lags)),'k--');
xlabel('Atrasos')
title('Cross-Correlation entre os sinais x(k) e y(k)') 

%% ITEM C
% tamanho para ID
N = length(x)/2;
x_identificacao = x(1:N);
y_identificacao = y(1:N);
% dados para validacao

x_validacao = x(N+1:end);
y_validacao = y(N+1:end);

%% ITEM 2 - ESTRUTURA DO MODELO

% matriz de regressores ordem 1
psi1(:,1) = y_identificacao(toeplitz(1:N-1,1));
psi1(:,2) = x_identificacao(toeplitz(1:N-1,1)); 
%parametros estimados (a1 e b1)
theta1 = pinv(psi1'*psi1)*psi1'*y_identificacao(2:N);  
% vetor de residuos
sigma1 = y_identificacao(2:N)-psi1*theta1;
% criterio de informacao de Akaike (ordem 1)
AIC_1 = N*log(sigma1.^2)+2*length(theta1);

% matriz de regressores ordem 2
psi2 = [y_identificacao(toeplitz(2:N-1,[2 1])) x_identificacao(toeplitz(2:N-1,[2 1]))];
% parametros estimados (a1, a2, b1 e b2)
theta2 = pinv(psi2'*psi2)*psi2'*y_identificacao(3:N);  
% vetor de residuos
sigma2 = y_identificacao(3:N)-psi2*theta2;
% criterio de informacao de Akaike (ordem 2)
AIC_2 = N*log(sigma2.^2)+2*length(theta2);

% matriz de regressores ordem 3
psi3 = [y_identificacao(toeplitz(3:N-1,[3 2 1])) x_identificacao(toeplitz(3:N-1,[3 2 1]))];
% parametros estimados a1, a2, a3, b1, b2 e b3
theta3 = pinv(psi3'*psi3)*psi3'*y_identificacao(4:N);  
% vetor de residuos
sigma3 = y_identificacao(4:N)-psi3*theta3;
% criterio de informacao de Akaike (ordem 3)
AIC_3 = N*log(sigma3.^2)+2*length(theta3);


% matriz de regressores ordem 4
psi4 = [y_identificacao(toeplitz(4:N-1,[4 3 2 1])) x_identificacao(toeplitz(4:N-1,[4 3 2 1]))];
% parametros estimados a1, a2, a3, a4, b1, b2, b3 e b4
theta4 = pinv(psi4'*psi4)*psi4'*y_identificacao(5:N);  
% vetor de residuos
sigma4 = y_identificacao(5:N)-psi4*theta4;
% criterio de informacao de Akaike (ordem 4)
AIC_4 = N*log(sigma4.^2)+2*length(theta4);

% matriz de regressores ordem 5
psi5 = [y_identificacao(toeplitz(5:N-1,[5 4 3 2 1]))...
    x_identificacao(toeplitz(5:N-1,[5 4 3 2 1]))];
% parametros estimados a1, a2, a3, a4, a5, b1, b2, b3, b4 e b5
theta5 = pinv(psi5'*psi5)*psi5'*y_identificacao(6:N);  
% vetor de residuos
sigma5 = y_identificacao(6:N)-psi5*theta5;
% criterio de informacao de Akaike (ordem 5)
AIC_5 = N*log(sigma5.^2)+2*length(theta5);

% matriz de regressores ordem 6
psi6 = [y_identificacao(toeplitz(6:N-1,[6 5 4 3 2 1]))...
    x_identificacao(toeplitz(6:N-1,[6 5 4 3 2 1]))];
% parametros estimados a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5 e b6
theta6 = pinv(psi6'*psi6)*psi6'*y_identificacao(7:N);  
% vetor de residuos
sigma6 = y_identificacao(7:N)-psi6*theta6;
% criterio de informacao de Akaike (ordem 6)
AIC_6 = N*log(sigma6.^2)+2*length(theta6);

% matriz de regressores ordem 7
psi7 = [y_identificacao(toeplitz(7:N-1,[7 6 5 4 3 2 1]))...
    x_identificacao(toeplitz(7:N-1,[7 6 5 4 3 2 1]))];
% parametros estimados a1, a2, a3, a4, a5, a6,a7, b1, b2, b3, b4, b5 b6 e
% b7
theta7 = pinv(psi7'*psi7)*psi7'*y_identificacao(8:N);  
% vetor de residuos
sigma7 = y_identificacao(8:N)-psi7*theta7;
% criterio de informacao de Akaike (ordem 7)
AIC_7 = N*log(sigma7.^2)+2*length(theta7);

% matriz de regressores ordem 8
psi8 = [y_identificacao(toeplitz(8:N-1,[8 7 6 5 4 3 2 1]))...
    x_identificacao(toeplitz(8:N-1,[8 7 6 5 4 3 2 1]))];
% parametros estimados a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5,
%b6, b7 e b8
theta8 = pinv(psi8'*psi8)*psi8'*y_identificacao(9:N);  
% vetor de residuos
sigma8 = y_identificacao(9:N)-psi8*theta8;
% criterio de informacao de Akaike (ordem 8)
AIC_8 = N*log(sigma8.^2)+2*length(theta8);

% matriz de regressores ordem 9
psi9 = [y_identificacao(toeplitz(9:N-1,[9 8 7 6 5 4 3 2 1]))...
    x_identificacao(toeplitz(9:N-1,[9 8 7 6 5 4 3 2 1]))];
% parametros estimados a1, a2, a3, a4, a5, a6, a7, a8, a9, b1, b2, b3, b4,
%b5, b6, b7, b8 e b9
theta9 = pinv(psi9'*psi9)*psi9'*y_identificacao(10:N);  
% vetor de residuos
sigma9 = y_identificacao(10:N)-psi9*theta9;
% criterio de informacao de Akaike (ordem 9)
AIC_9 = N*log(sigma9.^2)+2*length(theta9);


% matriz de regressores ordem 10
psi10 = [y_identificacao(toeplitz(10:N-1,[10 9 8 7 6 5 4 3 2 1]))...
    x_identificacao(toeplitz(10:N-1,[10 9 8 7 6 5 4 3 2 1]))];
% parametros estimados a1, a2, a3, a4, a5, a6, a7, a8, a9, b1, b2, b3, b4,
%b5, b6, b7, b8 e b9
theta10 = pinv(psi10'*psi10)*psi10'*y_identificacao(11:N);  
% vetor de residuos
sigma10 = y_identificacao(11:N)-psi10*theta10;
% criterio de informacao de Akaike (ordem 9)
AIC_10 = N*log(sigma10.^2)+2*length(theta10);

AIC = [ min(AIC_1) min(AIC_2) min(AIC_3) min(AIC_4) min(AIC_5) min(AIC_6)...
    min(AIC_7) min(AIC_8) min(AIC_9) min(AIC_10)]
figure
plot(AIC,'b-','linewidth',2);
xlabel('Ordem do Modelo')
ylabel('AIC_{min}')
title('AIC_{min} x ordem do modelo utilizando N = 5000')
%% ITEM 3 
%coeficientes
theta = theta4;

%% ITEM 4 - VALIDAÇÃO DO MODELO
% a
r = [sigma4]; % residuos do modelo ARX de 4a ordem calculados no item 2

% calculo da FAC de r(k)
[acf_r,lags_r,bounds_r] = autocorr(r,200);

figure;
plot(lags_r,acf_r,'b-','linewidth',2);
ylim([-1 1]);
xlim([-1 200]);
hold on
plot(lags_r,bounds_r(1)*ones(length(lags_r)),'k--','linewidth',0.05)
plot(lags_r,bounds_r(2)*ones(length(lags_r)),'k--','linewidth',0.05)
xlabel('atraso')
ylabel('correlação')
title('Autocorrelação dos resíduos')

%b autocorrelaçao entre entrada e residuos
[FCC_r,lags_r] = xcov(x_identificacao,r);
FCC_media_r = mean(FCC_r); % media da FCC entre xi e r
FCC_sd_r = std(FCC_r); % desvio padrao da FCC entre xi e r

figure 
plot(lags_r,FCC_r,'b-','linewidth',2); 
hold on 
plot(lags_r,(FCC_media_r-2*FCC_sd_r)*ones(1,length(lags_r)),'k--');
plot(lags_r,(FCC_media_r+2*FCC_sd_r)*ones(1,length(lags_r)),'k--');
xlabel('Atrasos')
title('FCC entre os sinais x(k) e r(k)')
%% ITEM 5 - VALIDAÇÃO COM DADOS PARA VALIDAÇÃO

psi = [y_validacao(toeplitz(4:N-1,[4 3 2 1])) x_validacao(toeplitz(4:N-1,[4 3 2 1]))]; 

% saida da simulacao
yv_1(1,1) = y_validacao(1);
yv_1(2,1) = y_validacao(2);
yv_1(3,1) = y_validacao(3);
yv_1(4,1) = y_validacao(4);
yv_1 = [yv_1; psi*theta];

figure
plot(y_validacao,'b-','linewidth',2);
hold on;
plot(yv_1,'r--','linewidth',2);
xlabel('t(s)')
legend('y(t) original','y(t) simulação um passo a frente')
title('Comparação da resposta real do sistema com a resposta da simulação um passo a frente')

% b) simulacao livre 

yv_2(1) = y_validacao(1);
yv_2(2) = y_validacao(2);
yv_2(3) = y_validacao(3);
yv_2(4) = y_validacao(4);
for i = 5:N
    yv_2(i) = [yv_2(i-1) yv_2(i-2) yv_2(i-3) yv_2(i-4)...
        x_validacao(i-1) x_validacao(i-2) x_validacao(i-3) ...
        x_validacao(i-4)]*theta;
end
yv_2 = yv_2';

figure
plot(y_validacao,'b-','linewidth',2); 
hold on;
plot(yv_2,'r--','linewidth',2);
xlabel('t(s)')
legend('y(t) original','y(t) simulação livre')
title('Comparação da resposta real do sistema com a resposta da simulação livre')

%% ITEM 5(C) - Resíduos e Root mean square

r_yv1 = y_validacao - yv_1; % residuos da simulacao um passo a frente
r_yv2 = y_validacao - yv_2; % residuos da simulacao livre 


figure
plot(r_yv1,'b-','linewidth',2); 
xlabel('k')
title('Resíduos da simulação um passo a frente dos dados de validação (Passo a frente)')

figure
plot(r_yv2,'b-','linewidth',2);
xlabel('k')
title('Resíduos da simulação livre dos dados de validação (Simulação Livre)')

% RMSE simulacao um passo a frente 
RMSE_1 = sqrt(sum((y_validacao-yv_1).^2))/sqrt(sum((y_validacao-mean(y_validacao)).^2));

% RMSE simulacao livre 
RMSE_2 = sqrt(sum((y_validacao-yv_2).^2))/sqrt(sum((y_validacao-mean(y_validacao)).^2));