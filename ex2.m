%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                               EXERCÍCIO 2                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Trabalho de técnicas de modelagem     
%%%         Francis Santos               
%%%         2012022167                    
%       Mínimos Quadrados                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ITEM 1
baroreflex=load('baroreflex.mat'); %carrega matriz de entrada

dados_RR = baroreflex.RR(1:100);
dados_SAP= baroreflex.SAP(1:100);
%Calcula coeficientes N=1
% matriz de regressores (ordem 1)
psi1 = [dados_RR(toeplitz(1:99,1))' dados_SAP(toeplitz(2:100,[2 1]))]; 
% parametros estimados (a1, b0 e b1)
theta_2_1 = pinv(psi1'*psi1)*psi1'*dados_RR(2:100)'; 

%Matriz regressores ordem 2 
psi2 = [dados_RR(toeplitz(2:99,[2 1])) dados_SAP(toeplitz(3:100,[3 2 1]))];
%coeficientes estimados a1, a2, b0, b1 e b2
theta_2_2 = pinv(psi2'*psi2)*psi2'*dados_RR(3:100)'; 

%Matriz regressores ordem 3
psi3 = [dados_RR(toeplitz(3:99,[3 2 1])) dados_SAP(toeplitz(4:100,[4 3 2 1]))];
%coeficientes estimados a1, a2,a3, b0, b1 b2 e b3
theta_2_3 = pinv(psi3'*psi3)*psi3'*dados_RR(4:100)'; 

%Matriz regressores ordem 4
psi4 = [dados_RR(toeplitz(4:99,[4 3 2 1])) dados_SAP(toeplitz(5:100,[5 4 3 2 1]))];
%coeficientes estimados a1, a2,a3,a4, b0, b1 b2 b3 e b4
theta_2_4 = pinv(psi4'*psi4)*psi4'*dados_RR(5:100)'; 

%Matriz regressores ordem 5
psi5 = [dados_RR(toeplitz(5:99,[5 4 3 2 1])) ...
    dados_SAP(toeplitz(6:100,[6 5 4 3 2 1]))];
%coeficientes estimados a1, a2,a3,a4,a5, b0, b1 b2 b3 b4 e b5
theta_2_5 = pinv(psi5'*psi5)*psi5'*dados_RR(6:100)'; 

%Matriz regressores ordem 6
psi6 = [dados_RR(toeplitz(6:99,[6 5 4 3 2 1])) ...
    dados_SAP(toeplitz(7:100,[7 6 5 4 3 2 1]))];
%coeficientes estimados a1, a2,a3,a4,a5 a6, b0, b1 b2 b3 b4 b5 e b6
theta_2_6 = pinv(psi6'*psi6)*psi6'*dados_RR(7:100)'; 

%Matriz regressores ordem 7
psi7 = [dados_RR(toeplitz(7:99,[7 6 5 4 3 2 1])) ...
    dados_SAP(toeplitz(8:100,[8 7 6 5 4 3 2 1]))];
%coeficientes estimados a1, a2,a3,a4,a5 a6,a7, b0, b1 b2 b3 b4 b5 e b6 e b7
theta_2_7 = pinv(psi7'*psi7)*psi7'*dados_RR(8:100)'; 

%Matriz regressores ordem 8
psi8 = [dados_RR(toeplitz(8:99,[8 7 6 5 4 3 2 1])) ...
    dados_SAP(toeplitz(9:100,[9 8 7 6 5 4 3 2 1]))];
%coeficientes estimados a1, a2,a3,a4,a5 a6 e a7, b0, b1 b2 b3 b4 b5 b6 b7 e
%b8
theta_2_8 = pinv(psi8'*psi8)*psi8'*dados_RR(9:100)'; 

%% ITEM 2 - PREDIÇÃO COM UM PASSO A FRENTE

% predicao um passo a frente ordem 1
psi2_1 = [(baroreflex.RR(toeplitz(1:1199,1)))' baroreflex.SAP(toeplitz(2:1200,[2 1]))]; 
y_2_1 = [baroreflex.RR(1); psi2_1*theta_2_1];

% predicao um passo a frente ordem 2
psi2_2 = [(baroreflex.RR(toeplitz(2:1199,[2 1]))) ...
    baroreflex.SAP(toeplitz(3:1200,[3 2 1]))]; 
y_2_2 = [baroreflex.RR(1);baroreflex.RR(2); psi2_2*theta_2_2];

% predicao um passo a frente ordem 3
psi2_3 = [(baroreflex.RR(toeplitz(3:1199,[3 2 1]))) ...
    baroreflex.SAP(toeplitz(4:1200,[4 3 2 1]))]; 
y_2_3 = [baroreflex.RR(1);baroreflex.RR(2);baroreflex.RR(3); psi2_3*theta_2_3];


% predicao um passo a frente ordem 4
psi2_4 = [(baroreflex.RR(toeplitz(4:1199,[4 3 2 1]))) ...
    baroreflex.SAP(toeplitz(5:1200,[5 4 3 2 1]))]; 
y_2_4 = [baroreflex.RR(1);baroreflex.RR(2);baroreflex.RR(3);baroreflex.RR(4); psi2_4*theta_2_4];

% predicao um passo a frente ordem 5
psi2_5 = [(baroreflex.RR(toeplitz(5:1199,[5 4 3 2 1]))) ...
    baroreflex.SAP(toeplitz(6:1200,[6 5 4 3 2 1]))]; 
y_2_5 = [baroreflex.RR(1);baroreflex.RR(2);baroreflex.RR(3);...
    baroreflex.RR(4);baroreflex.RR(5); psi2_5*theta_2_5];

% predicao um passo a frente ordem 6
psi2_6 = [(baroreflex.RR(toeplitz(6:1199,[6 5 4 3 2 1]))) ...
    baroreflex.SAP(toeplitz(7:1200,[7 6 5 4 3 2 1]))]; 
y_2_6 = [baroreflex.RR(1);baroreflex.RR(2);baroreflex.RR(3);...
    baroreflex.RR(4);baroreflex.RR(5); baroreflex.RR(6); psi2_6*theta_2_6];

% predicao um passo a frente ordem 7
psi2_7 = [(baroreflex.RR(toeplitz(7:1199,[7 6 5 4 3 2 1]))) ...
    baroreflex.SAP(toeplitz(8:1200,[8 7 6 5 4 3 2 1]))]; 
y_2_7 = [baroreflex.RR(1);baroreflex.RR(2);baroreflex.RR(3);...
    baroreflex.RR(4);baroreflex.RR(5); baroreflex.RR(6);...
    baroreflex.RR(7); psi2_7*theta_2_7];

% predicao um passo a frente ordem 7
psi2_8 = [(baroreflex.RR(toeplitz(8:1199,[8 7 6 5 4 3 2 1]))) ...
    baroreflex.SAP(toeplitz(9:1200,[9 8 7 6 5 4 3 2 1]))]; 
y_2_8 = [baroreflex.RR(1);baroreflex.RR(2);baroreflex.RR(3);...
    baroreflex.RR(4);baroreflex.RR(5); baroreflex.RR(6);...
    baroreflex.RR(7);baroreflex.RR(7); psi2_8*theta_2_8];


%% ITEM 3 - SIMULAÇÃO LIVRE
%ordem 1
y_3_1(1) = baroreflex.RR(1);
for i = 2:1200
    y_3_1(i) = [y_3_1(i-1) baroreflex.SAP(i) baroreflex.SAP(i-1)]*theta_2_1;
end
y_3_1 = y_3_1'; 

%ordem 2
y_3_2(1) = baroreflex.RR(1);
y_3_2(2) = baroreflex.RR(2);
for i = 3:1200
    y_3_2(i) = [y_3_2(i-1) y_3_2(i-2) baroreflex.SAP(i) baroreflex.SAP(i-1) ...
        baroreflex.SAP(i-2)]*theta_2_2;
end
y_3_2 = y_3_2'; 

%ordem 3
y_3_3(1) = baroreflex.RR(1);
y_3_3(2) = baroreflex.RR(2);
y_3_3(3) = baroreflex.RR(3);
for i = 4:1200
    y_3_3(i) = [y_3_3(i-1) y_3_3(i-2) y_3_3(i-3) baroreflex.SAP(i) baroreflex.SAP(i-1) ...
        baroreflex.SAP(i-2) baroreflex.SAP(i-3)]*theta_2_3;
end
y_3_3 = y_3_3'; 

%ordem 4
y_3_4(1) = baroreflex.RR(1);
y_3_4(2) = baroreflex.RR(2);
y_3_4(3) = baroreflex.RR(3);
y_3_4(4) = baroreflex.RR(4);
for i = 5:1200
    y_3_4(i) = [y_3_4(i-1) y_3_4(i-2) y_3_4(i-3) y_3_4(i-4) ...
        baroreflex.SAP(i) baroreflex.SAP(i-1) ...
        baroreflex.SAP(i-2) baroreflex.SAP(i-3) baroreflex.SAP(i-4)]*theta_2_4;
end
 y_3_4 = y_3_4'; 

%ordem 5
y_3_5(1) = baroreflex.RR(1);
y_3_5(2) = baroreflex.RR(2);
y_3_5(3) = baroreflex.RR(3);
y_3_5(4) = baroreflex.RR(4);
y_3_5(5) = baroreflex.RR(5);
for i = 6:1200
    y_3_5(i) = [y_3_5(i-1) y_3_5(i-2) y_3_5(i-3) y_3_5(i-4) y_3_5(i-5) ...
        baroreflex.SAP(i) baroreflex.SAP(i-1) ...
        baroreflex.SAP(i-2) baroreflex.SAP(i-3) baroreflex.SAP(i-4) ...
        baroreflex.SAP(i-5)]*theta_2_5;
end
y_3_5 = y_3_5'; 

%ordem 6
y_3_6(1) = baroreflex.RR(1);
y_3_6(2) = baroreflex.RR(2);
y_3_6(3) = baroreflex.RR(3);
y_3_6(4) = baroreflex.RR(4);
y_3_6(5) = baroreflex.RR(5);
y_3_6(6) = baroreflex.RR(6);
for i = 7:1200
    y_3_6(i) = [y_3_6(i-1) y_3_6(i-2) y_3_6(i-3) y_3_6(i-4) y_3_6(i-5) ...
        y_3_6(i-6) baroreflex.SAP(i) baroreflex.SAP(i-1) ...
        baroreflex.SAP(i-2) baroreflex.SAP(i-3) baroreflex.SAP(i-4) ...
        baroreflex.SAP(i-5) baroreflex.SAP(i-6)]*theta_2_6;
end
 y_3_6 = y_3_6'; 

%ordem 7
y_3_7(1) = baroreflex.RR(1);
y_3_7(2) = baroreflex.RR(2);
y_3_7(3) = baroreflex.RR(3);
y_3_7(4) = baroreflex.RR(4);
y_3_7(5) = baroreflex.RR(5);
y_3_7(6) = baroreflex.RR(6);
y_3_7(7) = baroreflex.RR(7);
for i = 8:1200
    y_3_7(i) = [y_3_7(i-1) y_3_7(i-2) y_3_7(i-3) y_3_7(i-4) y_3_7(i-5) ...
        y_3_7(i-6) y_3_7(i-7) baroreflex.SAP(i) baroreflex.SAP(i-1) ...
        baroreflex.SAP(i-2) baroreflex.SAP(i-3) baroreflex.SAP(i-4) ...
        baroreflex.SAP(i-5) baroreflex.SAP(i-6)  baroreflex.SAP(i-7)]*theta_2_7;
end
y_3_7 = y_3_7'; 

%ordem 8
y_3_8(1) = baroreflex.RR(1);
y_3_8(2) = baroreflex.RR(2);
y_3_8(3) = baroreflex.RR(3);
y_3_8(4) = baroreflex.RR(4);
y_3_8(5) = baroreflex.RR(5);
y_3_8(6) = baroreflex.RR(6);
y_3_8(7) = baroreflex.RR(7);
y_3_8(8) = baroreflex.RR(8);
for i = 9:1200
    y_3_8(i) = [y_3_8(i-1) y_3_8(i-2) y_3_8(i-3) y_3_8(i-4) y_3_8(i-5) ...
        y_3_8(i-6) y_3_8(i-7) y_3_8(i-8) baroreflex.SAP(i) baroreflex.SAP(i-1) ...
        baroreflex.SAP(i-2) baroreflex.SAP(i-3) baroreflex.SAP(i-4) ...
        baroreflex.SAP(i-5) baroreflex.SAP(i-6)  baroreflex.SAP(i-7)...
        baroreflex.SAP(i-8)]*theta_2_8;
end
 y_3_8 = y_3_8'; 
%% ITEM 4 - PLOT

% plot ordem 1
figure,
plot(1:1200,baroreflex.RR,'k-','linewidth',2); box off;
hold on
plot(1:1200,y_2_1,'b-','linewidth',2); box off;
plot(1:1200,y_3_1,'r-','linewidth',2); box off;
xlabel('tempo(milissegundos)')
legend('RR(k) original','RR(k) item 2','RR(k) item 3')
title('Sistema de ordem 1')

% plot ordem 2
figure,
plot(1:1200,baroreflex.RR,'k-','linewidth',2); box off;
hold on
plot(1:1200,y_2_2,'b-','linewidth',2); box off;
plot(1:1200,y_3_2,'r-','linewidth',2); box off;
xlabel('tempo(milissegundos)')
legend('RR(k) original','RR(k) item 2','RR(k) item 3')
title('Sistema de ordem 2')

% plot ordem 3
figure,
plot(1:1200,baroreflex.RR,'k-','linewidth',2); box off;
hold on
plot(1:1200,y_2_3,'b-','linewidth',2); box off;
plot(1:1200,y_3_3,'r-','linewidth',2); box off;
xlabel('tempo(milissegundos)')
legend('RR(k) original','RR(k) item 2','RR(k) item 3')
title('Sistema de ordem 3')

% plot ordem 4
figure,
plot(1:1200,baroreflex.RR,'k-','linewidth',2); box off;
hold on
plot(1:1200,y_2_4,'b-','linewidth',2); box off;
plot(1:1200,y_3_4,'r-','linewidth',2); box off;
xlabel('tempo(milissegundos)')
legend('RR(k) original','RR(k) item 2','RR(k)item 3')
title('Sistema de ordem 4')

% plot ordem 5
figure,
plot(1:1200,baroreflex.RR,'k-','linewidth',2); box off;
hold on
plot(1:1200,y_2_5,'b-','linewidth',2); box off;
plot(1:1200,y_3_5,'r-','linewidth',2); box off;
xlabel('tempo(milissegundos)')
legend('RR(k) original','RR(k) item 2','RR(k) item 3')
title('Sistema de ordem 5')

% plot ordem 6
figure,
plot(1:1200,baroreflex.RR,'k-','linewidth',2); box off;
hold on
plot(1:1200,y_2_6,'b-','linewidth',2); box off;
plot(1:1200,y_3_6,'r-','linewidth',2); box off;
xlabel('tempo(milissegundos)')
legend('RR(k) original','RR(k) item 2','RR(k) item 3')
title('Sistema de ordem 6')

% plot ordem 7
figure,
plot(1:1200,baroreflex.RR,'k-','linewidth',2); box off;
hold on
plot(1:1200,y_2_7,'b-','linewidth',2); box off;
plot(1:1200,y_3_7,'r-','linewidth',2); box off;
xlabel('tempo(milissegundos)')
legend('RR(k) original','RR(k) item 2','RR(k) item 3')
title('Sistema de ordem 7')

% plot ordem 8
figure,
plot(1:1200,baroreflex.RR,'k-','linewidth',2); box off;
hold on
plot(1:1200,y_2_8,'b-','linewidth',2); box off;
plot(1:1200,y_3_8,'r-','linewidth',2); box off;
xlabel('tempo(milissegundos)')
legend('RR(k) original','RR(k) item 2','RR(k) item 3')
title('Sistema de ordem 8')
%% ITEM 5 - ERRO QUADRATICO
% erro quadratico (ordem 1, considerando a predicao de passo a frente)
 MSE_1 = mean((y_2_1 - baroreflex.RR').^2);
 MSE_2 = mean((y_2_2 - baroreflex.RR').^2);
 MSE_3 = mean((y_2_3 - baroreflex.RR').^2);
 MSE_4 = mean((y_2_4 - baroreflex.RR').^2);
 MSE_5 = mean((y_2_5 - baroreflex.RR').^2);
 MSE_6 = mean((y_2_6 - baroreflex.RR').^2);
 MSE_7 = mean((y_2_7 - baroreflex.RR').^2);
 MSE_8 = mean((y_2_8 - baroreflex.RR').^2);
 
 % erro quadratico (ordem 1, considerando a predicao livre)
%  MSE_11 = mean((y_3_1 - baroreflex.RR).^2);
%  MSE_22 = mean((y_3_2 - baroreflex.RR).^2);
%  MSE_33 = mean((y_3_3 - baroreflex.RR).^2);
%  MSE_44 = mean((y_3_4 - baroreflex.RR).^2);
%  MSE_55 = mean((y_3_5 - baroreflex.RR).^2);
%  MSE_66 = mean((y_3_6 - baroreflex.RR').^2);
%  MSE_77 = mean((y_3_7 - baroreflex.RR).^2);
%  MSE_88 = mean((y_3_8 - baroreflex.RR').^2);
 
fprintf('\n');

disp('          N    MSE item 2        MSE item 3   ')
disp(sprintf('  1     %f          %.4f          ', MSE_1,MSE_11));
disp(sprintf('  2     %f          %.4f         ', MSE_2,MSE_22));
disp(sprintf('  3     %f          %.4f         ', MSE_3,MSE_33));
disp(sprintf('  4     %f          %.4f         ', MSE_4,MSE_44));
disp(sprintf('  5     %f          %.4f         ', MSE_5,MSE_55));
disp(sprintf('  6     %f          %.4f         ', MSE_6,MSE_66));
disp(sprintf('  7     %f          %.4f         ', MSE_7,MSE_77));
disp(sprintf('  8     %f          %.4f         ', MSE_8,MSE_88));
