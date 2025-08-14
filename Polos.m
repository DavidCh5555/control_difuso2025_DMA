% Matriz de estados A
A = [0 1 0;
     1 1 -1;
     0 -1 0];
% Matriz de entrada B
B = [1;
     1;
     0];
% Matriz de salida C 
C = [1 0 1];
% Matriz de transmisión directa D
D = 0;

sys = ss(A, B, C, D);

fprintf('Valores propios del sistema sin control (A):\n');
disp(eig(A));

co = ctrb(A, B);
fprintf('Rango de la matriz de controlabilidad:\n');
disp(rank(co));

p = [-1 -3 -2];

K = place(A, B, p);
fprintf('Ganancia de retroalimentación K:\n');
disp(K);

Ac = A - B * K;

fprintf('Valores propios del sistema con control (Ac):\n');
disp(eig(Ac));

sys_cl = ss(Ac, B, C, D);

t = 0:0.01:8;
x0 = [1; 0; 0];

[y, t, x] = initial(sys_cl, x0, t);

figure;
subplot(3, 1, 1);
plot(t, x(:, 1));
title('Estado 1');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(3, 1, 2);
plot(t, x(:, 2));
title('Estado 2');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

subplot(3, 1, 3);
plot(t, x(:, 3));
title('Estado 3');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;

sgtitle('Estados en lazo cerrado');