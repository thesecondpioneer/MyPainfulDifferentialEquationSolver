clearvars
clear all
% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

figure(1);

[X, Y, Z] = csvimport('output.csv', 'columns', [1, 2, 3], 'noHeader', true);
X = log10(X);
Y = log10(Y);

plot(X, Y);

hold on;

x = -3.5:0.5:-2;
y = 2 * x;

plot(x,y,'g');

title('Норма точной полной погрешности метода из варианта','interpreter','latex', 'FontSize', 20);
xlabel('$log_{10}(h)$','interpreter','latex', 'FontSize', 20);
ylabel('$log_{10}(||y - y_{n}||)$','interpreter','latex', 'FontSize', 20);
legend('Погрешность', 'Линия наклона 2');

figure(2);

Z = log10(Z);
plot(X, Z);

hold on;

x = -3.5:0.5:-2;
y = 3 * x;

plot(x,y,'g');

title('Норма точной полной погрешности метода-оппонента','interpreter','latex', 'FontSize', 20);
xlabel('$log_{10}(h)$','interpreter','latex', 'FontSize', 20);
ylabel('$log_{10}(||y - y_{n}||)$','interpreter','latex', 'FontSize', 20);
legend('Погрешность', 'Линия наклона 3');

figure(3);

[X, Y] = csvimport('output2.csv', 'columns', [1, 2], 'noHeader', true);
Y = log10(Y);

plot(X, Y);
title('Зависимость нормы точной полной погрешности на оптимальном шаге от независимой переменной, метод из варианта','interpreter','latex', 'FontSize', 20);
xlabel('$x$','interpreter','latex', 'FontSize', 20);
ylabel('$log_{10}(||y - y_{i}||)$','interpreter','latex', 'FontSize', 20);

figure(4);

[X, Y] = csvimport('output3.csv', 'columns', [1, 2], 'noHeader', true);
Y = log10(Y);

plot(X, Y);
title('Зависимость нормы точной полной погрешности на оптимальном шаге от независимой переменной, метод-оппонент','interpreter','latex', 'FontSize', 20);
xlabel('$x$','interpreter','latex', 'FontSize', 20);
ylabel('$log_{10}(||y - y_{i}||)$','interpreter','latex', 'FontSize', 20);

figure(5)

[X, Y1, Y2, Y3, Y4] = csvimport('output4.csv', 'columns', [1, 2, 3, 4, 5], 'noHeader', true);
plot(X, Y1);
hold on
plot(X, Y2);
hold on
plot(X, Y3);
hold on
plot(X, Y4);
title('График решения, метод из варианта, автоматический выбор шага','interpreter','latex', 'FontSize', 20);
xlabel('$x$','interpreter','latex', 'FontSize', 20);
ylabel('$y_{i}$','interpreter','latex', 'FontSize', 20);
legend('y_1', 'y_2', 'y_3', 'y_4')

figure(6)

[X, Y1, Y2, Y3, Y4] = csvimport('output8.csv', 'columns', [1, 2, 3, 4, 5], 'noHeader', true);
plot(X, Y1);
hold on
plot(X, Y2);
hold on
plot(X, Y3);
hold on
plot(X, Y4);
title('График решения, метод-оппонент, автоматический выбор шага','interpreter','latex', 'FontSize', 20);
xlabel('$x$','interpreter','latex', 'FontSize', 20);
ylabel('$y_{i}$','interpreter','latex', 'FontSize', 20);
legend('y_1', 'y_2', 'y_3', 'y_4')

figure(7)
[X, Y] = csvimport('output5.csv', 'columns', [1, 2], 'noHeader', true);
plot(X, Y,'linestyle','none','marker','o', 'Color', 'g')
hold on
[X, Y] = csvimport('output6.csv', 'columns', [1, 2], 'noHeader', true);
plot(X, Y,'linestyle','none','marker','x', 'Color', 'r')
title('Размеры шагов в зависимости от независимой переменной, метод из варианта, автоматический выбор шага','interpreter','latex', 'FontSize', 20);
xlabel('$x$','interpreter','latex', 'FontSize', 20);
ylabel('$h$','interpreter','latex', 'FontSize', 20);
legend('Принятый шаг', 'Непринятый шаг')

figure(8)
[X, Y] = csvimport('output9.csv', 'columns', [1, 2], 'noHeader', true);
plot(X, Y,'linestyle','none','marker','o', 'Color', 'g')
hold on
[X, Y] = csvimport('output10.csv', 'columns', [1, 2], 'noHeader', true);
plot(X, Y,'linestyle','none','marker','x', 'Color', 'r')
title('Размеры шагов в зависимости от независимой переменной, метод-оппонент, автоматический выбор шага','interpreter','latex', 'FontSize', 20);
xlabel('$x$','interpreter','latex', 'FontSize', 20);
ylabel('$h$','interpreter','latex', 'FontSize', 20);
legend('Принятый шаг', 'Непринятый шаг')

figure(9)
[X, Y] = csvimport('output7.csv', 'columns', [1, 2], 'noHeader', true);
plot(X, Y);
title('Зависимость нормы точной полной погрешности на оптимальном шаге от независимой переменной, метод из варианта, автоматический выбор шага','interpreter','latex', 'FontSize', 15);
xlabel('$x$','interpreter','latex', 'FontSize', 20);
ylabel('$log_{10}(||y - y_{i}||)$','interpreter','latex', 'FontSize', 20);

figure(10)
[X, Y] = csvimport('output11.csv', 'columns', [1, 2], 'noHeader', true);
plot(X, Y);
title('Зависимость нормы точной полной погрешности на оптимальном шаге от независимой переменной, метод-оппонент, автоматический выбор шага','interpreter','latex', 'FontSize', 15);
xlabel('$x$','interpreter','latex', 'FontSize', 20);
ylabel('$log_{10}(||y - y_{i}||)$','interpreter','latex', 'FontSize', 20);


figure(11)
[X, Y] = csvimport('output12.csv', 'columns', [1, 2], 'noHeader', true);
X = log10(X);
Y = log10(Y);
plot(X, Y);
title('Зависимость количества вычислений правой части от относительной точности, метод из варианта','interpreter','latex', 'FontSize', 20);
xlabel('$log_{10}(rtol)$','interpreter','latex', 'FontSize', 20);
ylabel('$log_{10}(N)$','interpreter','latex', 'FontSize', 20);

figure(12)
[X, Y] = csvimport('output13.csv', 'columns', [1, 2], 'noHeader', true);
X = log10(X);
Y = log10(Y);
plot(X, Y);
title('Зависимость количества вычислений правой части от относительной точности, метод-оппонент','interpreter','latex', 'FontSize', 20);
xlabel('$log_{10}(rtol)$','interpreter','latex', 'FontSize', 20);
ylabel('$log_{10}(N)$','interpreter','latex', 'FontSize', 20);

