import numpy as np
import matplotlib.pyplot as plt

# Загрузка данных из файла
data = np.loadtxt("brusselator_output.txt", skiprows=1)
t = data[:, 0]
x = data[:, 1]
y = data[:, 2]

# Построение графиков
plt.figure(figsize=(10, 6))
plt.plot(t, x, label='x(t)')
plt.plot(t, y, label='y(t)')
plt.xlabel('Время')
plt.ylabel('Концентрация')
plt.title('Моделирование Брюсселятора (RK4)')
plt.legend()
plt.grid()
plt.show()

# Фазовый портрет
plt.figure(figsize=(6, 6))
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Фазовый портрет Брюсселятора')
plt.grid()
plt.show()