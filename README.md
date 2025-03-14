# Брюсселятор
Брюсселятор — это теоретическая модель типа автокаталитической реакции. 
Он характеризуется реакциями

![image](https://github.com/user-attachments/assets/aaadd3df-dae2-455a-8cc2-2f4cfbd77d35)

В условиях, когда A и B находятся в большом избытке и, следовательно, могут быть смоделированы при постоянной концентрации, уравнения скорости становятся

![image](https://github.com/user-attachments/assets/7b4d9949-092c-4df6-ae26-9d09d73e8c5c)

где для удобства константы скорости были установлены равными 1.
Брюсселатор имеет фиксированную точку в

![image](https://github.com/user-attachments/assets/64acf41d-40da-499a-ae03-79d5a13251d7)

Фиксированная точка становится нестабильной, когда

![image](https://github.com/user-attachments/assets/16fd2351-5683-4238-a399-3e758a8b9da1)

что приводит к колебаниям в системе.

## Описание программы
Эта программа моделирует динамику системы Брюсселятора — химической реакции с автоколебаниями. Численное решение дифференциальных уравнений выполняется с помощью двух методов:
- **Метод Рунге-Кутты 4-го порядка (RK4)**
- **Метод Дорманда-Принса 5-го порядка (с адаптивным шагом)**

Программа сравнивает точность этих методов, вычисляя абсолютные и относительные ошибки.

## Файлы
- `main.cpp` — основной код программы.
- `brusselator_output.txt` — файл с результатами численного решения.
- `error1.txt` — файл с абсолютными ошибками.
- `error2.txt` — файл с относительными ошибками.
- `main.py` — скрипт на Python для построения графиков.

## Компиляция и запуск

### Компиляция
Компиляция выполняется с использованием g++:
```bash
 g++ main.cpp -o main.exe
```

### Запуск программы
```bash
 ./main
```

### Запуск скрипта для построения графиков
```bash
 python3 main.py
```

## Выходные данные
Программа создает и записывает результаты в файлы:
- `brusselator_output.txt` содержит временные ряды x(t) и y(t).
- `error1.txt` и `error2.txt` содержат ошибки между методами RK4 и Дорманда-Принса.

## Визуализация
После выполнения программы можно запустить `main.py`, который строит:
- Графики изменения концентраций x(t) и y(t).
- Фазовый портрет системы.
