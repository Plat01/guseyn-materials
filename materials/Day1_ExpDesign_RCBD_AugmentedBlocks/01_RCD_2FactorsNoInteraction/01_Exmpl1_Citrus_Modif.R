###############################################################################################
###         First session of International Theoretical & Practical Hands-on Course on       ###
###          "Using whole-genome data in plant breeding: GWAS and Genomic Selection         ###
###           to genetically characterize key traits in large germplasm collections"        ###
###                         23-27 July, 2018, VIR, St Petersburg, Russia                    ###
###                                       ----------                                        ###
###                         Day 1: Complete Randomized Factorial Design                     ###
###                                                                                         ###
###         Case Study #1: Variety and Shade effect on lemon tolerance to drought           ###
###                                          ----                                           ###
###           R Script edited by Pr Laurent Gentzbittel & Dr Cécile Ben, INP-ENSAT          ###
###############################################################################################

# Пример: влияние сорта Variety и затенения Shade на показатель Index.
# Модель содержит главные эффекты двух факторов без взаимодействия.
# Перед запуском в RStudio выберите папку этого скрипта через
# Session -> Set Working Directory -> Choose Directory.
# Дополнительные пакеты не нужны. Кодировка файла — UTF-8.
# cat() печатает название шага, print() — результат, str() — структуру.
# Явный print() обеспечивает вывод и при запуске через source():
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/source.html

# Закрываем графические устройства и показываем оставшиеся (NULL — нет устройств).
cat("\n=== Закрытие графиков ===\n")
graphics.off()
print(dev.list())

# Удаляем объекты рабочего окружения, как в исходном скрипте.
# character(0) после очистки означает, что список объектов пуст.
cat("\n=== Очистка рабочего окружения ===\n")
rm(list = ls())
print(ls())

# Задаём contr.sum для неупорядоченных и contr.poly для упорядоченных факторов.
cat("\n=== Настройка контрастов ===\n")
options(contrasts = c("contr.sum", "contr.poly"))
print(getOption("contrasts"))

# Читаем CSV: разделитель — точка с запятой, первая строка — имена столбцов.
# Явно преобразуем текстовые столбцы в факторы независимо от версии R.
# https://stat.ethz.ch/R-manual/R-devel/library/utils/html/read.table.html
cat("\n=== Чтение 01_Exmpl1_Citrus.csv: все наблюдения ===\n")
Lemons <- read.table("01_Exmpl1_Citrus.csv", sep = ";", header = TRUE,
                     stringsAsFactors = TRUE)
print(Lemons)

# Показываем типы столбцов и число наблюдений.
cat("\n=== Структура таблицы ===\n")
str(Lemons)

# Показываем исходный порядок уровней обоих факторов.
cat("\n=== Исходные уровни Variety ===\n")
print(levels(Lemons$Variety))
cat("\n=== Исходные уровни Shade ===\n")
print(levels(Lemons$Shade))

# Меняем порядок уровней Shade и выводим результат преобразования.
# При contr.sum эффект последнего уровня равен минус сумме остальных эффектов.
# Это отличается от опорного уровня с нулевым эффектом при contr.treatment.
cat("\n=== Порядок Shade: Sun, Shade, StrongShade ===\n")
Lemons$Shade <- factor(Lemons$Shade, levels = c("Sun", "Shade", "StrongShade"))
print(Lemons$Shade)
print(levels(Lemons$Shade))
print(contrasts(Lemons$Shade))  # Показываем кодирование уровней для модели.

# Оцениваем аддитивную модель с эффектами сорта и затенения.
cat("\n=== Модель: Index ~ Variety + Shade ===\n")
modele <- aov(Index ~ Variety + Shade, data = Lemons)
print(modele)

# Выводим таблицы оценённых эффектов факторов.
cat("\n=== Таблицы эффектов ===\n")
print(model.tables(modele, type = "effect"))

# Показываем матрицу модели: строки — наблюдения, столбцы — параметры.
cat("\n=== Матрица модели ===\n")
print(model.matrix(modele))

# Показываем таблицу дисперсионного анализа.
cat("\n=== Дисперсионный анализ ===\n")
print(summary(modele))

# Разбиваем вклад Shade по столбцам контрастов 1 и 2, как в исходном скрипте.
# SunvsStrongShade и ShadevsStrongShade — исходные подписи частей разбиения,
# а не результаты отдельных парных тестов средних.
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/summary.aov.html
cat("\n=== Разбиение вклада Shade по контрастам ===\n")
print(summary(modele, split = list(
  Shade = list(SunvsStrongShade = 1, ShadevsStrongShade = 2)
)))

# Сравниваем три вложенные модели. R0 содержит только общее среднее.
cat("\n=== Модель R0: только общее среднее ===\n")
R0 <- aov(Index ~ 1, data = Lemons)
print(R0)

# R1 добавляет эффект сорта; свободный член включён по умолчанию.
cat("\n=== Модель R1: общее среднее и Variety ===\n")
R1 <- aov(Index ~ Variety, data = Lemons)
print(R1)

# R2 добавляет затенение к модели с сортом; формула совпадает с modele.
cat("\n=== Модель R2: общее среднее, Variety и Shade ===\n")
R2 <- aov(Index ~ Variety + Shade, data = Lemons)
print(R2)

# RSS — сумма квадратов остатков e = Y - X * B.
# Из таблицы ANOVA берём строку Residuals и столбец Sum Sq по именам,
# вместо извлечения по числовым позициям из unlist(summary(...)).
cat("\n=== Остаточная сумма квадратов R0 ===\n")
RSS0 <- summary(R0)[[1]]["Residuals", "Sum Sq"]
print(RSS0)

# Аналогично получаем RSS модели с сортом.
cat("\n=== Остаточная сумма квадратов R1 ===\n")
RSS1 <- summary(R1)[[1]]["Residuals", "Sum Sq"]
print(RSS1)

# Получаем RSS модели с обоими факторами.
cat("\n=== Остаточная сумма квадратов R2 ===\n")
RSS2 <- summary(R2)[[1]]["Residuals", "Sum Sq"]
print(RSS2)

# Вычисляем уменьшение RSS при добавлении сорта.
cat("\n=== Уменьшение RSS при добавлении Variety: RSS0 - RSS1 ===\n")
print(RSS0 - RSS1)

# Вычисляем уменьшение RSS при добавлении затенения после сорта.
cat("\n=== Уменьшение RSS при добавлении Shade после Variety: RSS1 - RSS2 ===\n")
print(RSS1 - RSS2)
