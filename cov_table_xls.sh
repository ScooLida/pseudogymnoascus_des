#!/bin/bash

# ==========================================
# ПАПКА С НАШИМИ РЕЗУЛЬТАТАМИ
# ==========================================
OUT_DIR="./genes_report"

echo "=== Шаг 1: Поиск готовых таблиц покрытия ==="
# Ищем все файлы, которые заканчиваются на _fil_cov_*.tsv (наши готовые результаты)
# Игнорируем файлы, которые уже имеют суффикс _wide.tsv (на случай повторного запуска)
files_to_process=$(ls "$OUT_DIR"/*_fil_cov_*.tsv 2>/dev/null | grep -v "_wide.tsv")

if [ -z "$files_to_process" ]; then
    echo "Ошибка: Файлы *_fil_cov_*.tsv не найдены в папке $OUT_DIR."
    exit 1
fi


echo "=== Шаг 2: Конвертация в широкие матрицы ==="

# Перебираем каждый найденный файл
for INPUT_FILE in $files_to_process; do
    
    # Формируем имя для нового файла (добавляем _wide перед расширением .tsv)
    OUTPUT_FILE="${INPUT_FILE%.tsv}_wide.tsv"
    
    echo "Транспонирую: $(basename "$INPUT_FILE") -> $(basename "$OUTPUT_FILE") ..."
    
    # Транспонируем данные через awk
    awk -F'\t' '
    BEGIN { OFS="\t" }
    NR > 1 {
        sample = $1; pos = $2; depth = $3;
        
        # Сохраняем уникальные позиции
        if (!(pos in pos_seen)) {
            positions[++p_count] = pos;
            pos_seen[pos] = 1;
        }
        
        # Сохраняем уникальные образцы
        if (!(sample in samp_seen)) {
            samples[++s_count] = sample;
            samp_seen[sample] = 1;
        }
        
        # Сохраняем глубину в двумерный массив
        data[sample, pos] = depth;
    }
    END {
        # Печатаем шапку: "Sample" + все номера позиций
        printf "Sample"
        for (i = 1; i <= p_count; i++) {
            printf "\t%s", positions[i]
        }
        printf "\n"
        
        # Печатаем строки: имя образца + все его значения глубин
        for (j = 1; j <= s_count; j++) {
            s = samples[j];
            printf "%s", s
            
            for (i = 1; i <= p_count; i++) {
                p = positions[i];
                # Если в ячейке пусто (не было чтений), принудительно ставим 0
                val = (data[s, p] != "") ? data[s, p] : 0;
                printf "\t%s", val
            }
            printf "\n"
        }
    }' "$INPUT_FILE" > "$OUTPUT_FILE"
    
done

echo "========================================================"
echo "ВСЕ ГОТОВО!"
echo "Широкие матрицы для всех генов успешно созданы в папке: $OUT_DIR"
echo "========================================================"
