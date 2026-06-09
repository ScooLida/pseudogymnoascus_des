#!/bin/bash

# ==========================================
# ВХОДНЫЕ ДАННЫЕ И НАСТРОЙКИ
# ==========================================
OUT_DIR="./genes_report"
INPUT_BREADTH="$OUT_DIR/genes_breadth.tsv"

# (Порог MIN_DEPTH=3 уже заложен внутри таблицы genes_breadth.tsv)
MIN_PERCENT=60       # Минимум 60% длины гена
MIN_SUCCESS_GENES=2  # 2 и более гена

GENES=("Pd_18S" "Pd_ITS" "Pd_28S" "Pd_MCM7" "Pd_TEF1alpha" "Pd_RPB2")
# ==========================================

echo "=== Шаг 1: Проверка файлов ==="
if [ ! -f "$INPUT_BREADTH" ]; then
    echo "Ошибка: Файл $INPUT_BREADTH не найден! Запустите скрипт из правильной папки."
    exit 1
fi

# Подготавливаем чистые временные файлы для всех 6 генов
for gene in "${GENES[@]}"; do
    echo -e "Sample\tPosition\tDepth" > "$OUT_DIR/${gene}_temp.tsv"
done


echo "=== Шаг 2: Формирование списка target_samples ==="
# awk собирает все подходящие имена образцов в одну строку, разделенную пробелами
target_samples=$(awk -v thresh="$MIN_PERCENT" -v min_genes="$MIN_SUCCESS_GENES" '
BEGIN { FS="\t" }
NR > 1 {
    sample = $1; val = $3;
    if (val + 0 >= thresh + 0) { success_count[sample]++ }
}
END {
    result = ""
    for (s in success_count) {
        if (success_count[s] >= min_genes) { result = result s " " }
    }
    print result
}' "$INPUT_BREADTH")

# Превращаем текстовую строку target_samples в полноценный массив elite_samples
read -r -a elite_samples <<< "$target_samples"
elite_count=${#elite_samples[@]}

if [ "$elite_count" -eq 0 ]; then
    echo "Образцов по критерию (>= 2 генов на >= 60%) не найдено. Остановка."
    exit 0
fi

echo "Список сформирован! Успешно отобрано образцов: $elite_count"
echo "В список вошли: $target_samples"


echo "=== Шаг 3: Выгрузка таблиц для 6 генов ==="
counter=0
for target in "${elite_samples[@]}"; do
    ((counter++))
    echo "[$counter/$elite_count] Выгрузка координат для образца: $target"
    
    # Ищем BAM-файл для текущего образца
    target_bam=$(find ./My_grib_genes/my_genome -name "${target}*rescaled.bam" | head -n 1)
    
    if [ -z "$target_bam" ]; then
        echo "  [!] BAM-файл для $target не найден, пропускаем."
        continue
    fi
    
    # Записываем покрытие для ВСЕХ 6 генов, включая нулевые участки
    for gene in "${GENES[@]}"; do
        samtools depth -a -r "$gene" "$target_bam" | awk -v sample="$target" '{print sample"\t"$2"\t"$3}' >> "$OUT_DIR/${gene}_temp.tsv"
    done
done


echo "=== Шаг 4: Переименование файлов ==="
# Прикрепляем к имени файлов финальное число отобранных образцов (например, 66)
for gene in "${GENES[@]}"; do
    mv "$OUT_DIR/${gene}_temp.tsv" "$OUT_DIR/${gene}_fil_cov_${elite_count}.tsv"
done

echo "========================================================"
echo "РАБОТА ЗАВЕРШЕНА!"
echo "Файлы успешно сохранены в папку $OUT_DIR:"
for gene in "${GENES[@]}"; do
    echo "  -> ${gene}_fil_cov_${elite_count}.tsv"
done
echo "========================================================"
