#!/bin/bash

# ==========================================
# НАСТРОЙКИ РАСЧЕТА
# ==========================================
OUT_DIR="./genes_report"
mkdir -p "$OUT_DIR"

BREADTH_FILE="$OUT_DIR/genes_breadth.tsv"

# Главный критерий: позиция засчитывается, если на ней >= 3 чтений
MIN_DEPTH=3          

GENES=("Pd_18S" "Pd_ITS" "Pd_28S" "Pd_MCM7" "Pd_TEF1alpha" "Pd_RPB2")
# ==========================================

echo "=== Шаг 1: Подготовка таблицы ==="
# Создаем файл и записываем заголовки столбцов
echo -e "Sample\tGene\tBreadthPercentage" > "$BREADTH_FILE"


echo "=== Шаг 2: Поиск BAM-файлов ==="
total_files=$(ls ./My_grib_genes/my_genome/*/*rescaled.bam 2>/dev/null | wc -l)
if [ "$total_files" -eq 0 ]; then
    echo "Ошибка: BAM-файлы не найдены!"
    exit 1
fi
echo "Найдено образцов для анализа: $total_files"
echo "Считаем % позиций с глубиной >= $MIN_DEPTH"


echo "=== Шаг 3: Расчет покрытия через samtools ==="
counter=0
for bam in ./My_grib_genes/my_genome/*/*rescaled.bam; do
    [ -e "$bam" ] || continue
    ((counter++))
    
    # Извлекаем чистое имя образца
    sample_name=$(basename "$bam" | sed 's/\.rescaled\.bam//' | sed 's/\.my_genome\.bam//' | sed 's/\.bam//')
    
    echo "[$counter/$total_files] Обработка: $sample_name ..."
    
    # Индексируем BAM-файл (если он еще не проиндексирован)
    samtools index "$bam"
    
    # Запускаем цикл по всем 6 генам
    for gene in "${GENES[@]}"; do
        
        # samtools извлекает глубину, awk считает позиции >= MIN_DEPTH и выдает процент
        breadth=$(samtools depth -a -r "$gene" "$bam" | awk -v mindepth="$MIN_DEPTH" '
            {
                total++; 
                if ($3 >= mindepth) covered++
            } 
            END {
                if (total == 0) print 0;
                else print (covered / total) * 100
            }
        ')
        
        # Записываем результат в общую таблицу
        echo -e "$sample_name\t$gene\t$breadth" >> "$BREADTH_FILE"
        
    done
done

echo "========================================================"
echo "РАСЧЕТ ЗАВЕРШЕН!"
echo "Файл с процентами сохранен в: $BREADTH_FILE"
echo "Теперь вы можете запускать скрипт отбора (target samples)."
echo "========================================================"
