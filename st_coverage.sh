
echo -e "Sample\tnumreads\tcoverage\tmeandepth" > merg_st.csv
# Получаем список уникальных имен (SM)
SAMPLES=$(samtools view -H My_grib3_18s.my_genome.bam | grep "@RG" | sed -e 's/.*SM:\([^ \t]*\).*/\1/' | sed 's/[\.-][0-9]$//' | sort -u)

for base_name in $SAMPLES; do
    echo "Processing Group: $base_name ..."
    
    # Ключевое изменение здесь: grep -E "^@|RG:Z:${base_name}"
    # Он берет ВСЕ строки заголовка (начинаются с @) И риды этого образца
    DATA=$(samtools view -h My_grib3_18s.my_genome.bam | grep -E "^@|RG:Z:${base_name}" | samtools coverage - | tail -n 1)
    
    NUMREADS=$(echo "$DATA" | cut -f 4)
    COV=$(echo "$DATA" | cut -f 6)
    DEPTH=$(echo "$DATA" | cut -f 7)
    
    # Проверяем, что это не пустая строка и не заголовок
    if [[ ! -z "$NUMREADS" && "$NUMREADS" != "numreads" && "$NUMREADS" != "0" ]]; then
        echo -e "$base_name\t$NUMREADS\t$COV\t$DEPTH" >> merg_st.csv
    fi
done
cat  merg_st.csv | awk '{print $1","$2","$3","$4}' > merg_st.txt
rm merg_st.csv 
echo "Готово! Результаты в файле merg_st.txt"















