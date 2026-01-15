#проверка на наличие гриба
#!/bin/bash

for file in ~/Numenius/*_R1.fastq;
 do
  name=$(basename "${file%_R1.fastq}")  

  echo "Обработка файла: $file" | tee -a otchet.txt

   ~/kraken2-2.1.3/bin/kraken2 --db my_fungi_db \
     --threads 10 --paired \
     --report "${name}_rep.txt" \
     --output "${name}_out.txt" \
   ~/Numenius/"${name}_R1.fastq" ~/Numenius/"${name}_R2.fastq" 
  done > otchet.txt
touch all_rep.txt
> all_rep.txt
 for  rep in *_rep.txt; do
  Sname=$(basename "${rep%_rep.txt}")
  if grep  -q "Pseudogymnoascus destructans" "$rep"; then echo "$Sname" >> all_rep.txt
 fi
done

