#!/bin/bash

# Asegura que el script se ejecute desde el directorio donde reside
cd "$(dirname "$0")" || exit 1

if [ $# -eq 0 ]; then
    echo "Uso: $0 <testcase> [version_prueba] [ruta_csv]"
    exit 1
fi

testcase_raw=$1	# El numero sin el 0 antepuesto
testcase=$(printf "%02d" "$testcase_raw") # Formateo del parámetro

version_prueba=$2
ruta_csv=${3:-"resultados.csv"}

params="./params.case.$testcase"

if [ -f "$params" ]; then
    casedir="case$testcase"
    
    if [ ! -d "$casedir" ]; then
        mkdir -p "$casedir"
        echo "Created output folder $casedir"
        
        # Carga las variables del archivo de parámetros (como $cmd)
        source "$params"
        
        # Ruta al ejecutable (Ajustar si el ejecutable no se llama mspar o está en otro sitio)
        EXEC="../../ms"
        
        if [ ! -f "$EXEC" ]; then
            echo "Error: No se encontró el ejecutable en $EXEC"
            exit 1
        fi

        for i in {1..10}; do
            echo "-> Running iteration ${i} of 10 for test case $testcase"
            filename="${casedir}/case${testcase}.${i}"
            
            STARTTIME=$(date +%s%3N)
            
            /usr/bin/time -v -o "$filename.time" "$EXEC" ${cmd} > "$filename.out"
            
            ENDTIME=$(date +%s%3N)
            
            RESTIME=$(echo "scale=3 ; ($ENDTIME - $STARTTIME) / 1000" | bc -l) 
            
            echo "--> Done in ${RESTIME} seconds."
            
            CLEAN_PARAMS=$(echo "$cmd" | tr -s ' ' ',')
            
            echo "${version_prueba},case${testcase},${i},\"${CLEAN_PARAMS}\",${RESTIME}" >> "$ruta_csv"
            
        done
        
        echo "Test case ${testcase} done!"
    else
        echo "Output directory $casedir already exists."
    fi
else
    echo "File $params does not exist"
fi
