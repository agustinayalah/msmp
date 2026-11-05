#!/bin/bash

# Directorio del script
cd "$(dirname "$0")" || exit 1

# Parámetros del benchmark
VERSION="${1:-v1.0}"
CSV_FILE="resultados_${VERSION}.csv"

# Si se pasa la ruta a otro programa, usamos ese, caso contrario, usamos el default
if [ -z "$2" ]; then 
	EXEC="../../ms"
else
	EXEC="$2"
fi

# Si el csv no existe, lo creamos
if [ ! -f "$CSV_FILE" ]; then
    echo "VERSION_DE_PRUEBA,TESTCASE,ITERATION,PARAMS,TIEMPO" > "$CSV_FILE"
fi

echo "###########################################"
echo "Iniciando pruebas | Versión: $VERSION"
echo "Guardando en: $CSV_FILE"
echo "###########################################"

# Iterar sobre los 6 test cases
for case_num in {1..6}; do
    echo ""
    echo "--- Ejecutando Caso 0$case_num de 6 ---"
    ./testcase.sh "0$case_num" "$VERSION" "$CSV_FILE" "$EXEC"
done

echo ""
echo "###########################################"
echo "Pruebas completadas. Archivo generado: $CSV_FILE"
echo "###########################################"
