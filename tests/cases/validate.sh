#!/bin/bash

# La idea es:
# - Correr los dos programas con todos los test cases (1 al 6) utilizando la misma seed
# - Redirigir la salida de los programas a dos archivos .txt
# - Comparar esas salidas
# Se usa el path absoluto de los programas /home/.../programa 

# Validamos parámetros

if [ "$#" -ne 2 ]; then
	echo "Error: 2 arguments expected"
	echo "Usage: ./validate.sh <program1-path> <program2-path>"
	exit 1
fi

PROGRAM1="$1"
PROGRAM2="$2"

for P in "$PROGRAM1" "$PROGRAM2"; do
	if [ ! -x "$P" ]; then
		echo "Error: '$P' doesn't exists or isn't executable"
		exit 1
	fi
done

echo "######################################"
echo "Validando programas "
echo "######################################"

OUT1=$(mktemp)
OUT2=$(mktemp)

# Corremos todos los testcases, redirigimos la salida a un directorio temporal y luego comparamos sin las dos primeras lineas

for N in {1..6}; do
	PARAMS_FILE="params.case.0${N}"
	echo "-------------Caso 0${N}-----------------"

	source "$PARAMS_FILE"

	$PROGRAM1 $cmd > $OUT1

	$PROGRAM2 $cmd > $OUT2
	
	if diff <(tail -n +2 "$OUT1") <(tail -n +2 "$OUT2") > /dev/null; then
		echo "Case 0${N}:Programs are equivalent"
	else
		echo "Case 0${N}:Programs are different"
		diff <(tail -n +2 "$OUT1") <(tail -n +2 "$OUT2")
	fi
done

rm -f "$OUT1" "$OUT2"

echo "######################################"

exit 0
