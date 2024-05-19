#!/bin/bash
VALGRIND="valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --log-file=valgrind.out.%p"
mpirun -np $1 $VALGRIND $2 ${@:3}
