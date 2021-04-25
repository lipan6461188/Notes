

PCUT=(0.05 0.15)
modes=(0 1 2)

for mode in ${modes[@]}; do
    for pcut in ${PCUT[@]}; do
        python trRosetta_AlphaFold_constrait.py -m $mode -pd ${pcut} T0949.pickle T0949.torsions models/T0949_mode${mode}_pcut${pcut}_orient &
        python trRosetta_AlphaFold_constrait.py -m $mode -pd ${pcut} --no-orient T0949.pickle T0949.torsions models/T0949_mode${mode}_pcut${pcut}_noorient &
    done
    wait
done


