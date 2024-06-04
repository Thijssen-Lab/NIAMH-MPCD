This example shows how the checkpointer works. The simulation is of a diffusing colloid. 
    - To explore how the checkpointer works, give `inputPrepareSim.json' as the input file and make an output directory 'originalOut', i.e. $mkdir originalOut; mpcd.out -i inputPrepareSim.json -o originalOut 
        - Before the simulation is complete, cancel the simulation. 
    - To restart the simulation from the last checkpoint, run give `inputRestoreFromCheckpoint.json' as the input file and save the output in a new directory, i.e. $mkdir savedOut; mpcd.out -i inputRestoreFromCheckpoint.json -o savedOut 
    - Notice that the results in savedOut start from the last checkpoint. 
        - **THIS WILL OVERWRITE YOUR ORIGINAL DATA IF YOU GIVE IT THE SAME OUTPUT PATH**
        - You must manually merge the results