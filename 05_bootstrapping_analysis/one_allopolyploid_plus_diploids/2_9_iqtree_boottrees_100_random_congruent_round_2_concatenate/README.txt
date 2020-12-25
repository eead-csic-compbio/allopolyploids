See README.md

We move the new boottrees files to the new 2_9_iqtree_boottrees_100_random_congruent_round_2_concatenate directory
```
for d in *.relabelled;
do for f in *.boottrees;
do (cd $d && mv $f ../../2_9_iqtree_boottrees_100_random_congruent_round_2_concatenate);
done; done
```

We recover the 79 .boottrees files and rename the boottrees files
```
for file in *.boottrees;
do mv "$file" "$(basename "$file" .boottrees.relabelled.boottrees).boottrees";
done
```

