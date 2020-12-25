See README.md

We put together in a single boottrees file each of the hundred trees of each allele (it will be create into each directory the boottrees file):
```
for d in *.relabelled;
do for f in *.ph;
do (cd $d && cat $f | sed 's/;/;\n/g' > $d.boottrees);
done; done
```

We move the new boottrees files to a new directory 2_5_iqtree_boottrees_100_random_congruent_concatenate (run in the folder 2_4_iqtree_boottrees_100_random_congruent)
```
for d in *.relabelled;
do for f in *.boottrees;
do (cd $d && mv $f ../../2_5_iqtree_boottrees_100_random_congruent_concatenate);
done; done
```

We can rename boottrees files in 2_5_iqtree_boottrees_100_random_congruent_concatenate directory:

