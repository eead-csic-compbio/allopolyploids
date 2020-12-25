See README.md

## 2.2) Run 1000 non-parametric bootstrapping replicates (iqtree)

```
ls *.fna | parallel --gnu -j 60 iqtree-1.6.9-Linux/bin/iqtree -b 1000 -nt 1 -AICc -s {} :::
```
You can run it using a bash script and nohup:
```
nohup iqtree_non_parametric_bootstrapping_1000.sh > nohup_bootstrapping_1000.txt &
```

These .boottrees files will be saved in the folder 2_2_iqtree_boottrees_1000.

