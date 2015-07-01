## simLoci
The program `simLoci` can be used to simulate sequence data on an input topology with gene flow. It uses a coalescent simulator from the Python package `Egglib` and also requires `Numpy`. It allows fast simulation of data under a range of parameter values and with data output in a variety of formats:  

     .LOCI (used by pyRAD)  
     .PHY  (used by raxml)  
     .NEX  (used by BEAST)  
     .STR  (used by STRUCTURE)  
     .MIG  (used by migrate-N)    

### Example usage
```
simLoci -h

Here is how you can use this script

Usage: python /home/deren/local/bin/simLoci
        -o <str>   output file prefix name (required)
        -L <int>   Number of loci to simulate (default 100) 
        -l <int>   length of loci (default 100) 
        -u <float> per-site mutation rate (default: 1e-9) 
        -N <int>   effective population size (default: 1e5) 
        -i <int>   number of sampled inds per tip taxon (default 1) 
        -t <str>   ultrametric tree w br lens (file or Newick string)
        -s1 <int>  random seed 1 (default 123456)
        -s2 <int>  random seed 2 (default 987654)
        -p <str>   print verbose : prints all params to screen
        -f <str>   output file format {l=loci (default), 
                                       p=phy, n=nex,
                                       a=alleles, v=vcf,
                                       t=treemix, k=structure, 
                                       g=geno, m=migrate} 
        -m <str>   migration events in string format:
                     e.g., 5e-5 migs/gen from b->c 
                     (forward in time) from time=0.5 to 0.25:
                     -m [c,b,0.5,0.25,5e-5]


simLoci -o "test1" -L 500 -l 50 -u 1e-8 -f "lap" 

simLoci -o "test2" -L 5000 -l 100 -u 1e-8 -N 2e5 -i 5

simLoci -o "test3" -L 500 -l 50 -u 1e-8 -f "lap" -t bigtree.tre

simLoci -o "test4" -L 500 -l 50 -u 1e-8 -f "lap" -t "((a:1,b:1):1,c:2);"

```


