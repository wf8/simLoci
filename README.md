## simLoci (under development)
The program `simLoci` can be used to simulate sequence data on an input topology with gene flow. It uses a coalescent simulator from the Python package `Egglib` and also requires `Numpy`. It allows fast simulation of data under a range of parameter values and with data output in a variety of formats:  

     .LOCI (used by pyRAD)  
     .PHY  (used by raxml)  
     .NEX  (used by BEAST)  
     .MIG  (used by migrate-N)    
     .ALLELES (similar to .loci but with diploid alleles)

### Example usage
```
python simLoci.py -h

Here is how you can use this script

Usage: python simLoci.py
        -o <str>   output file prefix name (required)
        -L <int>   Number of loci to simulate (default 100) 
        -l <int>   length of loci (default 100) 
        -u <float> per-site mutation rate (default: 1e-9) 
        -N <int>   effective population size (default: 1e5) 
        -i <int>   number of sampled inds per tip taxon (default 1) 
        -t <str>   ultrametric tree w br lens (file or Newick string)
        -s1 <int>  random seed 1 (default 123456)
        -s2 <int>  random seed 2 (default 987654)
        -v         print version
        -f <str>   output file format {l=loci (default), 
                                       p=phy, n=nex,
                                       a=alleles, v=vcf,
                                       t=treemix, k=structure, 
                                       g=geno, m=migrate} 
        -m <str>   migration events in string format:
                     e.g., 5e-5 migs/gen from b->c 
                     (forward in time) from time=0.5 to 0.25:
                     -m [c,b,0.5,0.25,5e-5]


python simLoci.py -o test1 -L 500 -l 50 -u 1e-8 -f lap 

python simLoci.py -o test2 -L 5000 -l 100 -u 1e-8 -N 2e5 -i 5

python simLoci.py -o test3 -L 500 -l 50 -u 1e-8 -f pn -t bigtree.tre

python simLoci.py -o test4 -L 500 -l 50 -u 1e-8 -f lmap -t "((a:1,b:1):1,c:2);"

python simLoci.py -o test5 -m [b,c,0.5,0.25,5e-5]/[b,a,0.1,0.05,5e-5]

```


