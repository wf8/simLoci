###############################################################
##  Simulate sequence data on an input tree w/ introgression ##
##  output in .loci readable format                          ##
##  author:  Deren Eaton                                     ##
##  contact: deren.eaton@yale.edu                            ##
##  date:    10/9/14                                         ##
##  version: 1.0                                             ##
###############################################################

## load standard modules
import getopt
import itertools
import sys

## load external modules
try: 
    import egglib
except ImportError:
    print "Python package Egglib not found"
try:
    import numpy as np
except ImportError:
    print "Python package Numpy not found"

##
def parseopts(opts):
    """ parses options to the command line """

    params = {'Loci':   100,
              'length': 100,
              'N':      1e6,
              'u':      1e-9,
              'inds':   1,
              'tree':   0,
              'migs':   "",
              'outname':"",
              'outform':"",
              'seed':   123456,
              'verbose':False
              }

    for opt,arg in opts:

        if opt in ["-L"]:
            params["Loci"]    = int(float(arg))
        elif opt in ["-l"]:
            params["length"]  = int(float(arg))
        elif opt in ["-u"]:
            params["mu"]      = float(arg)
        elif opt in ["-N"]:
            params["N"]       = int(float(arg))
        elif opt in ["-i"]:
            params["inds"]    = int(arg)
        elif opt in ["-t"]:
            params["tree"]    = str(arg)
        elif opt in ["-m"]:
            params["migs"]    = str(arg)
        elif opt in ["-o"]:
            params["outname"] = str(arg)
        elif opt in ["-f"]:
            params["outform"] = str(arg)
        elif opt in ["-v"]:
            params["verbose"] = True
        elif opt in ["-s"]:
            params["seed"] = int(arg)


    return params


def checkopts(opts):
    """ print errors if options are incorrect"""
    TODO

    
def usage():
    """
    brief description of various flags and options for this script
    """
    print "\nHere is how you can use this script\n"
    print "Usage: python %s"%sys.argv[0]
    print "\t -L <int>   Number of loci to simulate (default 100)"
    print "\t -l <int>   length of loci (default 100)"
    print "\t -u <float> per-site mutation rate (default: 1e-9"
    print "\t -N <int>   effective population size (default: 1e6)"
    print "\t -i <int>   number of sampled individuals per tip taxon (default 1)"
    print "\t -t <str>   tree as a newick string with branch lengths. No polytomy."
    print "\t -m <str>   migration events in string format (see documentation) "
    print "\t -o <str>   output file prefix "
    print "\t -f <str>   output file format {m=migrate,l=loci,p=phy,n=nex} "
    print "\t -s <int>   random seed "
    print "\t -v         verbose output - prints all params to screen"
    


def simdata(params):
    "split multiple migrations"
    migscenarios = [m.lstrip("[").rstrip("]") for m in params["migs"].split('/')]

    "if tree, parse it. Else use fixed tree"
    if params["tree"]: 
        tre = egglib.Tree(string=params["tree"])
        tiptax = tre.all_leaves()
    else:
        tiptax = list("ABCDEFGHI")    
        "clade 1 = ((A,B),C)"
        "clade 2 = ((D,E),F)"
        "clade 3 = ((G,H),I)"
        nodeABCDEFGHI = 3.0 
        nodeABCDEF = 1.5  
        nodeGHI   =  0.75 
        nodeDEF   =  0.75 
        nodeABC   =  0.75 
        nodeGH    =  0.5 
        nodeDE    =  0.5 
        nodeAB    =  0.5 


    "function to return node height"
    def lengthtotip(node):
        dec = 0.
        while node.descendants():
            dec += node.branch_to(node.descendants()[0])
            node = node.descendants()[0]
        return dec


    "sets the two parameter classes"
    paramSet = egglib.simul.CoalesceParamSet(singleSamples=None,
                                             doubleSamples=[params["inds"]]*len(tiptax),
                                             M=0.0)

    "traverse tree fusing populations at each node"
    if params["tree"]:
        for node in tre:
            if node.descendants():
                date = lengthtotip(node)
                s1 = min([tiptax.index(l) for l in node.descendants()[0].leaves_down()])
                s2 = min([tiptax.index(l) for l in node.descendants()[1].leaves_down()])
                paramSet.populationFusion(date,s1,s2)
    else:
        "clade 1"
        paramSet.populationFusion(nodeABC, 0,2)     
        paramSet.populationFusion(nodeAB, 0,1)      
        "clade 2"
        paramSet.populationFusion(nodeDEF, 3,5)     
        paramSet.populationFusion(nodeDE, 3,4)      
        "clade 3"
        paramSet.populationFusion(nodeGHI, 6,8)   
        paramSet.populationFusion(nodeGH, 6,7)    
        "together and outgroup"
        paramSet.populationFusion(nodeABCDEFGHI, 0,6)  
        paramSet.populationFusion(nodeABCDEF, 0,3)     


    "sets migration patterns "
    if params["migs"]:
        for mig in migscenarios:
            p2,p3,s,e,m = mig.split(',')
            p2 = tiptax.index(p2)
            p3 = tiptax.index(p3)
            M = 4.*params["N"]*float(m)
            paramSet.changePairwiseMigrationRate(0.0, int(p2), int(p3), 0.)
            paramSet.changePairwiseMigrationRate(float(s), int(p2), int(p3), float(M))
            paramSet.changePairwiseMigrationRate(float(e), int(p2), int(p3), 0.)

    mutator = egglib.simul.CoalesceFiniteAlleleMutator(theta=params["theta"],
                                                       alleles= 4,
                                                       randomAncestralState=True)
    mutator.setSites(params["length"])
    aligns = egglib.simul.coalesce(paramSet,
                                   mutator,
                                   params["Loci"],
                                   random=np.random.seed(params["seed"]))
    return aligns, tiptax


def hetero(n1,n2):
    """
    returns IUPAC symbol for ambiguity bases,
    used for polymorphic sites.
    """
    D = {('G','A'):"R",
         ('G','T'):"K",
         ('G','C'):"S",
         ('T','C'):"Y",
         ('T','A'):"W",
         ('C','A'):"M"}
    a = D.get((n1,n2))
    b = D.get((n2,n1))
    if a:
        return a
    else:
        return b


def unstruct(amb):
    amb = amb.upper()
    " returns bases from ambiguity code"
    D = {"R":["G","A"],
         "K":["G","T"],
         "S":["G","C"],
         "Y":["T","C"],
         "W":["T","A"],
         "M":["C","A"],
         "A":["A","A"],
         "T":["T","T"],
         "G":["G","G"],
         "C":["C","C"],
         "N":["N","N"],
         "-":["-","-"]}
    return D.get(amb)


def makehetero(seq1,seq2):
    seq = ""
    for base in zip(seq1,seq2):
        if base[0] != base[1]:
            seq += hetero(base[0],base[1])
        else:
            seq += base[0]
    return seq


def snpstring(array):
    snpsite = " "*9
    for site in array.transpose():
        reals = site.tolist()
        if len(set(reals)) > 1:
            " convert ambigs to reals"
            for i in xrange(len(reals)):
                if reals[i] in list("RWMSYK"):
                    for j in unstruct(reals[i]):
                        reals.append(j)
            reals = [i for i in reals if i not in "RWMSYK"]
            if sorted([reals.count(i) for i in set(reals)], reverse=True)[1] > 1:
                snpsite += "*"
            else:
                snpsite += "-"
        else:
            snpsite += " "
    snpsite += "|"
    return snpsite


def makeloci(aligns,names,outname):
    " output file "
    outloci = open(outname+".loci",'w')
    
    " make dictionary with list of loci for each sample "
    D = {}
    for i in set(names):
        D[i] = []
    for samp in D:
        for i in aligns:
            l = []
            for j in i:
                if j.name == samp:
                    l.append(j.sequence)
            D[samp].append(l) 

    " print in .loci readable format used by pyRAD " 
    nn = D.keys()
    nn.sort()

    if max([len(i) for i in nn]) > 9:
        print 'error: name lengths too long'; sys.exit(2)
    
    for i in range(params["Loci"]):
        l = []
        for n in nn:
            l.append(list(makehetero(D[n][i][0],D[n][i][1])))
            print >>outloci, ">"+n + " "*(10-len(n))+ makehetero(D[n][i][0],D[n][i][1])
        N = np.array(l)
        print >>outloci, "//"+snpstring(N)
    outloci.close()


def makephynex(outname, names, phy, nex):
    " order names "
    names = list(names)
    names.sort()
    longname = max(map(len, names))
    
    " read in loci file "
    finalfile = open(outname+".loci").read() 

    " dict for saving the full matrix "
    F = {}
    for name in names:
        F[name] = []

    " remove empty column sites and append edited seqs to dict F "
    for loc in [i for i in finalfile.split("|")[:-1]]:
        anames = [i.split(" ")[0][1:] for i in loc.strip().split("\n")[:-1]]
        array = np.array([tuple(i.split(" ")[-1]) for i in loc.strip().split("\n")][:-1])
        ## which columns are empty
        emps = [i for i in range(len(array.T)) if \
                np.all([j in ['N','-'] for j in array.T[i]])]
        ## delete those columns
        narray = np.delete(array, emps, 1)
        ## append data to dict
        for name in names:
            if name in anames:
                F[name] += "".join(narray[anames.index(name)])
            else:
                F[name] += "".join(["N"]*len(narray[0]))

    
    if phy:
        " print out .PHY file "
        superout = open(outname+".phy",'w')
        print >>superout, len(F), len("".join(F[names[0]]))
        for name in names:
            print >>superout, name+(" "*((longname+3)-len(name)))+"".join(F[name])
        superout.close()


    if nex:
        "print out .NEX file"
        nexout = open(outname+".nex", 'w')

        ntax = len(F)
        nchar = len(F.values()[0])

        print >>nexout, "#NEXUS"
        print >>nexout, "BEGIN DATA;"
        print >>nexout, "  DIMENSIONS NTAX=%s NCHAR=%s;" % (ntax,nchar)
        print >>nexout, "  FORMAT DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=YES;"
        print >>nexout, "  MATRIX"

        n=0
        sz = 100
        while n<len(F.values()[0]):
            for tax in F:
                print >>nexout, "  "+tax+" "*((longname-len(tax))+3)+"".join(F[tax][n:n+sz])
            n += sz
            print >>nexout, ""
        print >>nexout, ';'
        print >>nexout, 'END;'
        nexout.close()


def makemigrate(outname,tiptax,inds):
    outfile = open(outname+".migrate", 'w')

    ## cleanup taxadict
    taxa = {}
    for group in tiptax:
        taxa[group] = []
        for samp in range(inds):
            taxa[group].append(group+str(samp))

    ## read in data to sample names
    loci  = open(outname+".loci",'r').read().strip().split("|\n")[:]

    ## print data to file
    print >>outfile, len(taxa), len(loci), "( npops nloci for sim data", outname+".loci",")"
    
    ## print all data for each population at a time
    done = 0
    for group in taxa:
        ## print a list of lengths of each locus
        if not done:
            loclens = [len(loc.split("\n")[0].split(" ")[-1]) for loc in loci]
            print >>outfile, " ".join(map(str,loclens))
            done += 1

        ## print a list of number of individuals in each locus
        indslist = []
        for loc in loci:
            samps = [i.split(" ")[0].replace(">","") for i in loc.split("\n") if ">" in i]
            inds = sum([i in samps for i in taxa[group]])
            indslist.append(inds)
        print >>outfile, " ".join(map(str,indslist)), group

        ## print sample id, spaces, and sequence data
        #for loc in range(len(keep)):
        for loc in range(len(loci)):
            seqs = [i.split(" ")[-1] for i in loci[loc].split("\n") if \
                    i.split(" ")[0].replace(">","") in taxa[group]]
            for i in range(len(seqs)):
                print >>outfile, group+"_i"+str(i)+(" "*(10-len(group+"_i"+str(i))))+seqs[i]
    outfile.close()


def makestructure(outname,names):
    outfile = open(outname+".str",'w')

    "randomly select one SNP from each variable locus"
    loci  = open(outname+".loci",'r').read().strip().split("|\n")[:]

    F = {}  # taxa dict
    S = {}  # SNP dict
    for name in list(names):
        F[name] = []
        S[name] = []
    longname = max(map(len,names))

    " for each locus select out the SNPs"
    for loc in loci:
        ns = []
        ss = []
        for line in loc.split("\n"):
            if ">" in line:
                ns.append(line.split(" ")[0].replace(">",""))
                ss.append(line.split(" ")[-1])
            if "//" in line:
                snp = [i[0]-(longname+9) for i in enumerate(line) if i[1] in list('*-')]

        " assign snps to S"
        if snp:
            rando = snp[np.random.randint(len(snp))]
            for tax in F:
                #print ss[ns.index(tax)][rando]
                S[tax].append(ss[ns.index(tax)][rando])

    " make structure"
    SF = F.keys()
    SF.sort()
    
    N = np.array([list(S[i]) for i in SF])
    namescol = list(itertools.chain( * [[i+(" "*(longname-len(i)+3)),
                                   i+(" "*(longname-len(i)+3))] for i in SF] ))
    " add blank columns "
    empty = np.array(["" for i in xrange(len(SF)*2)])
    OUT = np.array([namescol,empty,empty,empty,empty,empty,empty])
    for col in xrange(len(N[0])):
        l = N[:,col]
        h = [unstruct(j) for j in l]
        h = list(itertools.chain(*h))
        bases = list("ATGC")
        final = [bases.index(i) if i not in list("-N") else '-9' for i in h]
        OUT = np.vstack([OUT, np.array(final)])
    np.savetxt(outname+".str", OUT.transpose(), fmt="%s", delimiter="\t")
                              

def maketreemix(outname):
    None
    

if __name__=="__main__":

    "parse command-line options"
    argv = sys.argv[1:]
    smallflags = "L:l:N:u:i:t:m:o:f:s:v"
    bigflags = []
    try:
        opts, args = getopt.getopt(argv,smallflags,bigflags)
        if not opts:
            usage()
            sys.exit(2)
    except getopt.GetoptError:
        print "Incorrect options passed"
        usage()
        sys.exit()

    "get params"
    params = parseopts(opts)

    "print options to screen"
    lu = params["mu"]*params["length"]   ## per loc mutation rate/gen
    params["theta"] = 4*params["N"]*lu
    if params["verbose"] == True:
        for p in params:
            print p, params[p]

    "simulate data"
    aligns,tiptax = simdata(params)

    " append names to seqs"
    names = []
    for name in tiptax:
        for i in range(params["inds"]):
            names.append(name+str(i))
            names.append(name+str(i))

    " appends names for allele 1 vs. allele 2 for each diploid sample "
    for i in aligns:
        for s,j in zip(i,names):
            s.name = j

    "make default LOCI output "
    makeloci(aligns,names,params["outname"])

    "make NEX or PHY outputs "
    phy = "p" in params["outform"]
    nex = "n" in params["outform"]
    if phy or nex:
        makephynex(params["outname"],names,phy,nex)

    "make MIGRATE outputs "
    if 'm' in params["outform"]:
        makemigrate(params["outname"],tiptax,params["inds"])

    "make STRUCTURE outputs "
    if "k" in params["outform"]:
        makestructure(params["outname"],names)

    
