#!/usr/bin/env python2

"""
###############################################################
##  Simulate sequence data on an input tree w/ introgression ##
##  output in .loci format and others                        ##
##  author:  Deren Eaton                                     ##
##  contact: deren.eaton@yale.edu                            ##
##  date:    6/30/15                                         ##
##  version: 0.9a                                            ##
##                                                           ##
##  TODO: allow migration among non-tip taxa                 ##
###############################################################
"""

## load standard modules
import getopt
import itertools
import sys
import datetime

## load external modules
try: 
    import egglib
except ImportError:
    print "Python package Egglib not found"
try:
    import numpy as np
except ImportError:
    print "Python package Numpy not found"

def parseopts(opts):
    """ parses options to the command line """

    ## defaults
    params = {'loci':   100,
              'length': 100,
              'N':      int(1e5),
              'mu':     1e-9,
              'inds':   1,
              'tree':   "(((a:1.0,b:1.0):1.0,c:2.0):1.0,d:3.0);",
              'migs':   "",
              'outname': "out",
              'outform': "l",
              'seed1':  123456,
              'seed2':  987654,
              'verbose':False
              }

    for opt, arg in opts:
        if opt in ["-L"]:
            params["loci"] = int(float(arg))
        elif opt in ["-l"]:
            params["length"] = int(float(arg))
        elif opt in ["-u"]:
            params["mu"] = float(arg)
        elif opt in ["-N"]:
            params["N"] = int(float(arg))
        elif opt in ["-i"]:
            params["inds"] = int(arg)
        elif opt in ["-t"]:
            params["tree"] = str(arg)
        elif opt in ["-m"]:
            params["migs"] = str(arg)
        elif opt in ["-o"]:
            params["outname"] = str(arg)
        elif opt in ["-f"]:
            params["outform"] = str(arg)
        elif opt in ["-v"]:
            params["print"] = True
        elif opt in ["-1"]:
            params["seed1"] = int(arg)
        elif opt in ["-2"]:
            params["seed2"] = int(arg)
    return params


def checkopts(opts):
    """ print errors if options are incorrect"""
    ## should move tree checking functions here...
    print opts


def usage():
    """
    brief description of various flags and options for this script
    """
    print "\nHere is how you can use this script\n"
    print "Usage: python %s"%sys.argv[0]
    print "\t -o <str>   output file prefix name (required)"
    print "\t -L <int>   Number of loci to simulate (default 100) "
    print "\t -l <int>   length of loci (default 100) "
    print "\t -u <float> per-site mutation rate (default: 1e-9) "
    print "\t -N <int>   effective population size (default: 1e5) "
    print "\t -i <int>   number of sampled inds per tip taxon (default 1) "
    print "\t -t <str>   ultrametric tree w br lens (file or Newick string)"
    print "\t -s1 <int>  random seed 1 (default 123456)"
    print "\t -s2 <int>  random seed 2 (default 987654)"
    #print "\t -v         print version"
    #print "\t -p <str>   print verbose : prints all params to screen"
    print "\t -f <str>   output file format {l=loci (default), \n"+\
    "                                        p=phy, n=nex,\n"+\
    "                                        a=alleles, v=vcf,\n"+\
    "                                        t=treemix, k=structure, \n"+\
    "                                        g=geno, m=migrate} "
    print "\t -m <str>   migration events in string format:\n"+\
    "                      e.g., m=5e-5 from b->c \n"+\
    "                      (forward in time) from time=0.5 to 0.25:\n"+\
    "                      -m [c,b,0.5,0.25,5e-5]"
    #print "\t -pt    prints tree info to screen (requires ete2 or biopython)"
    

def lengthtotip(node):
    """ returns node height """
    dec = 0.
    while node.descendants():
        dec += node.branch_to(node.descendants()[0])
        node = node.descendants()[0]
    return dec

def is_ultrametric(tree):
    """ check if tree is ultrametric """
    tips = tree.get_terminal_nodes()
    dist_to_root = []
    for node in tips:
        leaf = node.get_label()
        edgesums = 0
        while node.ascendants():
            edgesums += node.branch_from(node.ascendants()[0])
            node = node.ascendants()[0]
        dist_to_root.append((leaf, edgesums))
    if len(set([i[1] for i in dist_to_root])) == 1:
        return 1
    else:
        print "\nError: tree is not ultrametric"
        for node in dist_to_root:
            print node[0], "dist_to_root:", node[1]
        sys.exit(2)


def simdata(params):
    "split multiple migrations"
    migscenarios = [m.lstrip("[").rstrip("]") for \
                        m in params["migs"].split('/')]

    try: 
        tre = egglib.Tree(string=params["tree"])
    except IOError:
        try:
            tre = egglib.Tree(fname=params["tree"])
        except IOError:
            print "\n\tInvalid newick tree input (-t)"
            sys.exit(2)
    ## check that tree has branch lengths
    try: 
        tre.total_length()
    except ValueError: 
        print "\n\tInvalid newick string, tree must have branch lengths"
        sys.exit(2)
    ## check that tree is ultrametric
    is_ultrametric(tre)

    ## get terminal names
    tiptax = tre.all_leaves()

    ## sets the two parameter classes
    paramset = egglib.simul.CoalesceParamSet(
                                singleSamples=None,
                                doubleSamples=[params["inds"]]*len(tiptax),
                                M=0.0)

    ## traverse tree fusing populations at each node
    if params["tree"]:
        for node in tre:
            if node.descendants():
                date = lengthtotip(node)
                tip1 = min([tiptax.index(l) for l in \
                            node.descendants()[0].leaves_down()])
                tip2 = min([tiptax.index(l) for l in \
                            node.descendants()[1].leaves_down()])
                paramset.populationFusion(date, tip1, tip2)


    ## sets migration patterns
    ## ADD MIGRATION BETWEEN ANCESTRAL BRANCHES???
    ## LOOKING FORWARD IN TIME
    if params["migs"]:
        for mig in migscenarios:
            recipient, donor, start, end, migration = mig.split(',')
            if float(start) < float(end):
                sys.exit("\nError: migration start time should be greater"+\
                          " than migration end time (present=0) ")
            irec = tiptax.index(recipient)
            idon = tiptax.index(donor)
            imig = 4.*params["N"]*float(migration)
            paramset.changePairwiseMigrationRate(0.0, int(irec), int(idon), 0.)
            paramset.changePairwiseMigrationRate(float(start), int(irec), 
                                                   int(idon), float(imig))
            paramset.changePairwiseMigrationRate(float(end), int(irec),
                                                   int(idon), 0.)

    mutator = egglib.simul.CoalesceFiniteAlleleMutator(
                                            theta=params["theta"],
                                            alleles=4,
                                            randomAncestralState=True)
    mutator.setSites(params["length"])
    aligns = egglib.simul.coalesce(paramset,
                                   mutator,
                                   params["loci"],
                                   random=(params["seed1"],
                                           params["seed2"]))

    return [aligns, tiptax]


def hetero(base1, base2):
    """ returns IUPAC symbol for ambiguity bases,
    used for polymorphic sites. """
    Basedict = {('G','A'):"R",
                ('G','T'):"K", 
                ('G','C'):"S", 
                ('T','C'):"Y", 
                ('T','A'):"W", 
                ('C','A'):"M"}
    ambig1 = Basedict.get((base1, base2)) 
    ambig2 = Basedict.get((base2, base1))
    if ambig1:
        return ambig1
    else:
        return ambig2


def unstruct(amb):
    " returns bases from ambiguity code"
    amb = amb.upper()
    ambdict = {"R":["G","A"],
               "K":["G","T"],
               "S":["G","C"],
               "Y":["T","C"],
               "W":["T","A"],
               "M":["C","A"]}
    if amb in ambdict:
        return ambdict.get(amb)
    else:
        return amb


def makehetero(seq1, seq2):
    """ take two alleles and return one seqstring with
    ambiguities for heterozygous sites"""
    seq = ""
    for base in zip(seq1, seq2):
        if base[0] != base[1]:
            seq += hetero(base[0], base[1])
        else:
            seq += base[0]
    return seq


def snpstring(array, count, longname):
    """ make a string of snps from an array for .loci output """
    snpsite = " "*(longname+2)
    for site in array.transpose():
        reals = site.tolist()
        if len(set(reals)) > 1:
            ## convert ambigs to reals"
            for base in xrange(len(reals)):
                if reals[base] in list("RWMSYK"):
                    for allele in unstruct(reals[base]):
                        reals.append(allele)
            reals = [i for i in reals if i not in "RWMSYK"]
            if sorted([reals.count(i) for i in set(reals)], reverse=True)[1] > 1:
                snpsite += "*"
            else:
                snpsite += "-"
        else:
            snpsite += " "
    snpsite += "|"+str(count)+"|"
    return snpsite


def makeloci(aligns, names, nloci, outname, outform):
    """ format data for output in .loci format """
    ## output file 
    outloci = open(outname+".loci", 'w')
    
    ## make dictionary with list of loci for each sample
    ldict = {i:[] for i in set(names)}
    for sampname in ldict:
        for loc in aligns:
            seqlist = []
            for seqobj in loc:
                if seqobj.name == sampname:
                    seqlist.append(seqobj.sequence)
            ldict[sampname].append(seqlist) 

    ## get names in sorted order
    dicnames = ldict.keys()
    dicnames.sort()

    ## get longest name length.
    #if max([len(n) for n in dicnames]) > 9:
        #print 'error: name lengths too long'
        #sys.exit(2)
    longname = max([len(n) for n in dicnames])
    
    ## print each locus to file with a snpstring
    for locus in range(nloci):
        seqstring = []
        for dicname in dicnames:
            ## make a single seq string with ambiguity chars
            seqstring.append(list(makehetero(ldict[dicname][locus][0],
                                             ldict[dicname][locus][1])))
            print >>outloci, ">"+dicname+" "*(longname-len(dicname)+3)+\
                                           "".join(seqstring[-1])

        ## print snpstring to final line of locus
        print >>outloci, "//"+snpstring(np.array(seqstring),
                                        locus, longname) 
    outloci.close()

    ## .alleles output
    if 'a' in outform:
        outalleles = open(outname+".alleles", 'w')

        for locus in range(nloci):
            seqstring = []
            for dicname in dicnames:
                seqstring.append(list(makehetero(ldict[dicname][locus][0],
                                                 ldict[dicname][locus][1])))
                print >>outalleles, ">"+dicname+"_0"+\
                                    " "*(longname-len(dicname)+3)+\
                                    ldict[dicname][locus][0]
                print >>outalleles, ">"+dicname+"_1"+\
                                    " "*(longname-len(dicname)+3)+\
                                    ldict[dicname][locus][1]

            ## print snpstring to final line of locus
            print >>outalleles, "//  "+snpstring(np.array(seqstring),
                                                 locus, longname) 
        outalleles.close()



def makephynex(outname, names, phy, nex):
    " order names "
    names = list(names)
    names.sort()
    longname = max([len(n) for n in names])
    
    ## read in loci file
    finalfile = open(outname+".loci").read() 

    ## dict for saving the full matrix
    fdict = {name:[] for name in names}

    ## remove empty column sites and append edited seqs to Fdict
    loci = finalfile.split("|\n")[:-1]
    for loc in loci:
        anames = [i.split()[0].replace(">", "") for \
                    i in loc.split("\n") if ">" in i]
        array = np.array(tuple([i.split()[-1].replace(">", "") for \
                    i in loc.split("\n") if ">" in i]))
        ## which columns are empty
        emps = [i for i in range(len(array.T)) if \
                np.all([j in ['N', '-'] for j in array.T[i]])]
        ## delete those columns
        if emps:
            narray = np.delete(array, emps, 1)
        else:
            narray = array
        ## append data to dict
        for name in set(names):
            if name in anames:
                fdict[name] += "".join(narray[anames.index(name)])
            else:
                fdict[name] += "".join(["N"]*len(array[0]))

    if phy:
        ## print out .PHY file
        superout = open(outname+".phy", 'w')
        print >>superout, len(fdict), len("".join(fdict[names[0]]))
        duplicate_names = []
        for name in names:
            if name not in duplicate_names:
                print >>superout, name+(" "*((longname+3)-\
                                  len(name)))+"".join(fdict[name])
                duplicate_names.append(name)
        superout.close()

    if nex:
        ## print out .NEX file
        nexout = open(outname+".nex", 'w')

        ntax = len(fdict)
        nchar = len(fdict.values()[0])

        print >>nexout, "#NEXUS"
        print >>nexout, "BEGIN DATA;"
        print >>nexout, "  DIMENSIONS NTAX=%s NCHAR=%s;" % (ntax, nchar)
        print >>nexout, "  FORMAT DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=YES;"
        print >>nexout, "  MATRIX"

        ## print 100 chars at a time (interleaved format)
        blockstart = 0
        block = 100
        while blockstart < len(fdict.values()[0]):
            duplicate_names = []
            for tax in names:
                if tax not in duplicate_names:
                    print >>nexout, "  "+tax+" "*((longname-len(tax))+3)+\
                                "".join(fdict[tax][blockstart:blockstart+block])
                    duplicate_names.append(tax)
            blockstart += block
            print >>nexout, ""
        print >>nexout, ';'
        print >>nexout, 'END;'
        nexout.close()


def makemigrate(outname, tiptax, inds):
    outfile = open(outname+".migrate", 'w')

    ## cleanup taxadict
    taxa = {}
    for group in tiptax:
        taxa[group] = []
        for samp in range(inds):
            taxa[group].append(group+str(samp))
    taxalist = taxa.keys()
    taxalist.sort()

    ## read in data to sample names
    loci = open(outname+".loci", 'r').read().strip().split("//")[:-1]

    ## print header to file
    print >>outfile, len(taxa), len(loci), "("+outname+".loci", ")"
    
    ## print all data for each population at a time
    done = 0
    for group in taxalist:
        ## print a list of lengths of each locus
        if not done:
            loclens = [len(loc.split("\n")[1].split()[-1]) \
                         for loc in loci if ">" in loc]
            print >>outfile, " ".join([str(i) for i in loclens])
            done += 1

        ## print a list of number of individuals in each locus
        indslist = []
        for loc in loci:
            samps = [i.split()[0].replace(">", "") for \
                       i in loc.split("\n") if ">" in i]
            inds = sum([i in samps for i in taxa[group]])
            indslist.append(inds)
        print >>outfile, " ".join([str(i) for i in indslist]), group

        ## print sample id, spaces, and sequence data
        #for loc in range(len(keep)):
        for loc in range(len(loci)):
            seqs = [i.split(' ')[-1] for i in loci[loc].split("\n") if \
                    i.split(' ')[0].replace(">", "") in taxa[group]]
            for i in range(len(seqs)):
                print >>outfile, group+"_i"+str(i)+\
                          (" "*(10-len(group+"_i"+str(i))))+seqs[i]
    outfile.close()


def makestructure(outname):
    """ format data for structure output file """
    print "\nstructure output not implemented yet"
    # ## randomly select one SNP from each variable locus
    # loci = open(outname+".loci", 'r').read().strip().split("//")[:-1]
    # taxa_dict = {}
    # snp_dict = {}
    # for name in list(names):
    #     taxa_dict[name] = []
    #     snp_dict[name] = []
    # longname = max([len(n) for n in names])

    # ## for each locus select out the SNPs"
    # for loc in loci:
    #     nameseq = []
    #     seqsseq = []
    #     for line in loc.split("\n"):
    #         if ">" in line:
    #             nameseq.append(line.split(" ")[0].replace(">", ""))
    #             seqsseq.append(line.split(" ")[-1])
    #         if "//" in line:
    #             snp = [j[0]-(longname+9) for i, j in \
    #                     enumerate(line) if j[1] in list('*-')]

    #     ## assign snps to S
    #     if snp:
    #         rando = snp[np.random.randint(len(snp))]
    #         for tax in F:
    #             #print ss[ns.index(tax)][rando]
    #             S[tax].append(ss[ns.index(tax)][rando])

    # " make structure"
    # SF = F.keys()
    # SF.sort()
    
    # N = np.array([list(S[i]) for i in SF])
    # namescol = list(itertools.chain( * [[i+(" "*(longname-len(i)+3)),
    #                 i+(" "*(longname-len(i)+3))] for i in SF] ))
    # " add blank columns "
    # empty = np.array(["" for i in xrange(len(SF)*2)])
    # OUT = np.array([namescol,empty,empty,empty,empty,empty,empty])
    # for col in xrange(len(N[0])):
    #     l = N[:,col]
    #     h = [unstruct(j) for j in l]
    #     h = list(itertools.chain(*h))
    #     bases = list("ATGC")
    #     final = [bases.index(i) if i not in list("-N") else '-9' for i in h]
    #     OUT = np.vstack([OUT, np.array(final)])
    # np.savetxt(outname+".str", OUT.transpose(), fmt="%s", delimiter="\t")
                              

def maketreemix(outname):
    """ format data for treemix output file """
    print '\ntreemix format not implemented yet'

def makegeno(outname):
    """ format data for treemix output file """
    print '\nGENO format not implemented yet'

def makeVCF(outname):
    """ format data for treemix output file """
    print '\nVCF format not implemented yet'


if __name__ == "__main__":

    ## parse command-line options
    ARGS = sys.argv[1:]
    SMALLFLAGS = "L:l:N:u:i:t:m:o:f:s1:s2:v:"
    try:
        OPTS, ARGS = getopt.getopt(ARGS, SMALLFLAGS)
        if not OPTS:
            usage()
            sys.exit(2)
    except getopt.GetoptError:
        print "Incorrect options passed"
        usage()
        sys.exit()

    ## get params
    PARAMS = parseopts(OPTS)

    ## print log file 
    LOGFILE = open(PARAMS["outname"]+".log", 'w')
    LU = PARAMS["mu"]*PARAMS["length"]   ## per loc mutation rate/gen
    PARAMS["theta"] = 4*PARAMS["N"]*LU
    PARAMKEYS = PARAMS.keys()
    PARAMKEYS.sort()
    LONGKEY = max([len(k) for k in PARAMKEYS])

    ## print header
    print >>LOGFILE, PARAMS["outname"]
    print >>LOGFILE, str(str(datetime.datetime.now()).split(".")[0])
    print >>LOGFILE, "------------------------------------------"    
    if PARAMS['verbose']:
        print >>sys.stderr, "\n--------------------------------------"    
        print >>sys.stderr, PARAMS["outname"]+"\t"+\
                            str(str(datetime.datetime.now()).split(".")[0])
        print >>sys.stderr, "--------------------------------------"    

    ## print params to file or screen
    for param in PARAMKEYS:
        if PARAMS["verbose"]:
            print >>sys.stderr, param+" "*(LONGKEY-len(param)+4)+" "+\
                                str(PARAMS[param])
        print >>LOGFILE, param+" "*(LONGKEY-len(param)+4)+" "+\
                         str(PARAMS[param])
    LOGFILE.close()

    ## simulate data
    ALIGNS, TIPTAX = simdata(PARAMS)

    ## append names to seqs"
    NAMES = []
    for tip in TIPTAX:
        if PARAMS["inds"] == 1:
            NAMES.append(tip)
            NAMES.append(tip)
        else:
            for ind in range(PARAMS["inds"]):
                NAMES.append(tip+str(ind))
                NAMES.append(tip+str(ind))

    ## appends names for each allele copy in diploid sample 
    for alignobj in ALIGNS:
        for copy, cname in zip(alignobj, NAMES):
            copy.name = cname

    ## make default LOCI output 
    makeloci(ALIGNS, NAMES, PARAMS["loci"], 
             PARAMS["outname"], PARAMS["outform"])

    ## make NEX or PHY outputs 
    PHY = "p" in PARAMS["outform"]
    NEX = "n" in PARAMS["outform"]
    if PHY or NEX:
        makephynex(PARAMS["outname"], NAMES, PHY, NEX)

    ## make MIGRATE outputs 
    if 'm' in PARAMS["outform"]:
        makemigrate(PARAMS["outname"], TIPTAX, PARAMS["inds"])

    ## make STRUCTURE outputs 
    if "k" in PARAMS["outform"]:
        makestructure(PARAMS["outname"])

    ## make TREEMIX outputs 
    if "t" in PARAMS["outform"]:
        maketreemix(PARAMS["outname"])

    ## make GENO outputs 
    if "g" in PARAMS["outform"]:
        makegeno(PARAMS["outname"])

    ## make GENO outputs 
    if "v" in PARAMS["outform"]:
        makeVCF(PARAMS["outname"])

    
