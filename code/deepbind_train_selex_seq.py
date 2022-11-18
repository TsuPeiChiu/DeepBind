import sys
import re
import os
import os.path
import argparse
import numpy as np
import numpy.random as npr
import gzip
import deepbind_util as util

import cPickle
import scipy
import shutil

datadir = "../data/selex_seq"


def main():
    util.enable_reversecomplement()

    args   = loadargs()
    models = loadmodels(args)
    trdata = None
    tedata = None
    tfid = args.id[0]
    flen = 8

    if "calib" in args.steps:
        print "-------------- calib:", tfid, "--------------"
        set_motif_lengths(args, models, flen)
        trdata = load_trdata(trdata, tfid, None, args)
        print "args.ncalib:" + str(args.ncalib)
        util.calibrate(models, trdata, args.calibdir, nfold=args.nfold, ncalib=args.ncalib, allfolds=True)

    if "train" in args.steps:
        print "-------------- train:", tfid, "--------------"
        print "args.ntrial:" + str(args.ntrial)
        for i in range(args.ntrial):
            set_motif_lengths(args, models, flen)
            trdata = load_trdata(trdata, tfid, None, args)
            finaldir = os.path.join(args.finaldir, str(i))
            util.train(models, trdata, args.calibdir, finaldir, nfold=1, ntrial=1, metric_key="pearson.r")

            if "test" in args.steps:
                print "-------------- test:", tfid, "--------------"
                tedata = load_tedata(tedata, tfid, args)
                util.save_metrics(tedata, "test", finaldir)

    print "-------------- collect pearson.r for trials:", tfid, "--------------"
    with open(os.path.join(args.finaldir, "pearson.r.txt"), 'w') as fw:
        for i in range(args.ntrial):
            metricdir = os.path.join(args.finaldir, str(i), "relKa", "metrics.txt")
            with open(metricdir, "r") as fr:
                lines = fr.readlines()
                print(lines[6].strip())
                fw.write(lines[6])

#########################################################################

def loadargs():
    global datadir
    args = argparse.ArgumentParser(description="Generate the SELEX_seq experiments.")
    args = util.parseargs("selex_seq", args)

    args.calibdir  = "/".join([args.outdir, args.id[0], "calib"])
    args.finaldir  = "/".join([args.outdir, args.id[0], "final"])
    args.reportdir = "/".join([args.outdir, args.id[0], "report"])
    args.testdir   = "/".join([args.outdir, args.id[0], "test"])
    
    args.ncalib    = 100 # default: 30
    args.ntrial    = 10 # default: 6

    return args

#################################

def loadmodels(args):
    models = util.loadmodels(args, "cfg/regression/maxpool")
    return models

#################################

def set_motif_lengths(args, models, flen):
    print "Using filter size %d" % flen
    for cfg in models.itervalues():
        cfg["model"].conv_seq[0].fsize = flen

#################################

def load_trdata(trdata, tfid, fold, args):
    if trdata is not None:
        return trdata  # Training data was already loaded in earlier step.
    maxrows = 10000 if args.quick else None
    
    if tfid == 'max':
        data = util.datasource.fromtxt("../data/selex_seq/r0_r1_max_selex_seq_10mer_150_onestrand_seq_train.txt", 
                None, "../data/selex_seq/r0_r1_max_selex_seq_10mer_150_onestrand_relKa_train.txt", targetcols=[0], maxrows=maxrows)
    
    elif tfid == 'mef2b':
        data = util.datasource.fromtxt("../data/selex_seq/r0_r1_mef2b_selex_seq_10mer_100_onestrand_seq_train.txt",
                None, "../data/selex_seq/r0_r1_mef2b_selex_seq_10mer_100_onestrand_relKa_train.txt", targetcols=[0], maxrows=maxrows)
    
    elif tfid == 'p53':
        data = util.datasource.fromtxt("../data/selex_seq/r0_r1_p53_selex_seq_10mer_300_onestrand_seq_train.txt",
                None, "../data/selex_seq/r0_r1_p53_selex_seq_10mer_300_onestrand_relKa_train.txt", targetcols=[0], maxrows=maxrows)
    
    return data

#################################

def load_tedata(tedata, tfid, args):
    if tedata is not None:
        return tedata  # Test data was already loaded in earlier step.
    maxrows = 10000 if args.quick else None

    if tfid == 'max':
        data = util.datasource.fromtxt("../data/selex_seq/r0_r1_max_selex_seq_10mer_150_onestrand_seq_val.txt",
                None, "../data/selex_seq/r0_r1_max_selex_seq_10mer_150_onestrand_relKa_val.txt", targetcols=[0], maxrows=maxrows)

    elif tfid == 'mef2b':
        data = util.datasource.fromtxt("../data/selex_seq/r0_r1_mef2b_selex_seq_10mer_100_onestrand_seq_val.txt",
                None, "../data/selex_seq/r0_r1_mef2b_selex_seq_10mer_100_onestrand_relKa_val.txt", targetcols=[0], maxrows=maxrows)

    elif tfid == 'p53':
        data = util.datasource.fromtxt("../data/selex_seq/r0_r1_p53_selex_seq_10mer_300_onestrand_seq_val.txt",
                None, "../data/selex_seq/r0_r1_p53_selex_seq_10mer_300_onestrand_relKa_val.txt", targetcols=[0], maxrows=maxrows)

    return data

#################################



if __name__=="__main__":
    #util.disable_nultiprocessing()
    main()

