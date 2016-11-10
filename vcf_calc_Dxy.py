#!/usr/bin/env python
"""
this will be a fork of travc script to calculate Dxy from a vcf file. All credit for the 
functions go to him.
"""
import sys
import os
import time
import argparse
import glob
import collections
import ConfigParser
import shlex
import pipes
import subprocess
import tempfile
import re
import errno
import fileinput

import numpy as NP
import vcf # see: http://pyvcf.readthedocs.org/en/latest/

### Constants #################################################################

DEFAULT_CONFIG = 'Dxy_test.cfg' # Default config filename; None for none

################################################################################

### ConfigParser helper stuff (doesn't normally need any changing) #############

def printArgs(args, out_file_handle=sys.stderr, prefix=''):
    """print out the args (Namespace) created by argparse
       'prefix' is just added to beginning of each line"""
    for k,v in vars(args).iteritems() :
        print >>out_file_handle, prefix+str(k),":",v

class ConfigFakeSecHead(object):
    def __init__(self, fp, section='DEFAULTS'):
        self.fp = fp
        self.sechead = '['+str(section)+']\n'
    def readline(self):
        if self.sechead:
            try: return self.sechead
            finally: self.sechead = None
        else: return self.fp.readline()

class CustomArgparseHelpFormatter(argparse.HelpFormatter):
    """Help message formatter for argparse
        combining RawTextHelpFormatter and ArgumentDefaultsHelpFormatter
    """

    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])

    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (default: %(default)s)'
        return help

## End: ConfigParser helper stuff
################################################################################


### Main #######################################################################
def exit(retcode=0):
    # cleanup and exit
    global g_start_tic
    print >>sys.stderr, "END", "%0.2f"%(time.time()-g_start_tic),"secs elapsed"
    sys.exit(retcode)


def Main(argv=None):
    #print >>sys.stderr, "START"
    global g_start_tic
    g_start_tic = time.time()

    ### Handling the arguments/options/config file ###
    # parse cfg_file argument
    conf_parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=CustomArgparseHelpFormatter,
                                        add_help=False) # turn off help so later parse (with all opts) handles it
    conf_parser.add_argument('-c', '--cfg-file', type=argparse.FileType('r'),
            help="Config file specifiying options/parameters.\nAny long option can be set by remove the leading '--' and replace '-' with '_'",
            default=DEFAULT_CONFIG)
    args, remaining_argv = conf_parser.parse_known_args(argv)

    # build the config (read config files)
    if args.cfg_file:
        cfg = ConfigParser.SafeConfigParser()
        cfg.optionxform=str # make the ConfigParser options case sensitive
        cfg.readfp(ConfigFakeSecHead(args.cfg_file))
        defaults = dict(cfg.items("DEFAULTS"))
        # special handling of paratmeters that need it like lists
        for k in ['x_samples', 'y_samples']:
            if( k in defaults ):
                defaults[k] = [ x for x in defaults[k].split('\n') if x and x.strip() and not x.strip()[0] in ['#',';'] ]
    else:
        defaults = {}

    # Parse rest of arguments with a new ArgumentParser
    aparser = argparse.ArgumentParser(description=__doc__, parents=[conf_parser], formatter_class=CustomArgparseHelpFormatter)

    # input and output
    aparser.add_argument('-i', '--in-vcf', type=argparse.FileType('r'), default=sys.stdin)

    # populations/groups
    aparser.add_argument('-x', '--x-samples', help='x sample name list')
    aparser.add_argument('-y', '--y-samples', help='y sample name list')

    # filtering otpions
    aparser.add_argument('--allow-non-snp', action='store_true', help='Allow variants other than just SNPs')
    aparser.add_argument('--allow-non-biallelic', action='store_true', help='Allow non-biallelic variants')
    aparser.add_argument('--max-domAF', type=float,
                help="Max allowable frequency of the dominant (most frequent) allele to include a locus", default=None)
    # per sample filtering options
    aparser.add_argument('-q', '--min-GQ', type=int,
                help="Min allowable genotype quality (GQ) score of a sample at a locus", default=50)
    aparser.add_argument('-d', '--min-depth', type=int,
                help="Min allowable coverage depth of a sample at a locus", default=2)
    aparser.add_argument('-D', '--max-depth', type=int,
                help="Max allowable coverage depth of a sample at a locus", default=1000)
    # windowing parameters (include overall start and end positions)
    aparser.add_argument('-w', '--window-size', help="Size of sliding window", default='10e3')
    aparser.add_argument('-s', '--window-step', help="Step between sliding windows", default='5e3')
    aparser.add_argument('--start', help="Ignore any variants at positions before this", default='1')
    aparser.add_argument('--end', help="Ignore any variants at positions after this", default=None)
    aparser.add_argument('-N', '--min-num-points-in-window', type=int, help="Ignore windows with fewer than this many point", default='1')

    # general option (like verbose)
    aparser.add_argument('-v', '--verbose', action='count', help="Increase verbosity level")

    aparser.set_defaults(**defaults) # applies any values read from the config file by setting them as defaults

    # process options/arguments
    args = aparser.parse_args(remaining_argv)

    # custom/more complex argument parsing errors
    # numbers which are given as strings so stuff like '1e5' works
    for arg_name in ['window_size', 'window_step', 'start', 'end', ]:
        if( getattr(args, arg_name) is not None ):
            if( getattr(args, arg_name).lower() == 'none' ):
                setattr(args, arg_name, None)
            else:
                setattr(args, arg_name, int(round(float(getattr(args, arg_name)))))


    # optionally output the configuration/arguments
    if( args.verbose > 0 ):
        print >>sys.stderr, "CONFIGURATION OPTIONS/ARGUMENTS:"
        printArgs(args, prefix='\t')

    ## Done: Handling the arguments/options/config file

    ### Actually do the work ###
    dxy_parser = DxyFromVcf(args) # reads and filters the vcf

    for win,data in sparseSlidingWindow(
                    dxy_parser,
                    size=args.window_size,
                    step=args.window_step,
                    start=args.start,
                    end=args.end):

        #print win, data
        dxy = NP.array([x[1] for x in data])

        # @TCC this multiple print is ugly, but at least explicit... try to find a better way
        print (win[0]+win[1]-1)/2.0, # midpoint of window
        print len(dxy), # num points

        if( len(data) < args.min_num_points_in_window ):
            print ' '.join(['nan']*3),
        else:
            print NP.mean(dxy),
            print NP.std(dxy),
            print NP.median(dxy),
        print
        sys.stdout.flush()

    ### Cleanup and end normally ###
    exit()


## functions ############################################################

def DxyFromVcf(args):
    """ iterater
    """
    # @TCC TODO Make make this a generator... that would be cool (and more efficient memory wise)
    # used from args... (listed here to be explicit)
    verbose = args.verbose
    in_vcf = "/scratch365/ssmall2/Wb_analysis/analysis/dxy/Wb.10k.54.recode.vcf"
    x_samples = args.x_samples
    y_samples = args.y_samples
    # per locus thresholds/filter values
    filter_non_snp = not args.allow_non_snp
    filter_non_biallelic = not args.allow_non_biallelic
    # per sample thresholds/filter values
    min_GQ = args.min_GQ
    min_depth = args.min_depth
    max_depth = args.max_depth
    max_domAF = args.max_domAF
    # constants which may be options later
    info_outfp = sys.stderr # where informative and summary info goes (must be set)

    # create the vcf.Reader object; pyVCF is way cool: http://pyvcf.readthedocs.org/en/latest/
    
    vcf_reader = vcf.Reader(open(in_vcf))
    print vcf_reader.samples
    # verify samples in x and y groups are actually in the vcf
    for samp in set(x_samples+y_samples):
        if( samp not in vcf_reader.samples ):
            raise Exception("sample '{}' not in vcf file '{}'".format(samp, "vcf"))

    # bookkeeping (including purely informative)
    num_sites = 0
    num_sites_kept = 0
    filter_counts = collections.OrderedDict([ # 1st col of values is count, 2nd col is just for informative output
        ('monomorphic', [0]),
        ('non_snp', [0]),
        ('non_biallelic', [0]),
        ('max_domAF', [0, max_domAF]),
        ('samp_no_GT_call', [0]),
        ('samp_min_depth', [0, min_depth]),
        ('samp_max_depth', [0, max_depth]),
        ('samp_min_GQ', [0, min_GQ]),
        ])
    total_filtered = 0

    for rec in vcf_reader:
        num_sites += 1

        # filtering
        if( rec.is_monomorphic ):
            filter_counts['monomorphic'][0] += 1
            total_filtered += 1
            continue
        if( filter_non_snp and not rec.is_snp ): # SNP only (@TCC make an option)
            filter_counts['non_snp'][0] += 1
            total_filtered += 1
            continue
        if( filter_non_biallelic and len(rec.alleles) != 2 ): # biallelic only (@TCC make an option?)
            # @TCC TODO rejects if biallelic but different from reference; We can be more clever
            filter_counts['non_biallelic'][0] += 1
            total_filtered += 1
            continue
        # allele frequencies (needed for domAF filter, ect.)
        allele_freq = list(rec.aaf)
        allele_freq.insert(0, 1-sum(rec.aaf))
        if( max_domAF is not None and max(allele_freq) > max_domAF ):
            filter_counts['max_domAF'][0] += 1
            total_filtered += 1
            continue

        # per-sample filters (@TCC requires all requested samples to pass all filters... Maybe relax that)
        filter_rec_flag = False
        for samp in set(x_samples+y_samples):
            c = rec.genotype(samp)
            if( not c.called ):
                filter_counts['samp_no_GT_call'][0] += 1
                filter_rec_flag = True
                break
            if( c['DP'] < min_depth ):
                filter_counts['samp_min_depth'][0] += 1
                filter_rec_flag = True
                break
            if( c['DP'] > max_depth ):
                filter_counts['samp_max_depth'][0] += 1
                filter_rec_flag = True
                break
            if( c['GQ'] < min_GQ ):
                filter_counts['samp_min_GQ'][0] += 1
                filter_rec_flag = True
                break
        if( filter_rec_flag ):
            total_filtered += 1
            continue

        # passed filtering...
        num_sites_kept += 1

        ## calculations ##

        ## Dxy ##############################################################
        # Frequency based Dxy formula provided Jacob Crawford
        # @TCC THIS ASSUMES BIALLLELIC DIPLOID ...
        # $$ D_{xy} = \sum_1^S{\frac{ (f_1 n_1 (1-f_2) n_2) + ((1-f_1) n_1 f_2 n_2) }{n_1 n_2}} $$
        # $f_i$ is major allele freq in pop $i$ (Jacob had `minor', but it is equivalent) \\
        # $n_i$ is number of dipolid samples in pop $i$ (for the site being evaluated) \\
        # $S$ is the num of sites with data for both pops... we are outputting Dxy per-site, so $S=1$
        # Actual formula being used is based on allele counts $c_i = 2 f_i n_i$
        # $$ D_{xy} = \sum{ \frac{ c_1 ( 2 n_2 - c_2 ) + (2 n_1 - c_1 ) c_2 }{ 4 n_1 n_2 } }$$

        if( len(rec.alleles) != 2 ): # ensure biallelic
            raise Exception("Can only compute Dxy using biallelic sites")
        num_pops = 2 # @TCC ONLY 2 WORKS (for now)
        pop_allele_count = [[]]*num_pops # per pop, per allele
        total_allele_count = [0]*(len(rec.alleles)) # per allele
        n = [0]*num_pops # $n_i$
        pop = [x_samples, y_samples] # pop$_i$
        for i in range(num_pops):
            pop_allele_count[i] = [0]*(len(rec.alleles))
            for samp in pop[i]:
                n[i] += 1 # DIPLOID
                if( len(rec.genotype(samp).gt_alleles) != 2 ):
                    raise Exception("Can only compute Dxy on diploid sites/samples")
                for allele in rec.genotype(samp).gt_alleles:
                    pop_allele_count[i][int(allele)] += 1
                    total_allele_count[int(allele)] += 1
        # determine which allele is `major' in these samples
        major_allele_idx = max( (v,i) for i,v in enumerate(total_allele_count) )[1]
        # get the count of the major allele occurance in each pop
        c = [ pop_allele_count[i][major_allele_idx] for i in range(num_pops) ]
        # finally, the Dxy formula
        Dxy = ( c[0]*(2*n[1]-c[1]) + (2*n[0]-c[0])*c[1] ) / float( 4*n[0]*n[1] )
        # the frequency (as opposed to count) based formula... (also works)
        #f = [ (c[i] / float(2*n[i])) for i in range(num_pops) ]
        #Dxy_freq = ( f[0]*n[0]*(1-f[1])*n[1] + (1-f[0])*n[0]*f[1]*n[1] ) / float( n[0]*n[1] )

        # @TCC TODO Implement a more robust (brute-force) Dxy calc

        yield (rec.POS, Dxy)
        #data_out.append([rec.POS, Dxy])
    # done going through vcf records #

    ##
    # informative/summary output
    print >>info_outfp, "Filtered (applied in order, so eg: non-SNP isn't tested for non-biallelic) --"
    # pretty output of the filter counts
    tmp = []
    for k,v in filter_counts.iteritems():
        if( len(v) > 1 ):
            tmp.append(["{}.({})".format(k,v[1]), v[0]])
        else:
            tmp.append([str(k), v[0]])
    tmp_width = max(len(x[0]) for x in tmp)
    for v in tmp:
        print >>info_outfp, '    '+v[0].ljust(tmp_width+2,'.')+':', propRepr(v[1], num_sites)
    print >>info_outfp, "    **Total Filtered** :", propRepr(total_filtered, num_sites)
    print >>info_outfp, "Num sites kept :", propRepr(num_sites_kept, num_sites)



## misc functions ########################################################
# @TCC some of these may be general enough to move to separate modules

def sparseSlidingWindow(data, size, step, start=0, end=None, key=lambda x:x[0]):
    """generator which produces chunks of the input data in windows of size 'size', stepping 'step' between windows
       windows are defined in the half-open manner "[st,st+size)"; for ones-based ints, just used start=1
       key(element) {default: element[0]} must be an index/position
       DATA MUST BE IN ORDER BY KEY
       if 'start' is None, start of first window will be the index/position of first element of data
       Note: it will not produce windows past what is in the data even if a 'end' is provided
    """
    window_idx = 0
    data_it = iter(data) # ensure we are using an iterator
    buf = [data_it.next()] # start with the first element... buf should always contain at least 1 element
    # discard any data before start (if we have a start value)
    if( start is not None ):
        while( key(buf[-1]) < start ):
            last_x = key(buf[-1]) # just for checking data is in order
            buf = [data_it.next()]
            assert key(buf[-1]) >= last_x, "Input data must be in order/sorted by key"
            # having data_it.next() raise StopIteration (if runs out of data) is what we want to happen
    else: # start is None; so use the x of the first element as the start
        start = key(buf[-1])
    # the main loop
    while buf[-1] is not None : # end of buf being None signifies all data read
        # generate the next window range
        window_start = start+step*window_idx
        window_end = window_start + size
        # add data to end of buf until we hit one with x >= window_end
        while( buf[-1] is not None and key(buf[-1]) < window_end ):
            try:
                last_x = key(buf[-1]) # just for checking data is in order
                buf.append(data_it.next())
                assert key(buf[-1]) >= last_x, "Input data must be in order/sorted by key"
                if( end is not None and key(buf[-1]) > end ): # end condition (if set)
                    buf[-1] = None
            except StopIteration:
                buf.append(None) # signifies the end of the list
        # remove data from front of buf which is outside window
        while( key(buf[0]) < window_start ):
            del buf[0]
        # next loop (if there is one) will be next window
        window_idx += 1
        # yield (return this window and wait until we are asked for the next to continue the main loop)
        yield ((window_start, window_end), buf[0:-1])


def propRepr(num, total, show_total=True, places=2, num_width='total', total_width=None):
    """Return a string representing a proportion... for eg:
       if num = 10 and total = 100, return '10 / 100 (10.00%)'
    """
    # figure the widths
    if( total_width is None ):
        total_width = len(str(total))
    if( num_width is None ):
        num_width = len(str(num))
    elif( num_width == 'total' ):
        num_width = total_width
    # num (and possibly total)
    if( show_total ):
        s = '{0} / {1} '.format(str(num).rjust(num_width), str(total).rjust(total_width))
    else:
        s = '{0} '.format(str(num).rjust(num_width))
    # percent part
    if( num is None ):
        s += '(None)'
    elif( total == 0 ):
        s += '(INF)'
    else:
        if( places <= 0 ):
            fmt = ':d'
        else:
            fmt = ':{0}.{1}f'.format(places+3, places)
        fmt = '({'+fmt+'}%)'
        s += fmt.format(100*num/float(total))
    return s


#########################################################################
# Main loop hook... if run as script run main, else this is just a module
if __name__ == "__main__":
    sys.exit(Main(argv=None))
