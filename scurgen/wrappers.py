import subprocess
from toolshed import nopen
from toolshed.files import ProcessException
import numpy as np


def bigwig_summary(bigwig, chrom, start, stop, bins, func='mean'):
    """
    Returns a NumPy array containing summarized bins.

    Assumes that you have installed bigWigSummary from UCSC.

    :bigwig:
        Filename or URL of a bigWig format file

    :chrom, start, stop:
        Genomic coords to summarize

    :bins:
        Number of bins to summarize; the length of the returned NumPy array
        will be this long.

    :func:
        Function to apply for each bin; default is "max".  Can be one of (mean,
        min, max, std, coverage).

    """
    cmds = [
        'bigWigSummary',
        bigwig,
        chrom,
        str(int(start)),
        str(int(stop)),
        str(bins),
        '-type=%s' % func]
    y = []
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    results, stderr = p.communicate()
    for i in results.split('\t'):
        try:
            y.append(float(i))
        except ValueError:
            y.append(np.NaN)
    return np.array(y)


def bigbed_summary(bigwig, chrom, start, stop, bins, func='mean'):
    """
    Returns a NumPy array containing summarized bins.

    Assumes that you have installed bigBedSummary from UCSC.

    :bigwig:
        Filename or URL of a bigWig format file

    :chrom, start, stop:
        Genomic coords to summarize

    :bins:
        Number of bins to summarize; the length of the returned NumPy array
        will be this long.

    :func:
        Function to apply for each bin; default is "max".  Can be one of (mean,
        min, max, std, coverage).

    """
    cmds = [
        'bigBedSummary',
        bigwig,
        chrom,
        str(int(start)),
        str(int(stop)),
        str(bins),
        '-type=%s' % func]

    y = []
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    results, stderr = p.communicate()
    for i in results.split('\t'):
        try:
            y.append(float(i))
        except ValueError:
            y.append(np.NaN)
    return np.array(y)


if __name__ == "__main__":
    import time
    import os
    t00 = time.time()
    from matplotlib import pyplot as plt
    plt.rcParams['font.size'] = 10
    bigbed_url = ('http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/'
                  'integration_data_jan2011/byDataType/peaks/jan2011/'
                  'spp/optimal/hub/spp.optimal.wgEncodeHaibTfbsHepg2R'
                  'ad21V0416101AlnRep0_VS_wgEncodeHaibTfbsHepg2Contro'
                  'lV0416101AlnRep0.bb')
    bigwig_url = ('http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/'
                  'integration_data_jan2011/byDataType/rna_signal/jan'
                  '2011/hub/wgEncodeCshlLongRnaSeqAg04450CellPapPlusR'
                  'awSigRep1.bigWig')
    chrom = 'chr1'
    start = 1e8
    stop = 2e8
    bins = 128 * 128
    funcs = ['mean', 'std', 'coverage', 'min', 'max']
    fig = plt.figure()
    nax = len(funcs) + 1
    x = np.linspace(start, stop, bins)
    t = 0
    for i, func in enumerate(funcs):
        t0 = time.time()
        y1 = bigwig_summary(bigwig_url, chrom, start, stop, bins, func=func)
        elapsed = time.time() - t0
        t += elapsed
        print '%.2fs to download %s bigWig values for %s' \
            % (elapsed, bins, func)
        if i == 0:
            ax = fig.add_subplot(nax, 1, i + 1)
        else:
            ax = fig.add_subplot(nax, 1, i + 1, sharex=fig.axes[0])
        ax.plot(x, np.log1p(y1))
        ax.set_ylabel(func)
    ax = fig.add_subplot(nax, 1, i + 2, sharex=fig.axes[0])
    t0 = time.time()
    y2 = bigbed_summary(bigbed_url, chrom, start, stop, bins, func='max')
    elapsed = time.time() - t0
    t += elapsed
    print '%.2fs to download %s bigBed values for %s' % (elapsed, bins, func)
    ax.plot(x, y2)
    ax.set_ylabel('max peaks\nin bin')
    fig.subplots_adjust(hspace=0.4)
    ax.axis('tight')
    print '(%.2fs total download time)' % t
    print '(%.2fs total everything time)' % (time.time() - t00)
    plt.show()
