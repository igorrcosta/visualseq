#!/usr/bin/env python
# coding: utf-8
# visualseq.py

'''This program plots a sequence conservation graph of the comparison between three multiple alignment fasta files.
Run with -h to see the options.'''

__author__ = 'Igor Rodrigues da Costa'
__copyright__ = "Copyright 2012-2013, Igor Rodrigues da Costa"
__contact__ = 'igor.bioinfo@gmail.com'
__licence__ = '''
Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.'''

import argparse
from os.path import isfile
from matplotlib.pyplot import show, plot, title, ylim, xlim, savefig, rcParams, clf


ID_NUC = {
('T', 'T') : 1, ('A', 'A') : 1, ('G', 'G') : 1, ('C', 'C') : 1,
('T', 'A') : 0, ('T', 'G') : 0, ('T', 'C') : 0,
('A', 'T') : 0, ('A', 'G') : 0, ('A', 'C') : 0,
('G', 'T') : 0, ('G', 'A') : 0, ('G', 'C') : 0,
('C', 'T') : 0, ('C', 'A') : 0, ('C', 'G') : 0,
('N', 'T') : 0, ('N', 'A') : 0, ('N', 'G') : 0, ('N', 'C') : 0}

BLOSUM = {
('W', 'F') : 1, ('L', 'R') : -2, ('S', 'P') : -1, ('V', 'T') : 0, 
('Q', 'Q') : 5, ('N', 'A') : -2, ('Z', 'Y') : -2, ('W', 'R') : -3, 
('Q', 'A') : -1, ('S', 'D') : 0, ('H', 'H') : 8, ('S', 'H') : -1, 
('H', 'D') : -1, ('L', 'N') : -3, ('W', 'A') : -3, ('Y', 'M') : -1, 
('G', 'R') : -2, ('Y', 'I') : -1, ('Y', 'E') : -2, ('B', 'Y') : -3, 
('Y', 'A') : -2, ('V', 'D') : -3, ('B', 'S') : 0, ('Y', 'Y') : 7, 
('G', 'N') : 0, ('E', 'C') : -4, ('Y', 'Q') : -1, ('Z', 'Z') : 4, 
('V', 'A') : 0, ('C', 'C') : 9, ('M', 'R') : -1, ('V', 'E') : -2, 
('T', 'N') : 0, ('P', 'P') : 7, ('V', 'I') : 3, ('V', 'S') : -2, 
('Z', 'P') : -1, ('V', 'M') : 1, ('T', 'F') : -2, ('V', 'Q') : -2, 
('K', 'K') : 5, ('P', 'D') : -1, ('I', 'H') : -3, ('I', 'D') : -3, 
('T', 'R') : -1, ('P', 'L') : -3, ('K', 'G') : -2, ('M', 'N') : -2, 
('P', 'H') : -2, ('F', 'Q') : -3, ('Z', 'G') : -2, ('X', 'L') : -1, 
('T', 'M') : -1, ('Z', 'C') : -3, ('X', 'H') : -1, ('D', 'R') : -2, 
('B', 'W') : -4, ('X', 'D') : -1, ('Z', 'K') : 1, ('F', 'A') : -2, 
('Z', 'W') : -3, ('F', 'E') : -3, ('D', 'N') : 1, ('B', 'K') : 0, 
('X', 'X') : -1, ('F', 'I') : 0, ('B', 'G') : -1, ('X', 'T') : 0, 
('F', 'M') : 0, ('B', 'C') : -3, ('Z', 'I') : -3, ('Z', 'V') : -2, 
('S', 'S') : 4, ('L', 'Q') : -2, ('W', 'E') : -3, ('Q', 'R') : 1, 
('N', 'N') : 6, ('W', 'M') : -1, ('Q', 'C') : -3, ('W', 'I') : -3, 
('S', 'C') : -1, ('L', 'A') : -1, ('S', 'G') : 0, ('L', 'E') : -3, 
('W', 'Q') : -2, ('H', 'G') : -2, ('S', 'K') : 0, ('Q', 'N') : 0, 
('N', 'R') : 0, ('H', 'C') : -3, ('Y', 'N') : -2, ('G', 'Q') : -2, 
('Y', 'F') : 3, ('C', 'A') : 0, ('V', 'L') : 1, ('G', 'E') : -2, 
('G', 'A') : 0, ('K', 'R') : 2, ('E', 'D') : 2, ('Y', 'R') : -2, 
('M', 'Q') : 0, ('T', 'I') : -1, ('C', 'D') : -3, ('V', 'F') : -1, 
('T', 'A') : 0, ('T', 'P') : -1, ('B', 'P') : -2, ('T', 'E') : -1, 
('V', 'N') : -3, ('P', 'G') : -2, ('M', 'A') : -1, ('K', 'H') : -1, 
('V', 'R') : -3, ('P', 'C') : -3, ('M', 'E') : -2, ('K', 'L') : -2, 
('V', 'V') : 4, ('M', 'I') : 1, ('T', 'Q') : -1, ('I', 'G') : -4, 
('P', 'K') : -1, ('M', 'M') : 5, ('K', 'D') : -1, ('I', 'C') : -1, 
('Z', 'D') : 1, ('F', 'R') : -3, ('X', 'K') : -1, ('Q', 'D') : 0, 
('X', 'G') : -1, ('Z', 'L') : -3, ('X', 'C') : -2, ('Z', 'H') : 0, 
('B', 'L') : -4, ('B', 'H') : 0, ('F', 'F') : 6, ('X', 'W') : -2, 
('B', 'D') : 4, ('D', 'A') : -2, ('S', 'L') : -2, ('X', 'S') : 0, 
('F', 'N') : -3, ('S', 'R') : -1, ('W', 'D') : -4, ('V', 'Y') : -1, 
('W', 'L') : -2, ('H', 'R') : 0, ('W', 'H') : -2, ('H', 'N') : 1, 
('W', 'T') : -2, ('T', 'T') : 5, ('S', 'F') : -2, ('W', 'P') : -4, 
('L', 'D') : -4, ('B', 'I') : -3, ('L', 'H') : -3, ('S', 'N') : 1, 
('B', 'T') : -1, ('L', 'L') : 4, ('Y', 'K') : -2, ('E', 'Q') : 2, 
('Y', 'G') : -3, ('Z', 'S') : 0, ('Y', 'C') : -2, ('G', 'D') : -1, 
('B', 'V') : -3, ('E', 'A') : -1, ('Y', 'W') : 2, ('E', 'E') : 5, 
('Y', 'S') : -2, ('C', 'N') : -3, ('V', 'C') : -1, ('T', 'H') : -2, 
('P', 'R') : -2, ('V', 'G') : -3, ('T', 'L') : -1, ('V', 'K') : -2, 
('K', 'Q') : 1, ('R', 'A') : -1, ('I', 'R') : -3, ('T', 'D') : -1, 
('P', 'F') : -4, ('I', 'N') : -3, ('K', 'I') : -3, ('M', 'D') : -3, 
('V', 'W') : -3, ('W', 'W') : 11, ('M', 'H') : -2, ('P', 'N') : -2, 
('K', 'A') : -1, ('M', 'L') : 2, ('K', 'E') : 1, ('Z', 'E') : 4, 
('X', 'N') : -1, ('Z', 'A') : -1, ('Z', 'M') : -1, ('X', 'F') : -1, 
('K', 'C') : -3, ('B', 'Q') : 0, ('X', 'B') : -1, ('B', 'M') : -3, 
('F', 'C') : -2, ('Z', 'Q') : 3, ('X', 'Z') : -1, ('F', 'G') : -3, 
('B', 'E') : 1, ('X', 'V') : -1, ('F', 'K') : -3, ('B', 'A') : -2, 
('X', 'R') : -1, ('D', 'D') : 6, ('W', 'G') : -2, ('Z', 'F') : -3, 
('S', 'Q') : 0, ('W', 'C') : -2, ('W', 'K') : -3, ('H', 'Q') : 0, 
('L', 'C') : -1, ('W', 'N') : -4, ('S', 'A') : 1, ('L', 'G') : -4, 
('W', 'S') : -3, ('S', 'E') : 0, ('H', 'E') : 0, ('S', 'I') : -2, 
('H', 'A') : -2, ('S', 'M') : -1, ('Y', 'L') : -1, ('Y', 'H') : 2, 
('Y', 'D') : -3, ('E', 'R') : 0, ('X', 'P') : -2, ('G', 'G') : 6, 
('G', 'C') : -3, ('E', 'N') : 0, ('Y', 'T') : -2, ('Y', 'P') : -3, 
('T', 'K') : -1, ('A', 'A') : 4, ('P', 'Q') : -1, ('T', 'C') : -1, 
('V', 'H') : -3, ('T', 'G') : -2, ('I', 'Q') : -3, ('Z', 'T') : -1, 
('C', 'R') : -3, ('V', 'P') : -2, ('P', 'E') : -1, ('M', 'C') : -1, 
('K', 'N') : 0, ('I', 'I') : 4, ('P', 'A') : -1, ('M', 'G') : -3, 
('T', 'S') : 1, ('I', 'E') : -3, ('P', 'M') : -2, ('M', 'K') : -1, 
('I', 'A') : -1, ('P', 'I') : -3, ('R', 'R') : 5, ('X', 'M') : -1, 
('L', 'I') : 2, ('X', 'I') : -1, ('Z', 'B') : 1, ('X', 'E') : -1, 
('Z', 'N') : 0, ('X', 'A') : 0, ('B', 'R') : -1, ('B', 'N') : 3, 
('F', 'D') : -3, ('X', 'Y') : -1, ('Z', 'R') : 0, ('F', 'H') : -1, 
('B', 'F') : -3, ('F', 'L') : 0, ('X', 'Q') : -1, ('B', 'B') : 4
}

ID_MATRIX = {}
for x in BLOSUM.keys():
    if x[0] == x[1]:
        ID_MATRIX[x] = 1
    else:
        ID_MATRIX[x] = 0

def argument_parser(h = False, args = []):
    '''visualseq.py -i f1.fasta f2.fasta f3.fasta -o outfile.png -t Title -m Matrix -a Filter -x Xlim
    '''
    
    parser = argparse.ArgumentParser(description = 'Prints an image comparing 3 multiple alignment files (1 group for file). \
                                     All sequences must have the same lenght.',\
                                     argument_default = None)
    parser.add_argument('-i', '--infile', nargs = 3, metavar = ('f1.fasta', 'f2.fasta', 'f3.fasta'), type = str,\
                        help = '3 multiple alignment fasta files (not necessarily the same number \
                        of sequences in every file, but every sequence must have the same size).')
    parser.add_argument('-o', '--outfile', nargs = '?', metavar = 'outfile', type = str, default = 'outfile.png',\
                        help = 'file where the image will be saved. (defaut: %(default)s)')
    parser.add_argument('-t', '--title', nargs = '?', metavar = '\'Title\'', type = str, default = 'Title',\
                        help = 'Title of the image. (defaut: %(default)s)')
    parser.add_argument('-a', '--alpha', nargs = '?', metavar = 'alpha', type = float, default = 0.05,\
                        help = 'Value of the alpha parameter of the low-pass filter, lower means smoother. (defaut: %(default)s)')
    parser.add_argument('-m', '--matrix', nargs = '?', metavar = 'matrix', type = str, choices = ['BLOSUM', 'ID_MATRIX'], default = 'BLOSUM',\
                        help = 'Comparison matrix to be used. (options: %(choices)s, defaut: %(default)s)')
    parser.add_argument('-x', '--xlim', nargs = 2, metavar = ('start', 'end'), type = float, default = [0, 0],\
                        help = 'Limit values for the start and end of the x axis.'\
                        + 'If not specified, it is set to automatic. Both values must be numbers. (defaut: %(default)s)')
    parser.add_argument('-y', '--ylim', nargs = 2, metavar = ('start', 'end'), type = float, default = [0, 0],\
                        help = 'Limit values for the start and end of the y axis.'\
                        + 'If not specified, it is set to automatic. Both values must be numbers. (defaut: %(default)s)')
    parser.add_argument('-f', '--font_size', nargs = '?', metavar = 'size', type = int, default = 18,\
                        help = 'Font size of the axis labels and title. (defaut: %(default)s)')
    parser.add_argument('-s', '--size', nargs = 2, metavar = ('width', 'height'), type = int, default = [16, 8],\
                        help = 'Size of the figure, Width Height. (defaut: %(default)s)')
    parser.add_argument('-n', '--no_lowpass', action = 'store_true', help = 'Disable the low_pass filter.')

    if h:
        args = parser.parse_args(['-h'])
    else:
        try:
            if args:
                args = parser.parse_args(args).__dict__
            else:
                args = parser.parse_args().__dict__
        except:
            args = parser.parse_args(['-h'])
            raise
    return args

def fasta_parser(fasta_file):
    if not isfile(fasta_file):
        print fasta_file, 'is not a file.'
        return []
    seq_dict = {}
    seq = ''
    seq_id = ''
    with open(fasta_file, 'r') as f:
        for seq_line in f.readlines():
            if seq_line[0] == '>':
                if seq:
                    seq_dict[seq_id] = seq
                seq_id = seq_line[1:-1]
                seq = ''
            else:
                seq += seq_line[:-1]
        if seq_id:
            seq_dict[seq_id] = seq
        else:
            return []
    return seq_dict.values()

def comparador(lseq1, lseq2, matrix):
    try:
        size = len(lseq1[0])
    except:
        size = len(lseq2[0])
    mut = [0]*size
    pon = [0]*size
    for seq1 in lseq1:
        for seq2 in lseq2:
            for n, i in enumerate(zip(seq1, seq2)):
                if '-' not in i:
                    try:
                        bmt = matrix[(i[0], i[1])]
                    except KeyError:
                        bmt = matrix[(i[1], i[0])]
                    mut[n] += bmt
                    pon[n] += 1
    p = max(pon)
    comp = []
    for m, n in zip(mut, pon):
        if n != 0:
            i = float(m)/n
            comp.append(i)
        else:
##            comp.append(None)
            try:
                comp.append(comp[-1])
            except IndexError:
                comp.append(0)
    return comp, pon


def avg_lowpass(dlist, a = 0.05):
    alist = lowpass(dlist, a)
    dlist.reverse()
    alist_reverse = lowpass(dlist, a)
    alist_reverse.reverse()
    alist_avg = []
    last = None
    for n, m in zip(alist, alist_reverse):
        if n is not None and m is not None:
            avg = (float(n) + m)/2
            last = avg
        else:
            avg = last
        alist_avg.append(avg)
        #alist_avg.append((n,m))
    return alist_avg

def old_lowpass(dlist, a = 0.05):
    if len(dlist) == 0:
        return dlist
    alist = [0] * len(dlist)
    alist[0] = dlist[0]
    for p, i in enumerate(dlist[1:]):
        n = p + 1
        alist[n] = i * a + (1 - a) * alist[n-1]
    return alist

def lowpass(dlist, a = 0.05):
    if len(dlist) == 0:
        return []
    nones = 0
    start = 0
    alist = [0] * len(dlist)
    for n, i in enumerate(dlist):
        if i is not None:
            alist[n] = i
            start = n + 1
            break
        alist[n] = None
    if alist[n] is None:
        return alist
    for p, i in enumerate(dlist[start:]):
        n = p + start
        if i is None:
            nones += 1
            alist[n] = None
        elif nones > 0:
            last = alist[n - (nones + 1)]
            last = last * pow((1 - a), nones)
            alist[n] = i * a + last
            nones = 0
        else:
            alist[n] = i * a + (1 - a) * alist[n-1]
    return alist

def run(path, matrix, alpha = 0.05, lowpass = True):
    if matrix == 'BLOSUM':
        matrix = BLOSUM
    elif matrix == 'ID_MATRIX':
        matrix = ID_MATRIX
    rr = []
    listaseq1 = fasta_parser(path[0])
    listaseq2 = fasta_parser(path[1])
    listaseq3 = fasta_parser(path[2])
    if len(listaseq1) > 0 and len(listaseq2) > 0:
        y1, n1 = comparador(listaseq1, listaseq2, matrix)
    else:
        y1 = []
        n1 = []
    if len(listaseq2) > 0 and len(listaseq3) > 0:
        y2, n2 = comparador(listaseq2, listaseq3, matrix)
    else:
        y2 = []
        n2 = []
    if len(listaseq1) > 0 and len(listaseq3) > 0:
        y3, n3 = comparador(listaseq1, listaseq3, matrix)
    else:
        y3 = []
        n3 = []
    if lowpass:
        y1 = avg_lowpass(y1, alpha)
        y2 = avg_lowpass(y2, alpha)
        y3 = avg_lowpass(y3, alpha)
    x1 = [x + 1 for x in range(len(y1))]
    x2 = [x + 1 for x in range(len(y2))]
    x3 = [x + 1 for x in range(len(y3))]
    plot1 = [(px1, py1) for px1, py1 in zip(x1, y1)]
    plot2 = [(px2, py2) for px2, py2 in zip(x2, y2)]
    plot3 = [(px3, py3) for px3, py3 in zip(x3, y3)]
    return plot1, n1, plot2, n2, plot3, n3   

def ploter(plot1, cv1, plot2, cv2, plot3, cv3, args):
    rcParams['figure.figsize'] = args['size'] #16, 8
    floor = 1
    factor = 0.9
    gap = 0.4
    points = len(plot1)
    limit_cv = float(max(cv1))
    for n, c in zip(xrange(points - 1), cv1):
        if c != 0:
            intensity = factor * (floor + ((c * (1-floor))/limit_cv))
        else:
            intensity = gap
        plot([plot1[n][0], plot1[n+1][0]], [plot1[n][1], plot1[n+1][1]], color = (0, intensity, intensity))
    points = len(plot2)
    limit_cv = float(max(cv2))
    for n, c in zip(xrange(points - 1), cv2):
        if c != 0:
            intensity = factor * (floor + ((c * (1-floor))/limit_cv))
        else:
            intensity = gap
        plot([plot2[n][0], plot2[n+1][0]], [plot2[n][1], plot2[n+1][1]], color = (intensity, 0, intensity))
    points = len(plot3)
    limit_cv = float(max(cv3))
    for n, c in zip(xrange(points - 1), cv3):
        if c != 0:
            intensity = factor * (floor + ((c * (1-floor))/limit_cv))
        else:
            intensity = gap
        plot([plot3[n][0], plot3[n+1][0]], [plot3[n][1], plot3[n+1][1]], color = (intensity, intensity, 0))
    #plot(x1, y1, '#009999', x2, y2, '#990099', x3, y3, '#999900', linewidth= 1.3)
    title(args['title'])
    xlimit = args['xlim']
    ylimit = args['ylim']
    if (xlimit[0] != 0 and xlimit[1] != 0) and (type(xlimit[0]) == type(xlimit[1]) == float):
        xlim(xlimit)
    else:
        xlimit = (1, max([len(plot1), len(plot2), len(plot3)]))
        xlim(xlimit)
    if (ylimit[0] != 0 and ylimit[1] != 0) and (type(ylimit[0]) == type(ylimit[1]) == float):
        ylim(ylimit)
##    else:
##        ny = [y for y in y1 + y2 + y3 if y is not None]
##        ylimit = (min(ny) - 1, max(ny) + 1)
##        ylim(ylimit)
    rcParams.update({'font.size': args['font_size']})
    savefig(args['outfile'], bbox_inches=0)
    show()
    clf()
    
if  __name__ == "__main__":
    #test: python visualseq.py -i ./test/p.fas ./test/f.fas ./test/m.fas -o ./test/out_test.png -t 'Test'
##    1 = green, blue (1, 2)
##    2 = blue, red (2, 3)
##    3 = green, red (1, 3)
    args = argument_parser()
    if not args['infile'] or len(args['infile']) < 2:
        print 'too few fasta files (need 3)'
        argument_parser(h = True)
    plot1, n1, plot2, n2, plot3, n3 = run(args['infile'], args['matrix'], args['alpha'], not args['no_lowpass'])
    ploter(plot1, n1, plot2, n2, plot3, n3, args)
