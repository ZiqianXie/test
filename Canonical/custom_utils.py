#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pickle as pkl
from random import randint
from collections import defaultdict


# Due to memory limit (google colab now only has 12G Mem), the chromosomes are
# read one at a time. The lenght of each exon segments are prefetched so the
# replacement process go through the whole genome twice, once for extracting
# random segments, once for overwriting the exons using extracted random
# segments. Some Cython optimization may further improve the speed, but it is
# too tedious to do so.
class Chromosome:
    def __init__(self, dataset='canonical_dataset.txt'):
        self.chr_list = ['chr'+str(x) for x in list(range(1, 23)) + ['X', 'Y']]
        self.dataset = dataset
        self.current_ch = None
        self.current_cl = None
        self.download()
        self.chr_len()
        self.exon_prefetch(dataset)

    def download(self):
        url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/' +\
              'chromosomes/{}.fa.gz'
        for ch in self.chr_list:
            gz_exist = os.path.isfile('{}.fa.gz'.format(ch))
            fa_exist = os.path.isfile('{}.fa'.format(ch))
            if not (gz_exist or fa_exist):
                os.system('wget '+url.format(ch))
            if not fa_exist:
                os.system('gzip -d {}.fa.gz'.format(ch))
                os.system('rm *.fa.gz')

    def chr_len(self):
        # get and save the length of each chromosome
        fname = 'chr_len.pkl'
        if not os.path.isfile(fname):
            self.chr_len = {}
            for ch in reversed(self.chr_list):
                # chr1 would be the current_ch after the loop
                self.chr_len[ch] = len(self.read_chr(ch))
            self.chr_len['total_len'] = sum(self.chr_len.values())
            pkl.dump(self.chr_len, open(fname, 'wb'))
        else:
            self.chr_len = pkl.load(open(fname, 'rb'))

    def read_chr(self, ch):
        if ch == self.current_ch:
            return self.current_cl
        self.current_cl = None  # deallocate space first to avoid memory error
        with open('{}.fa'.format(ch), 'r') as f:
            self.current_ch = ch
            self.current_cl = list(f.read()[len(ch)+2:].replace('\n', ''))
            return self.current_cl

    def rand_seg(self):
        rand_seg_dict = defaultdict(list)
        for orig_ch, segs in self.exon_dict.items():
            for orig_start, seg_len in segs:
                start = randint(0, self.chr_len['total_len']-seg_len*24)
                for ch in self.chr_list:
                    len_chr = self.chr_len[ch]
                    if start <= len_chr - seg_len:
                        rand_seg_dict[ch].append((start, seg_len, orig_ch,
                                                 orig_start))
                        break
                    start -= len_chr - seg_len
                # assert start + seg_len <= len_chr
        return rand_seg_dict

    def exon_prefetch(self, dataset='canonical_dataset.txt'):
        # prefetch the position and length of all the exon segments in the
        # dataset file
        if os.path.isfile('exon_dict.pkl'):
            self.exon_dict = pkl.load(open('exon_dict.pkl', 'rb'))
            return
        self.exon_dict = defaultdict(list)
        with open(dataset, 'r') as f:
            for line in f.readlines():
                info = line.strip().split('\t')
                ch = info[2]
                gene_start, gene_end = info[4:6]
                intron_start, intron_end = info[6:]
                exon_start = [int(gene_start)] +\
                             [int(x) for x in
                              intron_end.strip(',').split(',')]
                exon_end = [int(x) for x in
                            intron_start.strip(',').split(',')] +\
                           [int(gene_end)]
                for i, start in enumerate(exon_start):
                    end = exon_end[i]
                    self.exon_dict[ch].append((start, end-start))
        pkl.dump(self.exon_dict, open('exon_dict.pkl', 'wb'))

    def replace_exon(self):
        # replace exons with random gene segments to make this model predict
        # cryptic splicing.
        # first pass extracts random segments
        rand_seg_dict = self.rand_seg()
        seg_dict = defaultdict(list)
        for ch in self.chr_list:
            v = rand_seg_dict[ch]
            print('first pass, {}'.format(ch))
            cl = self.read_chr(ch)
            for start, seg_len, orig_ch, orig_start in v:
                seg_dict[orig_ch].append((orig_start, seg_len,
                                          cl[start: start+seg_len]))
        # second pass replaces the exons with random segments
        with open('hg19.dup{}.fa'.format(self.file_id), 'w') as f:
            for ch in reversed(self.chr_list):
                # reverse the order of list so that the last loaded
                # chromosome can be reused to save some time.
                print('second pass, {}'.format(ch))
                cl = self.read_chr(ch)
                for start, seg_len, seg in seg_dict[ch]:
                    for i in range(seg_len):
                        cl[start+i] = seg[i]
                s = '>{}\n'.format(ch)
                f.write(s+''.join(cl)+'\n')
