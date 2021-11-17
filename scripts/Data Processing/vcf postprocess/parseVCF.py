#!/usr/bin/env python3

import sys

output_opt=''
output_opt = sys.argv[2]
filename = sys.argv[1]

ignore = False
class Variant:
    def __init__(self, gene, variant,type,impact):
        self.gene = gene
        self.variant = variant
        self.variant_type =type
        self.impact=impact
        self.info_per_sample=dict()
        self.chrom = ''
        self.start = ''
        self.annotation = []
    def addCoordinate(self, chrom, start):
        self.chrom = chrom
        self.start = start
    def addInfoPerSample(self, sample_name, info):
        self.info_per_sample[sample_name]=info
    def printInfoPerSample(self, samples, option='af'):
        if option=='all':
            s = ''
            for i in range(len(samples)-1):
                #s = s+"{:.2f}".format(self.info_per_sample[samples[i]])+'\t'
                s=s+self.info_per_sample[samples[i]]+'\t'
            s=s+self.info_per_sample[samples[-1]]
            return s
        elif option=='allsep':
            s = ''
            for i in range(len(samples)-1):
                t='\t'.join(self.info_per_sample[samples[i]].split(':'))
                s=s+t+'\t'
            s=s+'\t'.join(self.info_per_sample[samples[-1]].split(':'))
            return s
        else:
            s = ''
            for i in range(len(samples)-1):
                af = self.info_per_sample[samples[i]].split(':')[-1]
                s=s+af+'\t'
            s=s+self.info_per_sample[samples[-1]].split(':')[-1]
            return s
    def addAnnotation(self, ann):
        self.annotation.append(ann)

samples = []
variants = []
cur_position=''
pre_position=''
for line in open(filename, 'r'):
    line=line.strip()
    if line[0:2]!= '##':
        a = line.split('\t')
        if a[0]=='#CHROM':
            samples = a[9:]
        else:
            cur_position=a[0]+':'+a[1]
            if pre_position=='' or (pre_position!='' and cur_position!= pre_position):
                # parse ANN
                info = a[7].split(';')
                ann=[]
                ann_first=[]
                for c in info:
                    if c[0:3]== 'ANN':
                        ann = c[4:].split(',')
                        ann_first=ann[0].split('|')
                v=Variant(ann_first[3],ann_first[9],ann_first[1],ann_first[2])
                v.addCoordinate(a[0], a[1])
                for j in ann:
                    v.addAnnotation(j)
                # parse Format
                Format = a[8].split(':')
                DP_index = Format.index("DP")
                AO_index = Format.index("AO")
                RO_index = Format.index("RO")
                for i in range(len(samples)):
                    if '.' in a[9+i:]:
                        ignore=True
                        break
                    else:
                        detail = a[9+i].split(":")
                        if detail[AO_index].isnumeric() and detail[RO_index].isnumeric():
                            if int(detail[AO_index])+int(detail[RO_index]) > 0:
                                AF = float(int(detail[AO_index])/(int(detail[AO_index])+int(detail[RO_index])))
                                information = detail[DP_index]+':'+detail[RO_index]+':'+detail[AO_index]+':'+"{:.2f}".format(AF)
                            else:
                                information = detail[DP_index]+':'+detail[RO_index]+':'+detail[AO_index]+':'+'NA'
                        else:
                            information = detail[DP_index]+':'+detail[RO_index]+':'+detail[AO_index]+':'+'NA'
                        v.addInfoPerSample(samples[i],information)
                if not ignore:
                    variants.append(v)
                ignore=False
            pre_position=cur_position

if output_opt=='allsep':
    samples_surfix=[]
    for t in samples:
        samples_surfix=samples_surfix+[t+'_DP',t+'_RO',t+'_AO',t+'_AF']
    print('Chrom'+'\t'+'Location'+'\t'+'Gene'+'\t'+'Variant'+'\t'+'Type'+'\t'+'Impact'+'\t'+'\t'.join(samples_surfix))
else:
    print('Chrom'+'\t'+'Location'+'\t'+'Gene'+'\t'+'Variant'+'\t'+'Type'+'\t'+'Impact'+'\t'+'\t'.join(samples))

# You can choose to print only AF (option=af) value or DP:RO:AO:AF (option=all)
for v in variants:
    print(v.chrom+'\t'+v.start+'\t'+v.gene+'\t'+v.variant+'\t'+v.variant_type+'\t'+v.impact+'\t'+v.printInfoPerSample(samples,option=output_opt))

