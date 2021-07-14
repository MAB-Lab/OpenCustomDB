import itertools as itt
import numpy as np
import gzip
import re
import collections
import Bio
import pyfaidx
import pickle
import sys
from pyfaidx import FastaVariant, Fasta
from Bio import SeqIO
from collections import defaultdict
from collections import Counter
from collections import OrderedDict
from operator import itemgetter
from OpenVar.openvar import SeqStudy, OpenVar, OPVReport

AAcode = {
    'Gly':'G','Glu':'E','Asp':'D','Ala':'A','Val':'V',
    'Arg':'R','Ser':'S','Lys':'K','Asn':'N','Thr':'T',
    'Met':'M','Ile':'I','Gln':'Q','His':'H','Pro':'P',
    'Leu':'L','Trp':'W','Cys':'C','Tyr':'Y','Ser':'S','Phe':'F','*':'*', 'Ter':'*'}

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

fields = ['OS', 'GN', 'TA', 'PA']

def get_synonyms_prot(fasta_prot):
    with open(fasta_prot,'r')as f:
        prot_syno = dict()
        for n,l in enumerate(f):
            if '>' in l:
                prot = l.split('|',1)[0].split('.',1)[0].strip('>')
                parsed = parse_fasta_header(l)
                syn = parsed['PA'].split(',')
                for x in syn:
                    prot_syno[x.split('.',1)[0]]=prot
    return prot_syno
    
def parse_fasta_header(h):
    fields = ['OS', 'GN', 'TA', 'PA']
    h = h.split()
    acc = h[0].split('|')[0]
    res = {}
    for f in h[1:]:
        for field in fields:
            if f[:2] == field:
                res[field] = f[3:]
    return res

def get_accession_seqio(record):

    parts = record.id.split("|")[0].split('.')
    return parts[0]

def get_fasta_dict(fasta_prot):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_prot, "fasta"), key_function=get_accession_seqio)
    return fasta_dict

def translate(mseq, frame=0):
    mseq = mseq[frame:]
    codons = [mseq[n:n+3] for n in range(0,len(mseq),3)]
    translation = ''
    if codons[0] != 'ATG':
        return 'start_lost'
    for codon in codons:
        if len(codon)<3:
            return translation
        elif gencode[codon] == '_':
            return translation
        else:
            translation += gencode[codon]
    return translation

def sort_sequences(hgvsc_snp):
    ordre = {'replace':0,'del':1,'ins':2,'dup':2}
    hgvsc_snp_sorted =sorted(hgvsc_snp, key=lambda x: (ordre[x['effect']],-x['pos']))
    
    return hgvsc_snp_sorted

def remove_duplicata_from_db(DB_custom):
    duplicata = list()
    seek_duplicate = Counter(DB_custom.values())
    for x,y in seek_duplicate.items():
        if y>1:
            if '@' in x:
                duplicata.append(x)
    for x in duplicata:
        del DB_custom[x]
        
    return DB_custom

def remove_fakevariant(seqname_seq,fasta_dict,prot_syno):

    duplicata=list()
    for acc in fasta_dict: 
        if acc in prot_syno:
            if str(fasta_dict[prot_syno[acc]].seq) in seqname_seq.values():
                duplicata.append([k for k,v in seqname_seq.items() if v == str(fasta_dict[prot_syno[acc]].seq)])
        else:
            if str(fasta_dict[acc].seq) in seqname_seq.values():
                duplicata.append([k for k,v in seqname_seq.items() if v == str(fasta_dict[acc].seq)])
    flat_duplicata = [seq for sublist in duplicata for seq in sublist]
    for x in flat_duplicata:
        del seqname_seq[x]
    return seqname_seq

def get_start_codon (tsv,tfasta):

    prot_tx_start = dict()
    with open(tsv,'r')as f:
        transcrit_prot_ref_start = list()
        Wtsequence = dict()
        Prot_gene = dict()
        for n,l in enumerate(f):
            if n<2:continue
            ls = l.split('\t')
            prot = ls[0]
            gene = ls[7]
            start = ls[15]
            transcrit = ls[12].split('.',1)[0]
            if transcrit[:4] == 'ENST':
                if prot[:3] == 'ENS' or prot[:3] == 'IP_' or prot[:3] == 'II_':
                    info_start = {'transcrit':transcrit,'prot':prot,'start':start}
                    transcrit_prot_ref_start.append(info_start)
            
    Ens_transcripts = Fasta(tfasta,key_function = lambda x: x.split('.')[0],duplicate_action="longest")        
        
    for x in transcrit_prot_ref_start:
        seq = Ens_transcripts[x['transcrit']][int(x['start'])-1:]
        prot_tx_start[str(x['prot'].split('.',1)[0]+'^'+x['transcrit'])] = seq
    
    return prot_tx_start  
    
def OpenVar_analysis(vcf_path, expname, specie='human'):
    filename = vcf_path.split('/')[-1]
    path = vcf_path.replace(filename,'')
    vcf = SeqStudy(
    data_dir   = path, # .../uploads
    file_name  = filename, # guid.vcf
    results_dir = filename.replace('.vcf','')+'_result', # .../results/guid/
    study_name = expname, # user input
    specie = specie, # user input
    genome_version = 'hg38', # user input
    annotation  = 'OP_Ensembl' # user input
    )

    opv = OpenVar(
    snpeff_path = '/home/xroucou_group/echange_de_fichiers/snpEff/',
    vcf         = vcf,
    )
    opv.run_snpeff_parallel_pipe()
    opvr = OPVReport(opv)
    out = list()
    for f in opvr.annOnePerLine_files:
        out.extend(list(opvr.parse_annOnePerLine(f, as_dict=True)))
    return out    

def parse_protvcf_file(protvariantfile): #get_variant_by_prot

    var_by_prot = defaultdict(list)
    transcrit_prot = dict()
    with open(protvariantfile, 'r') as f:
        for n,l in enumerate(f):
            if n<1:continue
            ls = l.split("\t")
            prot = ls[0]
            trx = ls[1]
            HGVS_P = ls[2]
            HGVS_C = ls[3]
            Error = ls[4]
            if prot.count('_')>1:continue
            if Error!='\n':continue
            elif {'HGVS_P':HGVS_P,'HGVS_C':HGVS_C,} not in var_by_prot[prot]:
                var_by_prot[prot].append({'HGVS_P':HGVS_P,'HGVS_C':HGVS_C,})
                transcrit_prot[prot] = trx
            
    return var_by_prot,transcrit_prot

def fastasynonymes(fasta_prot):
    prot_syno = get_synonyms_prot(fasta_prot)
    return prot_syno

def get_protvcf_file(parsenpff, expname, vcf_path):
    parsenpff =sorted(parsenpff, key=lambda x: x['ANN[*].FEATUREID'])
    filename = vcf_path.split('/')[-1]
    path = filename.replace('.vcf','')+'_result/'+expname+'_prot.tab'
    
    with open(path, 'w') as fwrite:
        fwrite.write('Prot'+'\t'+'Transcrit'+'\t'+'HGVS_P'+'\t'+'HGVS_C'+'\t'+'Potential_Error'+'\n')
        for n in parsenpff:
            if n['ANN[*].HGVS_P']:
                if '@' in n['ANN[*].FEATUREID']:
                    prot = n['ANN[*].FEATUREID'].split('@',1)[1]
                    trans = n['ANN[*].FEATUREID'].split('@',1)[0]
                    fwrite.write(prot+'\t'+trans+'\t'+n['ANN[*].HGVS_P']+'\t'+n['ANN[*].HGVS_C']+'\t'+n['ANN[*].ERRORS']+'\n')
                elif '^' in n['ANN[*].FEATUREID']:
                    prot = n['ANN[*].FEATUREID'].split('^',1)[1]
                    trans = n['ANN[*].FEATUREID'].split('^',1)[0]
                    fwrite.write(prot+'\t'+trans+'\t'+n['ANN[*].HGVS_P']+'\t'+n['ANN[*].HGVS_C']+'\t'+n['ANN[*].ERRORS']+'\n')
                else:continue
                    
    return path

def truncate(seq, length=80):
    return '\n'.join([seq[n:n+length] for n in range(0, len(seq), length)])