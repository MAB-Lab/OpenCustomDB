import itertools as itt
import numpy as np
import gzip
import re
import collections
import Bio
import pyfaidx
import pickle
import sys
import argparse
import logging
from pyfaidx import FastaVariant, Fasta
from Bio import SeqIO
from collections import defaultdict
from collections import Counter
from collections import OrderedDict
from operator import itemgetter
#from pysnpeff.pysnpeff import SnpEff, SnpEffParser




def parse_arguments():
    """Read arguments from a command line."""
    parser = argparse.ArgumentParser(description='Create a personalised fasta from snpffOP and kallisto output.')
    
    #nom de l'experience
     parser.add_argument('experiment_name', type=str, help='name of your experiment')
    
    #Choix de l'annotation
    parser.add_argument('-a',  "--anotation",dest='annotation', type=int, choices=[0, 1, 2, 3, 4],
                    help="Choice of your build, 0 : Ensemble-OpenProt,/n   1: RefSeq-OpenProt,/n   2 : Ensemble,/n  3 : Refseq,/n                            4 : custom")
    
    #Choix de l'input
    parser.add_argument('-if',  "--inputformat",dest='iformat', type=int, choices=[0, 1, 2],
                    help="Choice of your input, 0 : VCF,  1 : protVCF, 2 : Tab.file")
    parser.add_argument('-ui', '--userInput', dest='userInput', type=str, help='Path your input file')
    
    #Si build custom
 
    parser.add_argument('-f', '--fasta', dest='fasta', type=str, help='Path to your protein.fasta.')
    parser.add_argument('-t', '--transcritfasta', dest='tfasta', type=str, help='Path to your transcrit.fasta.')
    parser.add_argument('-tsv', '--tsvprotein', dest='tsv', type=str, help='Path to your protein.tsv.')
    parser.add_argument('-b', '--OpenVarbuild', dest='build', type=str, help='Path to your OpenVarbuild')

    #Si expression transcrit
    parser.add_argument("-te", "--trxexpression", dest='trxexpression', type=str, help='path to tpm file')
    parser.add_argument("-ft", "--filtre_treshold",  dest='trxnumber', type=int, help='number of protein selected default 100000')
    parser.add_argument("-ex","--exception", dest='trxsave', type=str, help='exception file containing transcript accsesion')
   
                        
    args = parser.parse_args()
    
    return args



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


#Fonction principal call par l'utilisateur, 3 path en fonction du format de l'input

def OpenCMS(expname,iformat,annotation, fasta=None,tsv=None,tfasta=None,build=None,userInput=None,trxexpression=None,trxnumber=None,trxsave=None)
    if iformat == 0:
        fasta,tsv,tfasta,build = Get_annotation(annotation,fasta,tsv,tfasta,build)
        byopenvar(userInput,fasta,tsv,tfasta,build,expname,trxexpression=None,trxnumber=None,trxsave=None):   

    elif iformat ==1:
        fasta,tsv,tfasta,build = Get_annotation(annotation,fasta,tsv,tfasta,build)
        byprotvcf(userInput,fasta,tsv,tfasta,expname,trxexpression=None,trxnumber=None,trxsave=None): 
 
    elif iformat ==2:
        fasta,tsv,tfasta = Get_annotation(annotation,fasta,tsv,tfasta,build)
        bytabfile(userInput,fasta,tsv,tfasta,expname,trxexpression=None,trxnumber=None,trxsave=None):  
            
        return 'done'
    
#Call par la fonction principal, permet de définir l'annotation

def Get_annotation(annotation, fasta=None,tsv=None,tfasta=None,build=None ):
    if annotation:
        if annotation == 0
            openprot_ensembl_fasta = '/home/'
            openprot_ensembl_tsv = '/home/' 
            ensemble_transcrit_fasta = '/home/'
            openensembl_openvar_directory = '/home/'
        return(fasta,tsv,tfasta,build)
    
        elif annotation == 1
            openprot_refseq_fasta = '/home/'
            openprot_refseq_tsv = '/home/' 
            refseq_transcrit_fasta = '/home/'
            openrefseq_openvar_directory = '/home/'
        return(fasta,tsv,tfasta,build)
        
        elif annotation == 2
            ensembl_fasta = '/home/'
            ensembl_tsv = '/home/' 
            ensembl_transcrit_fasta = '/home/'
            ensembl_openvar_directory = '/home/'
        return(fasta,tsv,tfasta,build)
    
        elif annotation == 3
            refseq_fasta = '/home/'
            refseq_tsv = '/home/' 
            refseq_transcrit_fasta = '/home/'
            refseq_openvar_directory = '/home/'
        return(fasta,tsv,tfasta,build)
                
        elif annotation == 4
            if fasta:
                custom_fasta = fasta
            if tsv:
                custom_tsv = tsv 
            if ftrx:
                custom_transcrit_fasta = tfasta
            if build:
                custom_openvar_directory = build
            return(fasta,tsv,tfasta,build)

    else:
        return "error"
    
#Création de la BD custom si l'user input est un tab file

def bytabfile(userInput,fasta,tsv,tfasta,expname,trxexpression=None,trxnumber=None,trxsave=None):
    prot_syno = get_synonyms_prot(fasta)
    var_by_prot,transcrit_prot = parse_protvcf_file(userInput)
    start_codon= get_start_codon(tsv,tfasta)
    seqname_seq= get_all_mut_sequences(var_by_prot,transcrit_prot,start_codon)
    if trxexpression:
        trx_allprot= get_mutated_protbytranscrit(seqname_seq,transcrit_prot)
        trx_allprot= append_wt_prot_to_transcrit_ENS(tsv_open_prot,trx_allprot)
        AllProtInMyDB = get_100_prot(trxexpression,prot_syno,trx_allprot,trxsave=None)
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"), key_function=get_accession_seqio)
        DB_custom = assembling_headers_sequences(AllProtInMyDB,seqname_seq,prot_syno,fasta_dict)
        write_Fasta_DB(DB_custom,path_to_fasta)
    else:
        m_wtprot=get_m_wtprot(seqname_seq)
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"), key_function=get_accession_seqio)
        DB_custom = assembling_headers_sequences(m_wtprot,seqname_seq,prot_syno,fasta_dict)
            
    
#Création de la BD custom si l'user input est un vcf avec information protéique
    
def byprotvcf(userInput,fasta,tsv,tfasta,expname,trxexpression=None,trxnumber=None,trxsave=None):
    from pysnpeff.pysnpeff import SnpEff, SnpEffParser
    parser = SnpEffParser(userInput)
    prot_syno = get_synonyms_prot(fasta)
    var_by_prot,transcrit_prot = parse_protvcf_file(userInput)
    start_codon= get_start_codon(tsv,tfasta)
    seqname_seq= get_all_mut_sequences(var_by_prot,transcrit_prot,start_codon)
    if trxexpression:
        trx_allprot= get_mutated_protbytranscrit(seqname_seq,transcrit_prot)
        trx_allprot= append_wt_prot_to_transcrit_ENS(tsv_open_prot,trx_allprot)
        AllProtInMyDB = get_100_prot(trxexpression,prot_syno,trxnumber=None,trxsave=None)
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"), key_function=get_accession_seqio)
        DB_custom = assembling_headers_sequences(AllProtInMyDB,seqname_seq,prot_syno,fasta_dict)
        write_Fasta_DB(DB_custom,path_to_fasta)
    else:
        m_wtprot=get_m_wtprot(seqname_seq)
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"), key_function=get_accession_seqio)
        DB_custom = assembling_headers_sequences(m_wtprot,seqname_seq,prot_syno,fasta_dict)
         
    return 'Fasta available at '+path_to_fasta

#Création de la BD custom si l'user input est un vcf, passe par openvar
    
def byopenvar(userInput,fasta,tsv,tfasta,build,expname,trxexpression=None,trxnumber=None,trxsave=None):
        
    from pysnpeff.pysnpeff import SnpEff, SnpEffParser
    s = SnpEff(
    expname, 
    userVCF,
    '/home/',
    snpeff_dir = build
)
    s.split_vcf()
    s.run()
    parser = SnpEffParser('/home/')
    parsedsnpf = parser.parse_annOnePerLine()
    get_protvcf_file(parsenpff,expname)
    protvariantfile = ('expname.tab')
    prot_syno = get_synonyms_prot(fasta)
    var_by_prot,transcrit_prot = parse_protvcf_file(protvariantfile)
    start_codon= get_start_codon(tsv,tfasta)
    seqname_seq= get_all_mut_sequences(var_by_prot,transcrit_prot,start_codon)
    if trxexpression:
        trx_allprot= get_mutated_protbytranscrit(seqname_seq,transcrit_prot)
        trx_allprot= append_wt_prot_to_transcrit_ENS(tsv_open_prot,trx_allprot)
        AllProtInMyDB = get_100_prot(trxexpression,prot_syno,trx_allprot,trxnumber=None,trxsave=None)
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"), key_function=get_accession_seqio)
        DB_custom = assembling_headers_sequences(AllProtInMyDB,seqname_seq,prot_syno,fasta_dict)
        write_Fasta_DB(DB_custom,path_to_fasta)
    else:
        m_wtprot=get_m_wtprot(seqname_seq)
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"), key_function=get_accession_seqio)
        DB_custom = assembling_headers_sequences(m_wtprot,seqname_seq,prot_syno,fasta_dict)
        return 'Fasta available at '+path_to_fasta   

#utilise le tsv pour connaitre le codon start de chaque prot
def get_start_codon (tsv,transcrit_fasta_v_29):
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
            
    Ens_transcripts = Fasta(transcrit_fasta_v_29,key_function = lambda x: x.split('.')[0],duplicate_action="longest")        
        
    for x in transcrit_prot_ref_start:
        seq = Ens_transcripts[x['transcrit']][int(x['start'])-1:]
        prot_tx_start[str(x['prot'].split('.',1)[0]+'^'+x['transcrit'])] = seq
    
    return(prot_tx_start)    
    
#crée un tab.file à partir du parser de snpf
def get_protvcf_file(parsenpff,expname):
    parsenpff =sorted(parsenpff, key=lambda x: x['ANN[*].FEATUREID'])
    with open('expname.tab', 'w') as fwrite:
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
    return()
#lis le tab.file
def parse_protvcf_file(protvariantfile):
    var_by_prot = defaultdict(list)
    transcrit_prot = dict()
    with open(protvariantfile, 'r') as f:
        for n,l in enumerate(f):
            if n<1:continue
            prot = l.split("\t")[0]
            trx = l.split("\t")[1]
            HGVS_P = l.split("\t")[2]
            HGVS_C = l.split("\t")[3]
            if {'HGVS_P':HGVS_P,'HGVS_C':HGVS_C,} not in var_by_prot[prot]:
                var_by_prot[prot].append({'HGVS_P':HGVS_P,'HGVS_C':HGVS_C,})
                transcrit_prot[prot] = trx
            
    return var_by_prot,transcrit_prot

#######MUTTION RELATED#######

def Parse_HGVS_C(HGVS_C):
    if'>' in HGVS_C:
        parsed_HGVS_C = missense(HGVS_C)
        return parsed_HGVS_C
    if'dup'in HGVS_C:
        parsed_HGVS_C = duplication(HGVS_C)
        return parsed_HGVS_C
    if 'ins' in HGVS_C and 'del' in HGVS_C and '+' not in HGVS_C:
        HGVS_C_del=re.split('ins',HGVS_C)[0]
        parsed_HGVS_C = deletion(HGVS_C_del)
        HGVS_C_ins = re.split('del',HGVS_C)[0]+'ins'+re.split('ins',HGVS_C)[1]
        ins_parsed_HGVS_c = insertion(HGVS_C_ins)
        parsed_HGVS_C.extend(ins_parsed_HGVS_c)
        return parsed_HGVS_C
    if'del'in HGVS_C:
        parsed_HGVS_C = deletion(HGVS_C)
        return parsed_HGVS_C    
    if'ins'in HGVS_C:
        parsed_HGVS_C = insertion(HGVS_C)
        return parsed_HGVS_C
    
    return parsed_HGVS_C
    
    
def get_all_mut_sequences(var_by_prot,transcrit_prot,
                         ):
    seqname_seq=dict()
    for acc,svar in var_by_prot.items():
        regroupement_HGVS_C = list()
        regroupement_HGVS_P = list()
        for var in svar:
            parsed_HGVS_C = Parse_HGVS_C(var['HGVS_C'])
            parsed_HGVS_P = var['HGVS_P'].split('.',1)[1]
            regroupement_HGVS_C.extend(parsed_HGVS_C)
            regroupement_HGVS_P.append(parsed_HGVS_P)
        sorted_HGVS_C = sort_sequences(regroupement_HGVS_C)
        seqname = acc+'@'+''.join(regroupement_HGVS_P)
        mutated_sequence = modify_transcript_sequence(sorted_HGVS_C,acc,transcrit_prot,start_codon)
        if mutated_sequence == 'synonymous_variant' or translate(mutated_sequence) == 'start_lost' or len(mutated_sequence)<7:continue
        if translate(mutated_sequence) in seqname_seq.values():continue
        seqname_seq[seqname]=translate(mutated_sequence)
    return seqname_seq

def modify_transcript_sequence(HGVS_C_snp_sorted,protacc,transcrit_prot,prot_tx_start):
    trx = str(transcrit_prot[protacc])
    acc = protacc.split('.',1)[0]
    mseq = [x for x in str(prot_tx_start[acc+'^'+trx])]
    for variant in HGVS_C_snp_sorted:
        if variant['effect'] == 'del':
            if mseq[int(variant['pos'])-1]:
                mseq[int(variant['pos'])-1] = ""
        if variant['effect'] == 'ins':
            mseq.insert(int(variant['pos']), variant['nt'])
        if variant['effect'] == 'replace':
            mseq[int(variant['pos'])-1] = variant['nt']
    mseq[:] = [x for x in mseq if x !='']
    if translate(''.join(mseq)) == translate(''.join([x for x in str(prot_tx_start[acc+'^'+trx])])):
        return 'synonymous_variant',acc
    if len(translate(''.join(mseq))) < 8:
        return 'toosmall',acc
    mseq = "".join(mseq)
    return mseq

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

####séléction proteins related ######
def append_wt_prot_to_transcrit_ENS (tsv_open_prot,trx_allprot):
    with open(tsv_open_prot,'r')as f:
        transcrit_prot_ref_start = list()
        for n,l in enumerate(f):
            if n<2:continue
            ls = l.split('\t')
            prot = ls[0]
            gene = ls[7]
            start = ls[15]
            transcrit = ls[12].split('.',1)[0]
            if transcrit[:4] == 'ENST':
                if prot[:3] == 'ENS' or prot[:3] == 'IP_' or prot[:3] == 'II_':
                    trx_allprot[transcrit].append(prot.split('.',1)[0])
    
    return(trx_allprot)

def fastasynonymes(fasta_prot):
    prot_syno = get_synonyms_prot(fasta_prot)
    return prot_syno

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
    return(prot_syno)

def get_100_prot(trx_expression,prot_syno,trx_allprot,trxnumber=None,trxsave=None):
    trx_expression_sorted = get_trxs_by_tpm_from_kallisto(trx_expression)
    AllProtInMyDB = set()
    if trxnumber:
        treshold = trxnumber
        else
            treshold = 100000
    counter = 0
    for x in trx_expression_sorted:
        for y in trx_allprot[x]:    
            if len(AllProtInMyDB)<treshold:
                if y in prot_syno:
                    if prot_syno[y] not in AllProtInMyDB:
                        AllProtInMyDB.add(prot_syno[y])
                else:
                    AllProtInMyDB.add(y)  
    if trxsave:
        with open(trxsave, 'r') as f:
            for n,l in enumerate(f):
                for y in trx_allprot[l]:
                    if y in prot_syno:
                        if prot_syno[y] not in AllProtInMyDB:
                        AllProtInMyDB.add(prot_syno[y])
                    else:
                        AllProtInMyDB.add(y)                                      
    return(AllProtInMyDB)

def get_trxs_by_tpm_from_kallisto (kallisto):
    trxp_tpms = {}
    with open(kallisto, 'r') as f:
        for n,l in enumerate(f):
            ls = l.strip().split('\t')
            if n==0:
                keys = ls
                continue
            line = dict(zip(keys, ls))
            trxp = line['target_id'].split('.')[0]
            tpm = float(line['tpm'])
            if tpm >0:
                trxp_tpms[trxp] = tpm
    trxp_tpms_sorted_pour_tresh = OrderedDict(sorted(trxp_tpms.items(), key=itemgetter(1), reverse=True))
    return (trxp_tpms_sorted_pour_tresh)

def get_mutated_protbytranscrit(seqname_seq,transcrit_prot):
    trx_allprot = defaultdict(list)
    for prot_mut in seqname_seq:
        for prot,transcrit in transcrit_prot.items():
            if prot in prot_mut:
                trx_allprot[transcrit].append(prot_mut)
    return trx_allprot

def get_m_wtprot(seqname_seq):
    m_wtprot = list()
    for prot_mut in seqname_seq:
        m_wtprot.append(prot_mut.split('@')[0]
        m_wtprot.append(prot_mut.split('@')[0]
    return m_wtprot
                    
##### Création de la DB ######
def assembling_headers_sequences(AllProtInMyDB,Msequence,prot_syno,fasta_dict):
    DB_custom = dict()
    for acc in AllProtInMyDB:
        regular_acc = acc.split('@')[0].split('.')[0]

        if acc in Msequence:
            if regular_acc in prot_syno:
                header = acc+'|'+str(fasta_dict[prot_syno[regular_acc]].description).split('|')[1]
                DB_custom[header] = Msequence[acc]
            else:
                header = acc+'|'+str(fasta_dict[regular_acc].description).split('|')[1]
                DB_custom[header] = Msequence[acc]
        else:
            if acc in prot_syno:
                header = acc+'|'+str(fasta_dict[prot_syno[regular_acc]].description).split('|')[1]
                DB_custom[header] = str(fasta_dict[prot_syno[regular_acc]].seq)
            else:
                header = acc+'|'+str(fasta_dict[regular_acc].description).split('|')[1]
                DB_custom[header] = str(fasta_dict[regular_acc].seq)
    
    return DB_custom

fields = ['OS', 'GN', 'TA', 'PA']
def parse_fasta_header(h):
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

def truncate(seq, length=80):
    return '\n'.join([seq[n:n+length] for n in range(0, len(seq), length)])

def write_Fasta_DB(DB_custom,path):
    with open (path, 'w') as f:
        for acc, prot_seq in DB_custom.items():
            f.write('>'+acc+'\n')
            f.write(truncate(prot_seq+'\n'))
         


if __name__ == '__OpenCMS__':
    args = parse_arguments()
    main(args.snpffdir,args.snpffsave,args.path,args.tsv,args.fasta_transcrit,args.kallisto,args.fasta_prot,args.out_dir)
