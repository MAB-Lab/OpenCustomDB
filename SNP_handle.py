import re

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

def deletion(HGVS_C):
    if'_' in HGVS_C:
          rangepos = range(int(re.split('(\d+)',HGVS_C)[1]),int(re.split('(\d+)',HGVS_C)[3])+1) 
          nts = [x for x in HGVS_C.split('del')[1]]
          parsed_HGVS_C = []
          for p,n in zip(rangepos,nts):
              pHGVSC = {'effect':'del','nt':n,'pos':int(p)}
              parsed_HGVS_C.append(pHGVSC)
    else:
        pos = int(re.split('(\d+)',HGVS_C)[1])
        nt = re.split('del',HGVS_C)[1]
        pHGVSC = {'effect':'del','nt':nt,'pos':pos}
        parsed_HGVS_C = []
        parsed_HGVS_C.append(pHGVSC)
    return parsed_HGVS_C

def insertion(HGVS_C):
    pos = int(re.split('(\d+)',HGVS_C)[1])
    nt = re.split('ins',HGVS_C)[1]
    pHGVSC = {'effect':'ins','nt':nt,'pos':pos}
    parsed_HGVS_C = []
    parsed_HGVS_C.append(pHGVSC)
    
    return parsed_HGVS_C

def duplication(HGVS_C):
    pos = int(re.split('(\d+)',HGVS_C)[1])
    nt = re.split('dup',HGVS_C)[1]
    pHGVSC = {'effect':'ins','nt':nt,'pos':pos}
    parsed_HGVS_C = []
    parsed_HGVS_C.append(pHGVSC)
    return parsed_HGVS_C

def missense(HGVS_C):
    pos = int(re.split('(\d+)',HGVS_C)[1])
    nt = re.split('>',HGVS_C)[1]
    pHGVSC = {'effect':'replace','nt':nt,'pos':pos}
    parsed_HGVS_C = []
    parsed_HGVS_C.append(pHGVSC)
    return parsed_HGVS_C