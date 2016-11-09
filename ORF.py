# -*- coding: utf-8 -*-
#Problem 4.6.1.1
import myBio as bio
import myBioDataSet as biodata

#Recupere la table standard ou celle du mycoplasme
def getGeneticCode(NCBI_ID):
    
    if NCBI_ID==1:
        bases = ['T', 'C', 'A', 'G']
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids_one = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, amino_acids_one))
        return codon_table
    
    if NCBI_ID==4:
        bases = ['T', 'C', 'A', 'G']
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids_one = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, amino_acids_one))
        return codon_table

#Bien qu'existant dans myBio j'ai du modifier ces fonctions pour fonctionner avec une string et non un dictionnaire
def translate(seq,codonTable=1):
    codonT=getGeneticCode(codonTable)
    seq_word=[]
    seq_translated=''
    for i in range(0,len(seq),3):
        seq_word.append(seq[i:i+3])
    for i in seq_word:
        if len(i)==3:
            seq_translated += codonT[i]
    return seq_translated

def revcomp(seq):
    comp_result = ''
    co_bases = {'A' : 'T', 'G' : 'C', 'C' : 'G', 'T' : 'A', 'R' : 'Y', 'Y' : 'R'}
    for letter in seq:
        if letter in co_bases.keys():
            comp_result += co_bases[letter]
        else:
            comp_result += letter
    rev_result=comp_result[::-1]
    return rev_result



def findORF(seq, threshold,codeTable):
    nseq=seq['data']
    POS_D,POS_F,a=[],[],0
    a=0
    for i in range(0,len(nseq)-len(nseq)%3,3):
        if nseq[i:i+3]=='ATG' and a==0:
            POS_D.append(i)
            a=1
        if (nseq[i:i+3]=='TAA' or nseq[i:i+3]=='TAG' or nseq[i:i+3]=='TGA') and a==1:
            POS_F.append(i)
            a=0
    if a==1:
        POS_D.pop()
        a=0
    for i in range(1,len(nseq)-(len(nseq)-1)%3,3):
        if nseq[i:i+3]=='ATG' and a==0:
            POS_D.append(i)
            a=1
        if (nseq[i:i+3]=='TAA' or nseq[i:i+3]=='TAG' or nseq[i:i+3]=='TGA') and a==1:
            POS_F.append(i)
            a=0
    if a==1:
        POS_D.pop()
        a=0
    for i in range(2,len(nseq)-(len(nseq)-2)%3,3):
        if nseq[i:i+3]=='ATG' and a==0:
            POS_D.append(i)
            a=1
        if (nseq[i:i+3]=='TAA' or nseq[i:i+3]=='TAG' or nseq[i:i+3]=='TGA') and a==1:
            POS_F.append(i)
            a=0
    if a==1:
        POS_D.pop()
        a=0
    rev_nseq=revcomp(nseq)
    for i in range(0,len(rev_nseq)-len(nseq)%3,3):
        if rev_nseq[i:i+3]=='ATG' and a==0:
            POS_D.append(-(len(rev_nseq)-i))
            a=1
        if (rev_nseq[i:i+3]=='TAA' or rev_nseq[i:i+3]=='TAG' or rev_nseq[i:i+3]=='TGA') and a==1:
            POS_F.append(-(len(rev_nseq)-i))
            a=0
    if a==1:
        POS_D.pop()
        a=0
    for i in range(1,len(rev_nseq)-(len(nseq)-1)%3,3):
        if rev_nseq[i:i+3]=='ATG' and a==0:
            POS_D.append(-(len(rev_nseq)-i-1))
            a=1
        if (rev_nseq[i:i+3]=='TAA' or rev_nseq[i:i+3]=='TAG' or rev_nseq[i:i+3]=='TGA') and a==1:
            POS_F.append(-(len(rev_nseq)-i-1))
            a=0
    if a==1:
        POS_D.pop()
        a=0
    for i in range(2,len(rev_nseq)-(len(nseq)-2)%3,3):
        if rev_nseq[i:i+3]=='ATG' and a==0:
            POS_D.append(-(len(rev_nseq)-i-2))
            a=1
        if (rev_nseq[i:i+3]=='TAA' or rev_nseq[i:i+3]=='TAG' or rev_nseq[i:i+3]=='TGA') and a==1:
            POS_F.append(-(len(rev_nseq)-i-2))
            a=0
    if a==1:
        POS_D.pop()
        a=0

    newD=[]
    newF=[]

    for i in range(len(POS_D)-1):
        if POS_F[i]-POS_D[i] > threshold:
            newD.append(POS_D[i])
            newF.append(POS_F[i])

    listORFs=[]
    for i in range(len(newD)):
        listORFs.append({'start':abs(newD[i])+1,'stop':abs(newF[i]),'frame':'','length':newF[i]-newD[i],'protein':""})
        if newF[i]%3==0 and newF[i]>=0:
            listORFs[i]['frame']=1
        if newF[i]%3==1 and newF[i]>=0:
            listORFs[i]['frame']=2
        if newF[i]%3==2 and newF[i]>=0:
            listORFs[i]['frame']=3
        if newF[i]%3==0 and newF[i]<0:
            listORFs[i]['frame']=-1
        if newF[i]%3==1 and newF[i]<0:
            listORFs[i]['frame']=-2
        if newF[i]%3==2 and newF[i]<0:
            listORFs[i]['frame']=-3

        orf=nseq[newD[i]:newF[i]]

        listORFs[i]['protein']=translate(orf,4)
    return listORFs


def getLengths(orf_list):
    result_list=[]
    for i in range(len(orf_list)):
        result_list.append(orf_list[i]['length'])
    return result_list

def getLongestORF(orflist):
    a=getLengths(orflist)
    max=a[0]
    for i in a:
        if i>max:
            max=i
    return max

def getTopLongestORF(orflist,value):
    a=getLengths(orflist)
    seuil=a[round(-(len(a)*value))]
    toplist=[]
    for i in a:
        if i>seuil:
            toplist.append(i)
    return toplist
    
                           

seqseq={'data':'ATGTAGTAGTAGTAGTAGTATATGTAGTATGTGATGTATGATAGTAGTATGATGTAGTAGCAT'}
a=findORF(seqseq,5,1)
b=getLengths(a)
c=getLongestORF(a)
d=getTopLongestORF(a,0.1)

print(a,b,c,d)


sequence=bio.readFASTA('./seq.txt')

a=findORF(sequence,180,4)
b=getLengths(a)
c=getLongestORF(a)
d=getTopLongestORF(a,0.1)

print(a,b,c,d)
