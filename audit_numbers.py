import re,sys
P=open('model_collapse_math_protocol.md').read()
def sec(start,end): 
    i=P.index(start); j=P.index(end,i+5) if end else len(P); return P[i:j]
S={'25.B':(sec('## 25.B','\n## 26. Block 1' if '\n## 26. Block 1' in P else None),['phase3_s25_analysis.txt']),
   '26.B':(sec('## 26.B Block 1 results','### 26.C'),['phase3_s26_analysis.txt']),
   '27.B':(sec('## 27.B Block 2 results','## 28.B'),['phase3_s27_analysis.txt']),
   '28.B':(sec('## 28.B Block 3 results','## 29.'),['phase3_s28_analysis.txt']),
   '30.B':(sec('## 30.B Confirmatory H1 results',None),['phase3_s30_analysis.txt'])}
num=re.compile(r'[-−+]?\d+\.\d+(?:e-?\d+)?')
for k,(txt,files) in S.items():
    src=' '.join(open(f).read() for f in files); vals=[float(x) for x in num.findall(src) if 'e' not in x]
    miss=[]; tot=0
    for m in num.findall(txt):
        if 'e' in m: continue
        s=m.replace('−','-'); d=len(s.split('.')[1]); x=abs(float(s)); tot+=1
        if not any(abs(round(abs(v),d)-x)<1e-9 for v in vals): miss.append(m)
    print(f'{k}: {tot} decimals in text, {len(miss)} not found verbatim (rounded) in analysis output; sample: {sorted(set(miss))[:25]}')
