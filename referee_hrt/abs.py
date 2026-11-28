import sys, re, time, urllib.request, html
for i in sys.argv[1:]:
    try:
        s=urllib.request.urlopen('https://arxiv.org/abs/'+i, timeout=60).read().decode()
    except Exception as e:
        print('ID',i,'ERR',e); continue
    t=re.search(r'<title>(.*?)</title>',s,re.S).group(1)
    a=re.findall(r'<meta name="citation_author" content="([^"]*)"',s)
    d=re.search(r'<meta name="citation_date" content="([^"]*)"',s)
    ab=re.search(r'<meta name="citation_abstract" content="([^"]*)"',s,re.S)
    print('ID',i,'|',html.unescape(t),'|','; '.join(a[:6]),'|',d.group(1) if d else '')
    if ab: print('  ABS:',' '.join(html.unescape(ab.group(1)).split()))
    print()
    time.sleep(1.5)
