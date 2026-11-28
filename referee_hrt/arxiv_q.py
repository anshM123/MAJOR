import sys, time, re, urllib.request, urllib.parse
def get(url, tries=6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                s = r.read().decode()
            if 'Rate exceeded' in s[:200]:
                raise Exception('rate')
            return s
        except Exception as e:
            time.sleep(5*(k+1))
    return ''
def show(s):
    for e in s.split('<entry>')[1:]:
        i=re.search(r'<id>(.*?)</id>',e).group(1)
        t=re.search(r'<title>(.*?)</title>',e,re.S)
        a=re.findall(r'<name>(.*?)</name>',e)
        p=re.search(r'<published>(.*?)</published>',e)
        ab=re.search(r'<summary>(.*?)</summary>',e,re.S)
        print('ID', i, '|', ' '.join(t.group(1).split()) if t else None, '|', ', '.join(a[:6]), '|', p.group(1)[:10] if p else '')
        if '-a' in sys.argv: print('   ', ' '.join(ab.group(1).split()) if ab else '')
mode=sys.argv[1]; arg=sys.argv[2]
if mode=='id':
    url='https://export.arxiv.org/api/query?id_list='+arg+'&max_results=50'
else:
    url='https://export.arxiv.org/api/query?search_query='+urllib.parse.quote(arg)+'&max_results=25&sortBy=submittedDate&sortOrder=descending'
show(get(url))
