import sys, re, time, urllib.request, urllib.parse, html
for q in sys.argv[1:]:
    url='https://arxiv.org/search/?query='+urllib.parse.quote(q)+'&searchtype=all&abstracts=hide&order=-announced_date_first&size=50'
    try:
        s=urllib.request.urlopen(url, timeout=60).read().decode()
    except Exception as e:
        print('ERR',q,e); continue
    print('=== ',q)
    for m in re.finditer(r'<p class="list-title[^"]*">\s*<a href="https://arxiv.org/abs/([^"]+)".*?<p class="title is-5 mathjax">(.*?)</p>\s*<p class="authors">(.*?)</p>',s,re.S):
        t=' '.join(html.unescape(re.sub('<[^>]+>','',m.group(2))).split())
        a=' '.join(html.unescape(re.sub('<[^>]+>','',m.group(3))).split())
        print(m.group(1),'|',t,'|',a[:120])
    time.sleep(2)
