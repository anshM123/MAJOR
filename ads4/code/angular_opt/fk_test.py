"""Is f_k = prod_{i<k} (N''-i)(J+1+i)/((N-i)(N+1+i)), N''=floor(N), a truncated Hausdorff moment sequence on [0,1] (k<=n)?
Check: Hankel H=[f_{i+j}] and [f_{i+j}-f_{i+j+1}] PSD (n=2m+1: [f_{i+j+1}], [f_{i+j}-f_{i+j+1}], i,j<=m)."""
import mpmath as mp, sys
mp.mp.dps = 400
def fseq(J, N, n, b=0):
    Nf = mp.floor(N); f=[mp.mpf(1)]
    # b-shifted version: J->J-b... use generic: prod (Nf-i)(J+1+i)/((N-i)(N+1+i)) with J,N given
    for k in range(n):
        f.append(f[-1]*(Nf-k)*(J+1+k)/((N-k)*(N+1+k)))
    return f
def hausdorff_min(f):
    n=len(f)-1
    if n%2==1:
        m=(n-1)//2
        H1=mp.matrix(m+1,m+1); H2=mp.matrix(m+1,m+1)
        for i in range(m+1):
            for j in range(m+1):
                H1[i,j]=f[i+j+1]; H2[i,j]=f[i+j]-f[i+j+1]
    else:
        m=n//2
        H1=mp.matrix(m+1,m+1); H2=mp.matrix(m,m)
        for i in range(m+1):
            for j in range(m+1): H1[i,j]=f[i+j]
        for i in range(m):
            for j in range(m): H2[i,j]=f[i+j+1]-f[i+j+2]
    # min eigen via cholesky attempt (relative)
    def mineig(H):
        E=mp.eigsy(H)[0]; return min(E)/max(E)
    return mineig(H1), mineig(H2)
for (J,h) in [(150,12.3),(300,20.7),(1000,30.5),(5000,50.5),(20000,100.7)]:
    N=mp.mpf(J+h-2)
    for n in [5,11,21]:
        f=fseq(J,N,n)
        a,b=hausdorff_min(f)
        print(J,h,n,mp.nstr(a,3),mp.nstr(b,3),flush=True)
