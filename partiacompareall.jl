include("graph.jl")
include("edgecore.jl")
#include("Algorithm.jl")
using LinearAlgebra
using Arpack
using Laplacians
t1=time();

fname = open("filename.txt", "r")
str   = readline(fname);
gn     = parse(Int, str);

for gphs=1:gn

str = readline(fname);
str = split(str);
G   = get_graph(str[1]);


#######

#Gc=findconnect(G);
#G=Gc;
fout = open("partialcompareall.txt", "a");
println(fout);
println(fout,str[1]," n,m= ",G.n," ",G.m," 1.du 2.approx 3.exact 4.lmd 5.bet 6.clo ,k=",100);
close(fout);
###
n=G.n;
k=100;

### end

L=lapsp(G);
selc1=zeros(k);
selc2=zeros(k);
selc3=zeros(k);
selc4=zeros(k);
selc5=zeros(k);
selc6=zeros(k);


###########max
d=zeros(n);
for i=1:n
    d[i]=L[i,i];
end
x=argmax(d);
dmax=d[x];

D=lapsp(G);
for i=1:n
    D[i,i]=0;
end
D=-D;
l,v=eigs(D,nev=1);
if sum(v)>0
    cl=v[:,1];
else
    cl=-v[:,1];
end




cd=d;
cb=BetweennessCentrality(G);
cc=ClosenessCentrality(G);


for i=1:k
    x=argmax(cd);
    selc1[i]=x;
    cd[x]=-1000;
    x=argmax(cl);
    selc4[i]=x;
    cl[x]=-1000;
    x=argmax(cb);
    selc5[i]=x;
    cb[x]=-1000;
    x=argmax(cc);
    selc6[i]=x;
    cc[x]=-1000;
end


#############max end
ansp=zeros(k);
s=union(1:n);
setdiff!(s,Int(selc1[1]));
exist=ones(G.n);
orig=zeros(G.n);
exist[Int(selc1[1])]=0;
xiab=zeros(G.n);
t=1;
for i=1:n
    if exist[i]>0
        xiab[i]=t;
        orig[t]=i;
        t+=1;
    end
end
selc2[1]=selc1[1];
for i=2:k

    Lg=L[s,s];
    f=approxchol_sddm(Lg,tol=1e-4);
    u=randn(n-i+1);
    v=randn(n-i+1);
    delt=1;
    while delt>1e-2
        u=v;
        v=f(u);
        u=u/norm(u);
        v=v/norm(v);
        delt=abs(1/(u'*f(u))-1/(v'*f(v)))*(u'*f(u));
    end
    ansp[i-1]=1/(u'*f(u));


    vect=zeros(n);
    for t=1:G.m
        x=G.u[t];
        y=G.v[t];
        if x in s && y in s
            vect[x]+=u[Int(xiab[x])]*u[Int(xiab[y])];
            vect[y]+=u[Int(xiab[x])]*u[Int(xiab[y])];
        end
    end
    x=argmax(vect);
    y=x;
    setdiff!(s,y);
    exist[Int(y)]=0;
    selc2[i]=y;
    xiab=zeros(G.n);
    t=1;
    for tt=1:n
        if exist[tt]>0
            xiab[tt]=t;
            orig[t]=tt;
            t+=1;
        end
    end
end

selc3[1]=selc2[1];
s=union(1:n);
setdiff!(s,selc3[1]);
for i=2:k

    Lg=L[s,s];
    tmp=0;
    loc=0;
    for j=1:n
        if j in s
            selc3[i]=j;
            ttmp=sp_lmd_LG(L,selc3[1:i],n,i);
            if ttmp>tmp
                tmp=ttmp;
                loc=j;
            end
        end
    end
    selc3[i]=loc;
    setdiff!(s,loc);
end


fout = open("partialcompareall.txt", "a");
for i=1:k
    println(fout,i,' ',sp_lmd_LG(L,selc1[1:i],n,i),' ',sp_lmd_LG(L,selc2[1:i],n,i),' ',sp_lmd_LG(L,selc3[1:i],n,i),' ',sp_lmd_LG(L,selc4[1:i],n,i),' ',sp_lmd_LG(L,selc5[1:i],n,i),' ',sp_lmd_LG(L,selc6[1:i],n,i),' ');
end
close(fout);

end
close(fname)
