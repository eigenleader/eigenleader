using LinearAlgebra
using SparseArrays
using Laplacians
using Random
using Arpack
##############chenche
#=
function betweenness()
	n=G.n;
	m=G.m;
	nbr=Array{Array{Int,1},1}(undef,n);
	for i=1:G.n
		nbr[i]=[];
	end
	for i=1:m
		x=G.u[i];
		y=G.v[i];
		push!(nbr[x],y);
		push!(nbr[y],x);
	end

	C_bet=zeros(n);

	for i=1:G.n
		pred=[];
		dist=zeros(G.n);
		sigma=zeros(G.n);
		dist[i]=0;
		sigma[i]=1;
		Q=union([]);
		push!(Q,i);
		S=union([]);
		while !empty(Q)
			v=pop(Q);
			push!(S,v);
			for w in neighbour(v)
				if dist[w]==0
					dist[w]=dist[v]+1
=#



function BetweennessCentrality(G)
    gg = zeros(Int, G.n)
    foreach(i -> gg[i] = i, 1 : G.n)
    g = Array{Array{Int32, 1}, 1}(undef, G.n)
    foreach(i -> g[i] = [], 1 : G.n)
    for i=1:G.m
        u=G.u[i];
        v=G.v[i];
        push!(g[u], v)
        push!(g[v], u)
    end
    C = zeros(G.n)
    p = Array{Array{Int32, 1}, 1}(undef, G.n)
    d = zeros(Int32, G.n)
    S = zeros(Int32, G.n+10)
    sigma = zeros(G.n)
    Q = zeros(Int32, G.n+10)
    delta = zeros(G.n)
    for s = 1 : G.n
        foreach(i -> p[i] = [], 1 : G.n)
        top = 0
        sigma .= 0
        sigma[s] = 1.0
        d .= -1
        d[s] = 0
        front = 1
        rear = 1
        Q[1] = s

        while front <= rear
            v = Q[front]
            front += 1
            top += 1
            S[top] = v
            for w in g[v]
                if d[w] < 0
                    rear += 1
                    Q[rear] = w
                    d[w] = d[v] + 1
                end
                if d[w] == (d[v] + 1)
                    sigma[w] += sigma[v]
                    push!(p[w], v)
                end
            end
        end

        delta .= 0

        while top > 0
            w = S[top]
            top -= 1
            for v in p[w]
                delta[v] += ((sigma[v] / sigma[w]) * (1 + delta[w]))
                if w != s
                    C[w] += delta[w]
                end
            end
        end

    end

    return C
end


function ClosenessCentrality(G)
    gg = zeros(Int, G.n)
    foreach(i -> gg[i] = i, 1 : G.n)
    g = Array{Array{Int32, 1}, 1}(undef, G.n)
    foreach(i -> g[i] = [], 1 : G.n)
    for i=1:G.m
        u=G.u[i];
        v=G.v[i];
        push!(g[u], v)
        push!(g[v], u)
    end
    C = zeros(G.n)
    d = zeros(Int32, G.n)
    Q = zeros(Int32, G.n+10)
    for s = 1 : G.n
        d .= -1
        d[s] = 0
        front = 1
        rear = 1
        Q[1] = s

        while front <= rear
            v = Q[front]
            front += 1
            for w in g[v]
                if d[w] < 0
                    rear += 1
                    Q[rear] = w
                    d[w] = d[v] + 1
                end
            end
        end

        C[s] = sum(d)
    end

    foreach(i -> C[i] = 1.0 / C[i], 1 : G.n)

    return C
end



function qpl(k,n,L);
	s=zeros(k);
	for i=1:k
		s[i]=i;
	end
	tmp=sp_lmd_LG(L,s,n,k);
	indx=k;
	while s[1]<=n-k+1
		tans=sp_lmd_LG(L,s,n,k);
		if tans>tmp
			tmp=tans;
		end
		s[k]+=1;
		if s[k]>n
			s[k]-=1;
			while (indx>=1) && (s[indx]==n+indx-k)
				indx-=1;
				if indx==0
					indx=1;
					break;
				end
			end
			s[indx]+=1;
			for i=indx+1:k
				s[i]=s[i-1]+1;
			end
			indx=k;
		end
	end
	return tmp;
end

function lmd_LG(L,ss,n)
	S=union(1:n);
	setdiff!(S,ss);
	LG=L[S,S];
	return eigmin(LG);
end

function sp_lmd_LG_norm(L,ss,n,t,dmax)
	S=union(1:n);
	setdiff!(S,ss);
	Lg=L[S,S];
	I=spzeros(n-t,n-t);
	for j=1:n-t
		I[j,j]=dmax;
	end
	ILg=I-Lg;
	lmd,u=eigs(ILg,nev=1);
	return dmax-lmd[1];
end

function sp_lmd_max(A,ss,n,t,d)
	S=union(1:n);
	I=spzeros(n-t,n-t);
	#setdiff!(S,ss);
	S=ss;
	AG=A[S,S];
	xx,dd=eigs(AG,nev=1);
	for i=1:n-t
		I[i,i]=xx[1];
	end
	#xx[1]=d;
	u=ones(n-t);
	v=ones(n-t);
	err=1e-3;
	delt=1;
	AG=I-AG;

#=
	while delt>1e-4
	#for i=1:100
		u=v;
		v=AG*u;
		u=u/norm(u);
		v=v/norm(v);
#		delt=abs((v'*AG*v-u'*AG*u)/(v'*AG*v-xx[1]));
		delt=norm(u-v);
	end
	=#
	l,r=eigs(AG,nev=1,tol=1e-4);
	return xx[1]-l[1];
end


function sp_lmd_LG(LL,ss,n,t)
	S=union(1:n);
	#setdiff!(S,ss);
	S=ss;
	Lgg=LL[S,S];
	err=1e-3;
	u=ones(n-t);
	v=ones(n-t);
	ff = approxchol_sddm(Lgg, tol=err);
	delt=1;
    while delt>1e-3
        u=v;
        v=ff(u);
        u=u/norm(u);
        v=v/norm(v);
        #delt=abs((u'*ff(u))/(u'*ff(u))-(u'*ff(u))/(v'*ff(v)));
		delt=norm(u-v);
    end
	return 1/(ff(v)'*v);
end



function lan(Lgg,n,t,err)
	u=ones(n-t);
	v=ones(n-t);
	ff = approxchol_sddm(Lgg, tol=err);
	delt=1;
    while delt>1e-3
        u=v;
        v=ff(u);
        u=u/norm(u);
        v=v/norm(v);
        delt=abs((u'*ff(u))/(u'*ff(u))-(u'*ff(u))/(v'*ff(v)));
		#delt=norm(u-v);
    end
	return 1/(ff(v)'*v);
end

function lap(G :: Graph)
    F = zeros(G.n, G.n);
    for i = 1 : G.m
        F[G.u[i], G.v[i]] -= 1
        F[G.v[i], G.u[i]] -= 1
        F[G.u[i], G.u[i]] += 1
        F[G.v[i], G.v[i]] += 1
    end
    return F
end

function lapsp(G :: Graph)
    F = spzeros(G.n, G.n)
    for i = 1 : G.m
        F[G.u[i], G.v[i]] -= 1
        F[G.v[i], G.u[i]] -= 1
        F[G.u[i], G.u[i]] += 1
        F[G.v[i], G.v[i]] += 1
    end
    return F
end

function eigcal(a::AbstractArray; nev = 1, tol = 0.0)
    f = approxchol_sddm(a)
    op = Laplacians.SqLinOp(true,1.0,size(a,1),f)
    e = eigs(op, which=:LM, nev=nev, tol=tol)
    e[1] .= 1 ./ e[1]
    return e
end
