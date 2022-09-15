struct Graph
    n :: Int # |V|
    m :: Int # |E|
    u :: Array{Int, 1}
    v :: Array{Int, 1} # uv is an edge
#    w :: Array{Float64, 1} # weight of each edge
    #nbr :: Array{Array{Int, 1}, 1}
end

include("edgecore.jl")

function get_graph(ffname)
    n = 0
    Label = Dict{Int32, Int32}()
    Origin = Dict{Int32, Int32}()
    #E = Set{Tuple{Int32, Int32, Float32}}()

    getID(x :: Int) = haskey(Label, x) ? Label[x] : Label[x] = n += 1

    fname = string("data/",ffname)
    fin = open(fname, "r")


    str = readline(fin)
    str = split(str)
    #n   = parse(Int, str[1])
    m   = parse(Int, str[3])
    u = Int[]
    v = Int[]
    #nbr=
    tot = 0
    for i = 1 : m
        str = readline(fin)
        str = split(str)
        x   = parse(Int, str[1])
        y   = parse(Int, str[2])
        if x!=y
            u1 = getID(x)
            v1 = getID(y)
            Origin[u1] = x
            Origin[v1] = y
            push!(u, u1)
            push!(v, v1)
            tot += 1
            #push!(nbr[u1],v1);
            #push!(nbr[v1],u1);
        end
    end
    close(fin)
    return Graph(n, tot, u, v)
end


function findconnect(G);
    n=G.n;
    bcj=zeros(G.n);
    idt=0;
    noc=zeros(G.n);
    L=lapsp(G);
    while minimum(bcj)==0
        idt=idt+1;
        label=0;
        for i=1:G.n
            if bcj[i]==0
                label=i;
                break;
            end
        end
        b=zeros(n);
        b[1]=label;
        bcj[label]=idt;
        noc[Int(idt)]=1;
        f=1;r=2;
        while f<r
            nbr=L[:,Int(b[Int(f)]):Int(b[Int(f)])].rowval;
            for i=1:size(nbr)[1]
                if bcj[Int(nbr[i])]==0
                    b[r]=nbr[i];
                    bcj[Int(nbr[i])]=idt;
                    noc[Int(idt)]+=1;
                    r=r+1;
                end
            end
            f=f+1;
        end
    end
    cop=argmax(noc);
    copset=union([]);


    nc = 0
    Label = Dict{Int32, Int32}()
    Origin = Dict{Int32, Int32}()
    #E = Set{Tuple{Int32, Int32, Float32}}()

    getID(x :: Int) = haskey(Label, x) ? Label[x] : Label[x] = nc += 1

    for i=1:n
        if bcj[i]==cop
            union!(copset,i);
        end
    end
    #nc=size(copset)[1];
    mc=0;
    uc = Int[]
    vc = Int[]
    for i=1:G.m
        if bcj[G.u[i]]==cop
            mc+=1;
            u1 = getID(G.u[i])
            v1 = getID(G.v[i])
            Origin[u1] =  G.u[i]
            Origin[v1] =  G.v[i]
            push!(uc,u1)
            push!(vc,v1)
        end
    end
    #u=uc;
    #v=vc;
    return Graph(nc,mc,uc,vc);
end
