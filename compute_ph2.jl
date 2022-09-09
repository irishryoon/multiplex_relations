using Distances
using Plots
using Ripserer
using Eirene
using PersistenceDiagrams
using Images
using LinearAlgebra
using JLD
using Random
include("extension_method.jl")
function getpath(id,time,celltype)
    if celltype == "Macrophage"
        return string("macrophage\\macrophage_id_",id,"_time_",time,".csv")
    end
    if celltype == "Tumour"
        return string("tumour\\tumour_id_",id,"_time_",time,".csv")
    end
    if celltype == "proTumourMacrophage"
        return string("protumourmacrophage\\protumourmacrophages_id_",id,"_time_",time,".csv")
    end
    if celltype == "antiTumourMacrophage"
        return string("antitumourmacrophage\\antitumourmacrophages_id_",id,"_time_",time,".csv")

    end
    if celltype == "Vessel"
        return string("vessel\\vessel_id_",id,"_time_",time,".csv")
    end
    if celltype == "Necrotic"
        return string("necrotic\\necrotic_id_",id,"_time_",time,".csv")
    end
end
function getNumberOf(id,time,celltype)
    file = getpath(id,time,celltype);
    if isfile(file)
        a = ezread(file);  
        numof = length(a[2:end-1,1]);
        return numof;
    else
        return 0;
    end
end
function getpos(id,time,celltype)
    file = getpath(id,time,celltype);
    if !isfile(file)
        return
    end
    a = ezread(file);                    # load csv file 
    b = a[2:end-1,1:2];                  # extract first 2 cols -> points_x | points_y
    b = convert(Array{Float64,2},b');

    ##### uncomment to randomly sample tumour cells #####
    #=
    if celltype == "Tumour"
        Random.seed!(3);
        ll = b;
        for i in 1:ceil(Int,0.8*size(b,2))
            row = rand((1:size(ll,2))) 
            #println("del row",row)
            #println(ll)
            ll = ll[:, setdiff(1:end, row)]
        end
        b = ll;
        #println("old length:",size(b,2))
        #println("new length:",size(ll,2))
    end
    =#
    ######
    return reinterpret(reshape,Tuple{Float64,Float64},b);  # correct form for ripser
end
function getdiagram(id,time,celltype)
    pointcloud = getpos(id,time,celltype);
    return ripserer(pointcloud);
end
function getimageform(sims,celltype,sz,var)
    all_diagrams = [last.([getdiagram(sim[1],sim[2],celltype) for sim in sims])];
    vec = reduce(vcat,all_diagrams);
    return PersistenceImage(vec;size=sz,sigma=var);
end
function getimageform0d(sims,celltype,sz,var)
    all_diagrams = [first.([getdiagram(sim[1],sim[2],celltype) for sim in sims])];
    vec = reduce(vcat,all_diagrams);
    return PersistenceImage(vec;size=[sz,sz],sigma=var);
end
function getimage(id,time,celltype,sz,var,sims)
    diagram = getdiagram(id,time,celltype)[2];
    image = getimageform(sims,celltype,sz,var);
    return image(diagram)
end
function getdist(id1,time1,id2,time2,celltype,sz,var,sims)
    image1 = getimage(id1,time1,celltype,sz,var,sims);
    image2 = getimage(id2,time2,celltype,sz,var,sims);
    return norm(image1 - image2);
end
function plotProportions(celltype1,celltype2,id)
    times = [250,300,350,400,450,500]
    n1 = [getNumberOf(id,time,celltype1) for time in times]
    n2 = [getNumberOf(id,time,celltype2) for time in times]
    plot(times,n1./(n1+n2),label=string("N1: ",celltype1),lw=3)
    plot!(times,n2./(n1+n2),label=string("N2: ",celltype2),lw=3)
    xlabel!("Time")
    ylabel!("Proportion")
end

function cleanSims(sims,celltype)
    todel = [];
    for sim in sims
        if !isfile(getpath(sim[1],sim[2],celltype))
            append!(todel,[sim])
        end
    end
    deleteat!(sims, findall(x->x in todel,sims));
    ##

    ## delete sims that have 2 or less cells 
    positions = [getpos(sim[1],sim[2],celltype) for sim in sims]
    todel2=[]
    for i in 1:length(sims)
        if (length(positions[i]) <= 2)
            sim = sims[i]
            append!(todel2,[sim])
        end
    end
    deleteat!(sims, findall(x->x in todel2,sims));
    N = length(sims);
    return sims
end
function getBinaryLabels(sims)
    labels = []
    for sim in sims
        id = sim[1]
        time = sim[2]
        M2 = getNumberOf(id,time,"proTumourMacrophage");
        M1plusM2 = getNumberOf(id,time,"Macrophage");
        if (M2 == nothing || M1plusM2 == nothing)
            append!(labels,0)
        elseif (M2 / M1plusM2 > 0.5)
            append!(labels,1)
        else
            append!(labels,0)
        end
    end
    return labels
end
function getConcentrationLabels(sims)
    labels = []
    for sim in sims
        id = sim[1]
        time = sim[2]
        M2 = getNumberOf(id,time,"proTumourMacrophage");
        M1plusM2 = getNumberOf(id,time,"Macrophage");
        ratio = M2 / M1plusM2
        append!(labels, ratio)
    end
    return labels
end


function getTumourEliminationLabels(sims)
    labels = []
    N = length(sims)
    for i in 1:N
        ID = sims[i][1]
        TIME = sims[i][2]
        final_time_that_exists = 0
        early_time_that_exists = 501
        for t in [500,450,400,350,300,250]
            if (t > final_time_that_exists && isfile(getpath(ID,t,celltype)))
                final_time_that_exists = t
            end
        end
        for t in [250,300,350,400,450]
            if (t < early_time_that_exists && isfile(getpath(ID,t,celltype)))
                early_time_that_exists = t
            end
        end
        #println("Early:",early_time_that_exists)
        #println("Final:",final_time_that_exists)
        # now final_time_that_exists holds final time available
        final_tumour_cell_count = getNumberOf(ID,final_time_that_exists,"Tumour")
        early_tumour_cell_count = getNumberOf(ID,early_time_that_exists,"Tumour")

        #threshold = 0 # meaning we check if final tumour size is less than 1-threshold * current
        # check if final / current < 1

        ratio = final_tumour_cell_count / early_tumour_cell_count
        append!(labels,ratio)
    end
    return labels
end

function getRipsFeatures(sims, celltype, size, variance)
    imform = getimageform(sims,celltype,size,variance);
    vector_features = []
    N = length(sims)
    for i in 1:N
        #println(i)
        sim = sims[i]
        append!(vector_features,[vec(imform(getdiagram(sim[1],sim[2],celltype)[2]))])
    end
    full_matrix = mapreduce(permutedims, vcat, vector_features)
    return full_matrix
end
function getRipsFeaturesH0(sims, celltype, size, variance)
    imform = getimageform0d(sims,celltype,size,variance);
    vector_features = []
    N = length(sims)
    for i in 1:N
        #println(i)
        sim = sims[i]
        append!(vector_features,[vec(imform(getdiagram(sim[1],sim[2],celltype)[2]))])
    end
    full_matrix = mapreduce(permutedims, vcat, vector_features)
    return full_matrix
end
function cleanSimsDowker(sims,celltype1,celltype2)
    for celltype in [celltype1,celltype2]
        todel = [];
        for sim in sims
            if !isfile(getpath(sim[1],sim[2],celltype))
                append!(todel,[sim])
            end
        end
        deleteat!(sims, findall(x->x in todel,sims));
        ##

        ## delete sims that have 2 or less cells - may have to change for H0
        positions = [getpos(sim[1],sim[2],celltype) for sim in sims]
        todel2=[]
        for i in 1:length(sims)
            if (length(positions[i]) <= 2)
                sim = sims[i]
                append!(todel2,[sim])
            end
        end
        deleteat!(sims, findall(x->x in todel2,sims));
    end
    return sims
end
function getDowkerfeatures(dimension, sims, celltype1, celltype2, sz, var)
    N = length(sims)
    W_diags_all = []
    for i in 1:N
        sim = sims[i]    
        P1_ = getpos(sim[1],sim[2],celltype1); 
        P2_ = getpos(sim[1],sim[2],celltype2);
        P1 = zeros(Float64,length(P1_),2);
        P2 = zeros(Float64,length(P2_),2);
        for i in 1:length(P1_)
            P1[i,1] = P1_[i][1]
            P1[i,2] = P1_[i][2]
        end
        for i in 1:length(P2_)
            P2[i,1] = P2_[i][1]
            P2[i,2] = P2_[i][2]
        end
        #P1S[i] = P1
        #P2S[i] = P2
        if length(P1) < length(P2)
            D_P1_P2 =  Images.ImageDistances.pairwise(Euclidean(), P1, P2, dims = 1)
        else
            D_P1_P2 =  Images.ImageDistances.pairwise(Euclidean(), P2, P1, dims = 1)
        end
        #println("starting witness persistence, i= ",i," out of: ",N)
        W_P1_P2 = ext.compute_Witness_persistence(D_P1_P2, maxdim = 1)
        
        W_barcode = Eirene.barcode(W_P1_P2["eirene_output"], dim = dimension);
        n = size(W_barcode,1)
        W_diag = PersistenceDiagram([(W_barcode[i,1], W_barcode[i,2]) for i = 1:n])
        append!(W_diags_all,[W_diag])
    end
    PI = PersistenceDiagrams.PersistenceImage(W_diags_all,
    sigma = var,
    size = sz)
    vector_features = []
    for i in 1:N
        #println(i," out of ",N)
        d = W_diags_all[i]
        append!(vector_features,[vec(PI(d))])
    end
    full_matrix = mapreduce(permutedims, vcat, vector_features)
    return full_matrix
end
