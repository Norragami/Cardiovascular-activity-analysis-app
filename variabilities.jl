#Функции расчета вариабельности
function variabilityPeaks(Peaks_x_updt::Vector{Int64}) 
    P_P= fill(0.0, length(Peaks_x_updt))
    ind = 1
    n=2
    for i in collect(1:length(Peaks_x_updt))
        if i < length(Peaks_x_updt)
        
            P_P[ind]= abs(Peaks_x_updt[i+1]-Peaks_x_updt[i])
            ind+=1
        end
        if i == length(Peaks_x_updt)
            break
        end
    end
    M=mean(P_P)
    while n < length(P_P)
        if P_P[n]/P_P[n-1] > 10
            P_P[n]=P_P[n-1]
        end
        n+=1
    end
    P_P=P_P./fs
    bar(P_P,fillcolor=:blue,title="Ритмограмма интервалов между пиками ПВ",xlabel="N интервала",ylabel="Интервал между пиками, с",legend=false)
end
    
function variabilityReachTime(R_x_end::Vector{Int64},Mins_x_updt::Vector{Int64}) 
    ReachTime= fill(0,length(R_x_end))
    i=1
    n=1
   for p in collect(1:length(R_x_end)-1)
        if i>length(R_x_end) || n > length(R_x_end)
            break
        end
       if (R_x_end[i]-Mins_x_updt[n]) > 600 #Пропущен R
            n+=1
       end 
       if (R_x_end[i]-Mins_x_updt[n]) < -600  #Пропущен min ПВ
        i+=1  
       end
       ReachTime[p]=abs(R_x_end[i]-Mins_x_updt[n])
       
      i+=1
      n+=1
      
    end

    bar(ReachTime,fillcolor=:blue,title="Ритмограмма времени распространения ПВ",xlabel="N интервала",ylabel="Время распространения, мс",legend=false)
end

function variabilityR_R(R_x_end::Vector{Int64}) 
    R_R= fill(0.0, length(R_x_end)-1)
    ind = 1
    for i in collect(1:length(R_x_end))
        if i < length(R_x_end)
             R_R[ind]= abs(R_x_end[i]-R_x_end[i+1])
            ind+=1
         end
        if i == length(R_x_end)
             break
        end
    end
    R_R=R_R./fs
    mRR=mean(R_R) *1000 #Сразу перейдем к миллисекундам
    SDRR=std(R_R) * 1000
    R=collect(1:length(R_R))
    S=0
    k=0
    S1=0
    MSD=0
    rMSSD=0
    pNN50=0
    for i in R
        if i < length(R_R)
            S+= abs(R_R[i]-R_R[i+1])
            if abs(R_R[i]-R_R[i+1]) > 0.050
                 k+=1
            end
            S1+=(abs(R_R[i]-R_R[i+1]))^2
        end
        MSD=(S/(length(R_R)-1))
        rMSSD=sqrt(S1/(length(R_R)-1))
        pNN50=(k/(length(R_R)-1))
    end

    MSD=MSD*1000 # В миллисекундах
    rMSSD= rMSSD*1000 #В миллисекундах
    pNN50=pNN50*100 # В процентах
    return mRR,SDRR,MSD,rMSSD,pNN50
end

function variabilityAP(ap_Peaks_x_updt_end::Vector{Int64},ap_Mins_x_updt_end::Vector{Int64},ap_Peaks_y_updt_end::Vector{Float64},ap_Mins_y_updt_end::Vector{Float64}) 
    SAP = fill(0.0,length(ap_Mins_x_updt_end))
    DAP = fill(0.0,length(ap_Mins_x_updt_end))
    ind=1
    for i in collect(1:length(ap_Mins_x_updt_end)-1)
        if abs(ap_Mins_x_updt_end[i]-ap_Peaks_x_updt_end[i+1]) < 300
            if i < length(ap_Mins_x_updt_end)
                SAP[ind]=ap_Peaks_y_updt_end[i+1]
                DAP[ind]=ap_Mins_y_updt_end[i]
                ind+=1
            end
            if i == length(ap_Mins_x_updt_end)
                break
            end
        end
    end
    ind21=0
    Rem1=Vector{Int64}[]
    for i in collect(1:length(SAP))
        if SAP[i]==0
            insert!(Rem1,ind21+1,[i])
        end
    end

    Rem1=collect(Iterators.flatten(Rem1))
    for i in Rem1
        deleteat!(SAP,i)
        deleteat!(DAP,i)
    end
    PulseAP=fill(0.0,length(SAP))
    MeanAP=fill(0.0,length(SAP))
    for i in collect(1:length(SAP))
        PulseAP[i]=SAP[i]-DAP[i]
        MeanAP[i]=(1/3)*SAP[i]+(2/3)*DAP[i]
    end
    return SAP,DAP,PulseAP,MeanAP
end