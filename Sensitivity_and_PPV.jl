#Функция для расчета чувствительности и положительной предсказательной доли

function Sense_PPV(Test::Vector{Int64},Ref::Vector{Any},d::Float64)    #Подаем координаты x найденных и референтных пиков/минимумов и допуск d (с)
    d=d*fs
    Refbeg=Ref.-d
    Refend=Ref.+d
    Testbeg=Test.-d
    Testend=Test.+d
    TP=0
    FP=0
    k=0
    for n in Test
    for i in collect(1:length(Refbeg))
        range78= collect(Refbeg[i]:1:Refend[i])
        if n in range78
            TP+=1
            k+=1
        end
    end
        if k==0
            FP+=1
        end
        k=0
    end
    
    v=0
    FN=0
    for m in Ref
        for i in collect(1:length(Testbeg))
            range79= collect(Testbeg[i]:1:Testend[i])
            if m in range79
                v+=1
            end
        end
            if v==0
                FN+=1
            end
            v=0
        end
        
    
    Sensitivity=TP/(TP+FN)*100
    PPV= TP/(TP+FP)*100
    
    return Sensitivity,PPV
end
    