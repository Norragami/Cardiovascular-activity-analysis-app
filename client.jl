using HTTP

# Send a GET request to the server
response = HTTP.get("http://127.0.0.1:8080/data")
  # Convert response body to string
data = JSON.parse(String(response.body))



dataToSend = Dict("path" => "Мельникова_Елизавета_Дмитриевна2_21-04-22_13-02-11_")
# dataToSend = Dict("path" => "Мельникова_Елизавета_Дмитриевна_21-04-22_11-43-20_")
json_data = JSON.json(dataToSend)
response2 = HTTP.post("http://127.0.0.1:8080/getRrIntervals",headers = Dict("Content-Type" => "application/json"),body = json_data)
ownedData = JSON.parse(String(response2.body))

#Для тестирования работы запросов — нужное раскомментировать
outputRR = ownedData["RrIntervals"]
mean(outputRR)


outputECG = ownedData["outputECG"]
outputQ_x = ownedData["outputQ_x"]
outputR_x = ownedData["outputR_x"]
outputS_x = ownedData["outputS_x"]
outputQ_y = ownedData["outputQ_y"]
outputR_y = ownedData["outputR_y"]
outputS_y = ownedData["outputS_y"]

# outputPPG = ownedData["outputPPG"]
# outputPPGPeaks_x = ownedData["outputPPGPeaks_x"]
# outputPPGPeaks_y = ownedData["outputPPGPeaks_y"]
# outputPPGMins_x = ownedData["outputPPGMins_x"]
# outputPPGMins_y = ownedData["outputPPGMins_y"]


# outputAP = ownedData["outputAP"]
# outputAPPeaks_x = ownedData["outputAP_Peaks_x"]
# outputAPPeaks_y = ownedData["outputAP_Peaks_y"]
# outputAPMins_x = ownedData["outputAP_Mins_x"]
# outputAPMins_y = ownedData["outputAP_Mins_y"]

Xcoordinate = range(1,length=length(outputECG),step=1) # Для построения графика задать начальную точку

output = ownedData["outputDecimatedECG"]
output2 = ownedData["outputECG"]

plotly()
plot(outputRR,layout=(1,1),legend=false)
# scatter!(outputR_x, outputR_y)
# scatter!(outputQ_x, outputQ_y)
# scatter!(outputS_x, outputS_y)
# scatter!(outputPPGPeaks_x, outputPPGPeaks_y)
# scatter!(outputPPGMins_x, outputPPGMins_y)

scatter!(outputAPPeaks_x, outputAPPeaks_y)
scatter!(outputAPMins_x, outputAPMins_y)