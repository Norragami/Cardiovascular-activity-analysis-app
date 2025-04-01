using HTTP

# Send a GET request to the server
response = HTTP.get("http://127.0.0.1:8080/data")
  # Convert response body to string
data = JSON.parse(String(response.body))



dataToSend = Dict("path" => "Мельникова_Елизавета_Дмитриевна_21-04-22_11-43-20_", "startPoint" => 100000, "endPoint" => 1000000)
json_data = JSON.json(dataToSend)
response2 = HTTP.post("http://127.0.0.1:8080/getPPG",headers = Dict("Content-Type" => "application/json"),body = json_data)
ownedData = JSON.parse(String(response2.body))

# outputECG = ownedData["outputECG"]
# outputQ_x = ownedData["outputQ_x"]
# outputR_x = ownedData["outputR_x"]
# outputS_x = ownedData["outputS_x"]
# outputQ_y = ownedData["outputQ_y"]
# outputR_y = ownedData["outputR_y"]
# outputS_y = ownedData["outputS_y"]

outputPPG = ownedData["outputPPG"]
outputPPGPeaks_x = ownedData["outputPPGPeaks_x"]
outputPPGPeaks_y = ownedData["outputPPGPeaks_y"]
outputPPGMins_x = ownedData["outputPPGMins_x"]
outputPPGMins_y = ownedData["outputPPGMins_y"]
Xcoordinate = range(100000,stop=1000000,length=length(outputPPG))
plotly()
 plot(Xcoordinate,outputPPG,layout=(1,1),legend=false)
#  scatter!(outputR_x, outputR_y)
#  scatter!(outputQ_x, outputQ_y)
#  scatter!(outputS_x, outputS_y)
scatter!(outputPPGPeaks_x, outputPPGPeaks_y)
scatter!(outputPPGMins_x, outputPPGMins_y)