using HTTP

# Send a GET request to the server
response = HTTP.get("http://127.0.0.1:8080/data")
  # Convert response body to string
data = JSON.parse(String(response.body))



dataToSend = Dict("a" => 1, "b" => 2)
json_data = JSON.json(dataToSend)
response2 = HTTP.post("http://127.0.0.1:8080/add",headers = Dict("Content-Type" => "application/json"),body = json_data)
ownedData = JSON.parse(String(response2.body))
