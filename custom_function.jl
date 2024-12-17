using JSON
function addAB(req::HTTP.Request)
    
    data = JSON.parse(String(req.body))
    a = data["a"]
    b = data["b"]
    return a+b
     
end