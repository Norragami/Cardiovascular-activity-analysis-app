using HTTP
using JSON
include("Ecg_analysis.jl")
# Define the server's host and port
const HOST = "127.0.0.1"  # Localhost
const PORT = 8080         # Port for the server

# Function to handle incoming HTTP requests
function request_handler(req::HTTP.Request)
    method = String(req.method)
    path = String(req.target)

    println("Request received: Method = ", method, ", Path = ", path)  # Debugging log

    if method == "GET" && path == "/"
        println("Serving '/' endpoint")
        return HTTP.Response(200, "Hello from the Julia server! Use /data for JSON response.")

    elseif method == "GET" && path == "/data"
        println("Serving '/data' endpoint")  # Debugging log
        data = Dict("message" => "Welcome!", "status" => "success", "code" => 200)
        json_response = JSON.json(data)  # Convert Dict to JSON
        println("JSON response: ", json_response)  # Verify JSON content
        return HTTP.Response(200, json_response)

    elseif method == "POST" && path == "/getECG"
        # body = String(req.body)
        # println("Serving '/add' endpoint. Received body: ", body)
        result = formHttpResponse(req)
        # data = Dict("received_data" => body, "message" => "Data received successfully!")
        
        return HTTP.Response(200, result)

    else
        println("Serving 404 for path: ", path)
        return HTTP.Response(404, "Not Found: Check your endpoint!")
    end
end

# Start the HTTP server
println("Starting Julia HTTP server at http://$HOST:$PORT...")
server_task = @async HTTP.serve(request_handler, HOST, PORT)
 sleep(5)
schedule(server_task, InterruptException())  # Gracefully stop the server
# println("Server stopped.")
