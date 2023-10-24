def get_access_token(client_id, client_secret, idt_username, idt_password):
    # Construct the HTTP request
    authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
    request_headers = { "Content-Type" : "application/x-www-form-urlencoded",
                        "Authorization" : "Basic " + authorization_string }
                    
    data_dict = {   "grant_type" : "password",
                    "scope" : "test",
                    "username" : idt_username,
                    "password" : idt_password }
    request_data = parse.urlencode(data_dict).encode()

    post_request = request.Request("https://www.idtdna.com/Identityserver/connect/token", 
                                    data = request_data, 
                                    headers = request_headers,
                                    method = "POST")

    # Transmit the HTTP request and get HTTP response
    response = request.urlopen(post_request)

    # Process the HTTP response for the desired data
    body = response.read().decode()
    
    # Error and return the response from the endpoint if there was a problem
    if (response.status != 200):
        raise RuntimeError("Request failed with error code:" + response.status + "\nBody:\n" + body)
    
    body_dict = json.loads(body)
    return body_dict["access_token"]

def get_data_from_IDT(seq, token):
    conn = http.client.HTTPSConnection("www.idtdna.com")

    payload = json.dumps({
        "Sequence": seq,
        "NaConc": 50,
        "MgConc": 3,
        "DNTPsConc": 0.8,
        "OligoConc": 0.25,
        "NucleotideType": "DNA"
    })

    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer ' + token
    }

    conn.request("POST", "/restapi/v1/OligoAnalyzer/Analyze", payload, headers)
    res = conn.getresponse()
    data = res.read()
    
    # Parse the JSON response
    response_data = json.loads(data.decode("utf-8"))
    
    # Print only the "MeltTemp" value
    melt_temp = response_data["MeltTemp"]
    return(melt_temp)   
    
def get_mismatch_from_IDT(seq, comp_seq, token):
    conn = http.client.HTTPSConnection("www.idtdna.com")

    payload = json.dumps({
  "Settings": {
    "Sequence": seq,
    "NaConc": 50,
        "MgConc": 3,
        "DNTPsConc": 0.8,
        "OligoConc": 0.25
  },
  "Sequence2": comp_seq
})

    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer ' + token
    }

    conn.request("POST", "/restapi/v1/OligoAnalyzer/TmMisMatch", payload, headers)
    res = conn.getresponse()
    data = res.read()
    
    # Parse the JSON response
    response_data = json.loads(data.decode("utf-8"))
    
    # Print only the "MeltTemp" value
    missmatch_tm = response_data["MeltTemp"]
    return(missmatch_tm)   
    
def get_hairpin_data_from_IDT(seq, token):
  
    conn = http.client.HTTPSConnection("www.idtdna.com")
    payload = json.dumps({
        "Sequence": seq,
        "NaConc": 50,
  "FoldingTemp": 25,
  "MgConc": 3,
  "NucleotideType": "DNA"
    })

    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer ' + token
    }

    conn.request("POST", "/restapi/v1/OligoAnalyzer/Hairpin", payload, headers)
    res = conn.getresponse()
    data = res.read()
    
    # Parse the JSON response
    response_data = json.loads(data.decode("utf-8"))
    
    # Print only the "deltaG" value
    delta_G = str(response_data[0]["deltaG"])
    return(delta_G)   
    

def get_selfdimer_data_from_IDT(seq, token):
    url = "https://www.idtdna.com/restapi/v1/OligoAnalyzer/SelfDimer"
    headers = {
        'Accept': 'application/json',
        'Authorization': f'Bearer {token}'
    }

    payload = {
        'primary': seq
    }

    response = requests.post(url, headers=headers, json=payload)

    if response.status_code == 200:
        response_data = response.json()
        # Check if the response is a list and not empty
        if isinstance(response_data, list) and response_data:
            # Extract the DeltaG from the first analysis result
            delta_G = response_data[0].get("DeltaG", "DeltaG value not found")
            return delta_G
        else:
            return "No self-dimer analysis results found."
    else:
        return f"Request failed with status code {response.status_code}"
