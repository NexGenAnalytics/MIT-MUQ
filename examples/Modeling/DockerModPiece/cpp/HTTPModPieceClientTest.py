import requests

host = "https://demo.linusseelinger.de"

print("Requesting input sizes...")
r = requests.get(f"{host}/GetInputSizes")
print(r.text)

print("Requesting output sizes...")
r = requests.get(f"{host}/GetOutputSizes")
print(r.text)

print("Requesting evaluation")
r = requests.post(f"{host}/Evaluate", json={"level":0, "input0":[0.0] })
print(r.text)
