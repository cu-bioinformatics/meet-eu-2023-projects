import requests
import json
import os

# dictionary format: {cabs_dock_job_id : peptid_number}

# Peptid 9 is missing

jobs={  "95939f1c72fa88b":1 ,"bbf8c26bd46a8af":2 ,"51ef516a7e5ffd3":3 ,"c1f4cd4961b98ff":4 ,"adccf2fb63bc8cc":5,
	"b0a641eaa67c880":6 ,"5b7d5fb6fc2a2ff":7 ,"213789e7f16caa5":8 ,                     "76a2206191b292c":10,
	"622d6442766a3f3":11,"b0e6835bfeb73a3":12,"e3b564a858e219f":13,"b71d28b84ba6142":14,"8413fcd0417cc6f":15,
	"3b889f60a89a98a":16,"e89ed893d7f9c4":17,"ec32fffa955c314":18,"8136ed7a580caa6":19,"e523254b59de882":20,
	"51ba0430d8b9d3f":21,"2193616c5936a6b":22,"14795f4d7f87a2":23,"245adb08c682b60":24,"fd52b20059dd172":25,
	"93338a3db5067ef":26,"19551ee79b06c42":27}

# path to the directory to which we save downloaded models; it must already exist
path = "/media/mateusz/dysk_E/IPZ/wyniki_dokowania_baza/cabs_dock/"

# template of directories names
templ="6zslA_K"


check=[]
for job in list(jobs.keys()):
	url = f'https://biocomp.chem.uw.edu.pl/CABSdock/REST/status/{job}'
	dic=f"{path}{templ}{jobs[job]}"
	if not os.path.exists(dic):
		os.mkdir(dic)
		print(f"Creating directory: {dic}")
	response = requests.get(url).json()
	print(f"protein nr {jobs[job]}\n{response}")
	if response['status'] =='done':
		check.append(jobs[job])
		url2=f'https://biocomp.chem.uw.edu.pl/CABSdock/REST/job_results/{job}'
		response = requests.get(url2).json()
		for i,m in enumerate(response["models"]):
			f=open(f"{dic}/model_{i+1}.pdb","w")
			f.write(m["model_data"])
			f.close()

print(f"Successfuly downloaded peptids:\n{check.sort()}")
			
