import json
import requests

#s=requests.get('https://api.figshare.com/v2/articles/'+ str(8273102))
s=requests.get('https://api.figshare.com/v2/articles/'+ str(27921984))
r=json.loads(s.text)
#print(r)

for rf in r["files"]:
    #print(rf['name'],rf['download_url'])
    print("wget --content-disposition "+rf['download_url']+" -O data/"+rf['name'] )
