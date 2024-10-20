import json
import requests

s=requests.get('https://api.figshare.com/v2/articles/'+ str(8273102))
r=json.loads(s.text)
print(r)
