# -*- coding: utf-8 -*-
"""
Created on Thu May 20 14:57:14 2021

@author: norab
"""

import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'UniProtKB ID',
'to': 'ENSEMBL_ID',
'format': 'tab',
'query': 'P40925 P40926 O43175 Q9UM73 P97793'
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
print(response.decode('utf-8'))