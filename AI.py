import os
import google.genai as genai
from google.genai import types
import time

from Bio import Entrez

key = os.getenv("GEMINI_API_KEY2")
client = genai.Client(api_key=key)

Entrez.email = "your_email@example.com"

def fetch_pubmed_abstracts(query, max_results=5):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]

    handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
    abstracts = handle.read()
    handle.close()
    return abstracts

query = "CRISPR gene editing"
pubmed_text = fetch_pubmed_abstracts(query)

chat = client.chats.create(model="gemini-2.5-pro")  # pick the current model

MAX_REQUESTS = 10
request_count = 0

if request_count < MAX_REQUESTS:
    chat = client.chats.create(model="gemini-2.5-pro")  # pick the current model
    request_count += 1
else:
    print("Request limit reached. No more requests can be made.")

MAX_REQUESTS = 10
request_count = 0

if request_count < MAX_REQUESTS:
    chat = client.chats.create(model="gemini-2.5-pro")  # pick the current model
    request_count += 1
    response1 = chat.send_message("Here are some PubMed abstracts on CRISPR gene editing:\n<your abstracts here>")
    time.sleep(30) 
else:
    print("Request limit reached. No more requests can be made.")
    time.sleep(30) 

'''
# Set your API key as an environment variable before running the code
key = os.getenv("GEMINI_API_KEY")
client =genai.Client()

response = client.models.generate_content(
    model="gemini-2.5-flash", contents="Explain how AI works in a few words"
)
print(response.text)'''