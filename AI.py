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
    if not ids:
        return "No articles found."

    handle = Entrez.efetch(
        db="pubmed",
        id=",".join(ids),
        rettype="abstract",
        retmode="text"
    )
    abstracts = handle.read()
    handle.close()

    return abstracts

# Fetch PubMed data
query = "Brain Cancer Research 2026"
pubmed_text = fetch_pubmed_abstracts(query)

# Create ONE chat (use free-tier model)
chat = client.chats.create(
    model="models/gemini-flash-lite-latest"
)

# Send PubMed abstracts
response = chat.send_message(
    f"These are PubMed abstracts about {query}:\n\n{pubmed_text}\n\n"
    "Summarize the key findings in simple langauge." 
    "Make it around 200 words."
)

print(response.text)