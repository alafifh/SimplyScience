import os
import json
import google.genai as genai
from Bio import Entrez

Entrez.email = "alafifhams@gmail.com"  # NCBI requires an email

key = os.getenv("GEMINI_API_KEY2")
client = genai.Client(api_key=key)

chat = client.chats.create(model="models/gemini-flash-lite-latest")

CLAIM_DB = {}  # store structured claims

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

def extract_claims(pubmed_text, query, chat):
    """
    Send abstracts to Gemini and store (high/moderate) claims in CLAIM_DB.
    """
    prompt = f"""
Analyze the following PubMed abstracts related to {query}:
{pubmed_text}

Instructions:
- Only include studies directly related to {query}.
- Group articles by theme.
- Summarize findings using clear bullet points with several sentences (~150 words each).
- Each bullet must contain only ONE primary claim.
- Use language that the average high school student can understand.

Evidence rules:
- Use only peer-reviewed studies from reputable institutions.
- Prefer clinically meaningful sample sizes.
- Do not infer beyond the abstracts.

Output JSON ONLY, list of claims with fields:
- claim_id
- claim_text
- evidence_strength (High / Moderate / Preliminary)
- pmids (list of PMIDs used)

[
  {{
    "claim_id": "unique_id_here",
    "claim_text": "Main claim text from the abstract",
    "category": "Molecular / Clinical / Epidemiological",
    "evidence_strength": "High / Moderate / Preliminary",
    "pmids": ["list_of_PMIDs_used"]
  }}
]
Do not include ```json in output

"""
    response = chat.send_message(prompt)

    try:
        claims = json.loads(response.text)
    except json.JSONDecodeError:
        print("Error. Raw response:")
        print(response.text)
        return

    # Store only high-evidence claims
    for c in claims:
        if c.get("evidence_strength") == "High" or "Moderate":
            CLAIM_DB[c["claim_id"]] = c

def get_facts(category=None):
    """
    Return a list of claims for display. No sources here.    """
    output = []
    for claim in CLAIM_DB.values():
        if category and claim["category"] != category:
            continue
        output.append({
            "id": claim["claim_id"],
            "text": claim["claim_text"],
            "category": claim["category"]
        })
    return output

def get_sources(claim_id):
    claim = CLAIM_DB.get(claim_id)
    if not claim:
        return "Claim not found"
    sources_str = ", ".join(claim["pmids"])
    return f"Sources: {sources_str}"