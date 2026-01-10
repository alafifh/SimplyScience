import os
import google.genai as genai
from google.genai import types
import time
import json
from Bio import Entrez

key = os.getenv("GEMINI_API_KEY2")
client = genai.Client(api_key=key)

Entrez.email = "your_email@example.com" #email for NCBI database access

CLAIM_DB = {} #database for claims

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
    Send abstracts to Gemini and store structured HIGH-EVIDENCE claims in CLAIM_DB
    """
    prompt = f"""
Analyze the following PubMed abstracts related to {query}:
{pubmed_text}

Instructions:
- Only include studies directly related to {query}.
- Group articles by theme.
- Summarize findings using clear bullet points (~200 words total).
- Each bullet must contain only ONE primary claim.

Evidence rules:
- Use only peer-reviewed studies from reputable institutions.
- Prefer clinically meaningful sample sizes.
- Only include HIGH-confidence findings.
- Do not infer beyond the abstracts.

For each bullet point:
- Label as [Molecular], [Clinical], or [Epidemiological]
- End with: Articles Referenced: N

Output JSON ONLY. Do NOT include backticks or extra text.
Return as a list of objects with fields:
- claim_id
- claim_text
- category
- evidence_strength (High / Moderate / Preliminary)
- pmids (list of PMIDs used)
"""

    response = chat.send_message(prompt)

    try:
        claims = json.loads(response.text)
        if isinstance(claims, str):
            claims = json.loads(claims)
        if not isinstance(claims, list):
            print("Unexpected Gemini format, expected a list of claims.")
            print(response.text)
            return
    except json.JSONDecodeError:
        print("Error parsing Gemini output. Raw response:")
        print(response.text)
        return

    for c in claims:
        if not isinstance(c, dict):
            continue
        if c.get("evidence_strength") == "High":
            CLAIM_DB[c["claim_id"]] = c



def get_facts(category=None):
    output = []
    for claim in CLAIM_DB.values():
        if category and claim["category"] != category:
            continue
        if claim.get("evidence_strength") != "High":
            continue
        output.append({
            "id": claim["claim_id"],
            "text": claim["claim_text"],
            "category": claim["category"],
            "evidence": claim.get("evidence_strength", "Unknown")  # <-- add this line
        })
    return output

def get_sources(claim_id):
    """
    Return PMIDs for a specific claim
    """
    claim = CLAIM_DB.get(claim_id)
    if not claim:
        return {"error": "Claim not found"}
    return {
        "claim_text": claim["claim_text"],
        "pmids": claim["pmids"]
    }