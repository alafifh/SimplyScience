from flask import Flask, request, jsonify
from flask_cors import CORS

from AI_Func import (
    client,
    fetch_pubmed_abstracts,
    extract_claims,
    get_facts,
    CLAIM_DB,
)

app = Flask(__name__)
CORS(app, origins=["https://alafifh.github.io"])

@app.get("/health")
def health():
    return jsonify(ok=True)

@app.post("/search")
def search():
    data = request.get_json(silent=True) or {}
    query = (data.get("query") or "").strip()
    max_results = int(data.get("max_results") or 5)

    if not query:
        return jsonify(ok=False, error="No query provided"), 400

    pubmed_text = fetch_pubmed_abstracts(query, max_results=max_results)

    # ✅ Create Gemini chat session (THIS fixes your crash)
    chat = client.chats.create(model="models/gemini-flash-lite-latest")

    # ✅ Prevent mixing results between requests
    CLAIM_DB.clear()

    extract_claims(pubmed_text, query, chat)
    facts = get_facts()

    return jsonify(ok=True, facts=facts)
